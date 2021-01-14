
library(dplyr)
library(tidyr)
library(bioformatr)
options(stringsAsFactors = F)

pfam.cols <- c("gene", "aln.start", "aln.end", "env.start", "env.end", "hmm.acc", "hmm.name", "type", "hmm.start", "hmm.end", "hmm.length", "bit.score", "e.value", "significance", "clan")
blast.cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

taxa.synonyms <- list(
	Saccharibacteria = "Patescibacteria",
	Dojkabacteria    = "Patescibacteria",
	WWE3             = "Patescibacteria",

	Tenericutes      = "Firmicutes",
	Firmicutes_A     = "Firmicutes",
	Firmicutes_B     = "Firmicutes",
	Firmicutes_C     = "Firmicutes",
	Firmicutes_D     = "Firmicutes",
	Firmicutes_E     = "Firmicutes",
	Firmicutes_F     = "Firmicutes",
	Firmicutes_G     = "Firmicutes",

	Dictyoglomi      = "Thermocalda",
	Dictyoglomota    = "Thermocalda",
	Thermotogae      = "Thermocalda",
	Thermotogota     = "Thermocalda",
	Caldiserica      = "Thermocalda",
	Caldisericota    = "Thermocalda",

	Chloroflexota    = "Chloroflexi",

	Actinobacteriota = "Actinobacteria"
)

proteinortho <- read.table("targets/targets.8.proteinortho.tsv", sep = "\t", col.names = c("n_species", "n_genes", "alg_conn", "gene_repr")) %>%
	mutate(orthogroup = sprintf("g%05d", 1:n())) %>%
	separate_rows(gene_repr, sep = ",") %>%
	select(-n_species)
proteinortho.min.scores <- read.table("targets/targets.8.proteinortho-graph", sep = "\t", col.names = c("gene_repr", "b", "evalue_ab", "bitscore_ab", "evalue_ba", "bitscore_ba")) %>%
	left_join(proteinortho, by = "gene_repr") %>%
	mutate(bitscore = pmin(bitscore_ab, bitscore_ba)) %>%
	group_by(orthogroup) %>%
	summarize(min_bitscore = 0.8 * min(bitscore))

gtdbtk_taxonomy  <- Sys.glob("databases/*_taxonomy_r95.tsv") %>%
	lapply(read.table, sep = "\t", col.names = c("assembly", "TAXON")) %>%
	bind_rows %>%
	extract(1, into = "assembly", regex = "GC[AF]_([0-9]+)")
gtdbtk_data <- read.table("databases/assembly_summary_genbank.txt", sep = "\t", quote = "", fill = T, header = T, skip = 1, comment.char = "", na.strings = c("", "na")) %>%
	extract(1, into = "assembly", regex = "GC[AF]_([0-9]+)") %>%
	extract(wgs_master, into = "wgs", regex = "^([A-Z]+)[0-9.]+") %>%
	select(assembly, wgs) %>%
	left_join(gtdbtk_taxonomy, by = "assembly") %>%
	filter(!is.na(TAXON)) %>%
	extract(TAXON, into = c("d", "p"), regex = "d__([^;]+);p__([^;]+)") %>%
	{bind_rows(select(., id = 1, d, p), select(., id = 2, d, p))} %>%
	filter(!is.na(id)) %>%
	mutate(TAXON.gtdbtk = ifelse(d == "Archaea", d, p))

#gtdbtk <- Sys.getenv("GTDBTK_DATA_PATH") %>%
#	file.path("metadata/genome_metadata.tsv") %>%
#	read.table(sep = "\t", quote = "", header = T, fill = T, na.strings = c("", "none")) %>%
#	mutate(assembly = paste(accession, ncbi_genbank_assembly_accession, sep = ";") %>% gsub("GCF_", "GCA_", .)) %>%
#	extract(assembly, into = "ID.assembly", regex = "(GCA_[0-9.]+)", remove = F) %>%
#	extract(ncbi_wgs_master, into = c("ID.wgs.base", "ID.wgs.ver"), regex = "^([A-Z]+)[0-9]+[.]([0-9]+)", remove = F) %>%
#	mutate(ID.wgs = sprintf("%s%02d", ID.wgs.base, as.integer(ID.wgs.ver))) %>%
#	separate_rows(gtdb_taxonomy, sep = ";") %>%
#	mutate(TAXON = recode(gtdb_taxonomy, !!!taxa.synonyms.gtdbtk)) %>%
#	mutate(TAXON = ifelse(is.na(TAXON) & grepl("p__", gtdb_taxonomy), gtdb_taxonomy, TAXON))
#	group_by(accession) %>%
#	mutate(TAXON = TAXON[!is.na(TAXON)])

taxfiles <- Sys.glob("*-*.txt")
taxdata <- lapply(taxfiles, readLines) %>%
	lapply(data.frame) %>%
	setNames(taxfiles) %>%
	bind_rows(.id = "fname") %>%
	rename(ID = 2) %>%
	extract(fname, into = c("TYPE", "TAXON"), regex = "(.+)-(.+)\\.txt")
bkg <- pipe("parallel --tagstring {/.} grep -c \"'^>'\" ::: background/*.faa") %>%
	read.table(col.names = c("ID","gene.num")) %>%
	left_join(taxdata, by = "ID") %>%
	group_by(TAXON) %>%
	mutate(taxon_gene.num = sum(gene.num))
pfamfiles <- Sys.glob("background/*.pfam")
bkg.pfam <- lapply(pfamfiles, read.table, col.names = pfam.cols) %>%
	setNames(pfamfiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = "ID", regex = "/(.+)\\.pfam") %>%
	mutate(clan = ifelse(clan == "No_clan", NA, clan)) %>%
	mutate(label = hmm.acc, label_type = "pfam") %>%
	distinct(ID, gene, label = hmm.acc, clan, label_type = "pfam")

bkg.clan <- filter(bkg.pfam, !is.na(clan)) %>%
	distinct(ID, gene, label = clan, label_type = "clan")

bkg.orthogroup <- pipe("parallel --tag sort -uk1,1 ::: background/*.blast") %>%
	read.table(sep = "\t", col.names = c("fname", blast.cols)) %>%
    extract(fname, into = "ID", regex = "/(.+)\\.blast") %>%
	left_join(proteinortho, by = c(sseqid = "gene_repr")) %>%
	left_join(proteinortho.min.scores, by = "orthogroup") %>%
	filter(evalue < 1e-10) %>%
	mutate(gene = qseqid, label = orthogroup, label_type = "orthogroup")

bkg.labels <- bind_rows(bkg.pfam, bkg.clan, bkg.orthogroup) %>%
	select(ID, gene, label, label_type) %>%
	group_by(ID, label, label_type) %>%
	summarize(bkg.num = n()) %>%
	group_by(ID) %>%
	mutate(gene.num = n()) %>%
	ungroup

#bkg.labels.bak <- bind_rows(
#	bkg.pfam, bkg.clan, bkg.orthogroup
#) %>% select(ID, gene, label, label_type) %>%
#	left_join(bkg, by = "ID") %>%
#	group_by(TAXON, label, label_type) %>%
#	summarize(bkg.num = n(), taxon_gene.num = unique(taxon_gene.num), bkg.freq = bkg.num / taxon_gene.num)
	
contigs <- pipe("find hmmsearch-1 -name '*.match' -size +1c | parallel --tagstring {/.} cat") %>%
	read.table(header = F, col.names = c("ID", "contig", "desc"), sep = "\t", row.names = NULL, fill = T) %>%
	mutate(contig = sub("_[0-9]+$", "", contig)) %>%
	distinct(ID, contig) %>%
	extract(ID, into = "assembly", regex = "GC[AF]_([0-9]+)", remove = F) %>%
	extract(ID, into = "wgs", regex = "^([A-Z]+)[0-9.]+", remove = F) %>%
	left_join(taxdata, by = "ID") %>%
	mutate(id = ifelse(is.na(wgs), assembly, wgs)) %>%
	left_join(gtdbtk_data, by = "id") %>%
	mutate(TAXON = recode(TAXON, !!!taxa.synonyms), TAXON.gtdbtk = recode(TAXON.gtdbtk, !!!taxa.synonyms))

pfamfiles <- Sys.glob("targets/*/*.pfam")
pfam <- lapply(pfamfiles, read.table, col.names = pfam.cols) %>%
	setNames(pfamfiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = "target", regex = "/(.+)\\.pfam") %>%
	mutate(clan = ifelse(clan == "No_clan", NA, clan))

helio30 <- read.cdhit.clstr("targets/phylogeny/targets_Heliorhodopsin_uniprot.clust30.comb.clstr") %>%
	select(cluster30 = Cluster, gene.T = Seq.Name)

clusters <- read.cdhit.clstr("targets/targets.8.clstr") %>%
	select(cluster = Cluster, Seq.Len, gene = Seq.Name, gene_repr = Representative) %>%
	left_join(proteinortho, by = "gene_repr") %>%
	mutate(orthogroup = ifelse(is.na(orthogroup), gene_repr, orthogroup))

gfffiles <- Sys.glob("targets/*/*.gff")
gff <- lapply(gfffiles, read.table, sep = "\t", col.names = c("contig", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")) %>%
	lapply(select, contig, start, end, strand, attributes) %>%
	setNames(gfffiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = c("hmm", "gene.T"), regex = "targets/(.+)/(.+)\\.gff") %>%
	extract(attributes, into = "gene", regex = "ID=([A-Za-z0-9_.-]+)") %>%
	left_join(contigs, by = "contig") %>%
	left_join(clusters, by = "gene") %>%
	left_join(helio30, by = "gene.T") %>%
	group_by(gene.T) %>%
	mutate(T.pos = which(gene == gene.T), T.strand = strand[T.pos], to.T = ifelse(T.strand == "-", 1, -1) * (T.pos - 1:n()), to.T.abs = abs(to.T), to.T.sign = sign(to.T)) %>%
	mutate(to.min = min(to.T), to.max = max(to.T), to.total = to.max - to.min) %>%
	filter(to.min < 0 | to.max > 0) %>%
	arrange(gene.T, to.T) %>%
	mutate(swap = strand != ifelse(to.T > 0, lag(strand), lead(strand))) %>%
	mutate(swap.num = sapply(min(to.T):max(to.T), function(x) which(to.T.abs <= abs(x) & to.T.sign == sign(x)) %>% `[`(swap,.) %>% sum)) %>%
	mutate(co = swap.num == 0 | swap.num == 1 & to.T < 0) %>%
	mutate(cluster1 = cluster[which(to.T == 0)], cluster2 = ifelse(to.min < 0, cluster[which(to.T == -1)], NA), cluster3 = ifelse(to.max > 0, cluster[which(to.T == 1)], NA)) %>%
	arrange(-to.total) %>%
	{ write.table(., "gff.tsv", sep = "\t", row.names = F); (.) } %>%
	group_by(cluster1, cluster2) %>%
	fill(cluster3, .direction = "downup") %>%
	group_by(cluster1, cluster3) %>%
	fill(cluster2, .direction = "downup") %>%
	ungroup %>%
	distinct(cluster1, cluster2, cluster3, to.T, .keep_all = T) %>%
	filter(gene.T %in% gene)

labels <- bind_rows(
	mutate(pfam, label = hmm.acc,    label_type = "pfam") %>% left_join(gff, ., by = "gene"),
	mutate(pfam, label = clan,       label_type = "clan") %>% filter(!is.na(clan)) %>% left_join(gff, ., by = "gene") %>% distinct(gene, label, .keep_all = T),
	mutate(gff,  label = orthogroup, label_type = "orthogroup")
) %>% filter(to.T != 0, !is.na(label)) %>% ungroup

grouping.types <- c("TAXON", "TAXON.gtdbtk", "cluster30")
labels.group <- lapply(grouping.types, function(grouping.type) {
	data <- select(labels, grouping = !!grouping.type, hmm, ID, gene.T, label, label_type, to.T, to.T.abs, to.T.sign, to.total, co) %>%
		filter(!is.na(grouping))
	data.grouping <- distinct(data, ID, grouping)
	data.bkg <- left_join(bkg.labels, data.grouping, by = "ID") %>%
		filter(!is.na(grouping)) %>%
		group_by(grouping, label, label_type) %>%
		summarize(bkg.num = sum(bkg.num), gene.num = sum(gene.num), bkg.freq = bkg.num / gene.num, .groups = "drop")
	arrange(data, to.T.abs) %>%
		group_by(hmm, grouping) %>%
		mutate(T.num = n_distinct(gene.T), to.total.num = sum(to.total[!duplicated(gene.T)])) %>%
		group_by(hmm, grouping, label, label_type) %>%
		summarize(
			  T.num = first(T.num), to.total.num = first(to.total.num),
			  label.num.T = n_distinct(gene.T), label.num = n(), co.rate = sum(co)/n(), to.T.abs.med = median(to.T.abs), to.T.abs.avg = mean(to.T.abs), to.T.abs.sd = sd(to.T.abs),
			  .groups = "drop"
		) %>%
		left_join(data.bkg, by = c("grouping", "label", "label_type"))
}) %>% setNames(grouping.types) %>% bind_rows(.id = "grouping.type") %>%
	mutate(gene.freq = label.num / to.total.num, gene.T.freq = label.num.T / T.num)

labels.overlap <- mutate(labels, label = paste(label, label_type, clan, sep = ":")) %>%
	distinct(gene, label) %>%
	group_by(gene) %>%
	summarize(label1 = paste(label, collapse = ","), num = 1) %>%
	group_by(label1) %>%
	summarize(num = sum(num), label2 = label1) %>%
	filter(num > 1) %>%
	separate_rows(label1, sep = ",") %>%
	separate_rows(label2, sep = ",") %>%
	filter(label1 != label2) %>%
	group_by(label1, label2) %>%
	summarize(num = sum(num)) %>%
	separate(label1, into = c("label1", "label_type_1", "clan1"), sep = ":", fill = "left") %>%
	separate(label2, into = c("label2", "label_type_2", "clan2"), sep = ":", fill = "left") %>%
	mutate(clan1 = ifelse(clan1 == "NA", NA, clan1), clan2 = ifelse(clan2 == "NA", NA, clan2))

write.table(labels.group,   file = "labels.group.txt",   row.names = F)
write.table(labels.overlap, file = "labels.overlap.txt", row.names = F)
