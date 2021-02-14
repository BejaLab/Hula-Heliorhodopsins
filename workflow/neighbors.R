
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

taxfiles <- c("gbk","wgs") %>% paste0("-*.txt") %>% Sys.glob
taxdata <- lapply(taxfiles, readLines) %>%
	lapply(data.frame) %>%
	setNames(taxfiles) %>%
	bind_rows(.id = "fname") %>%
	rename(ID = 2) %>%
	extract(fname, into = c("TYPE", "Taxon"), regex = "(.+)-(.+)\\.txt")
bkg <- pipe("parallel --tagstring {/.} grep -c \"'^>'\" ::: background/*.faa") %>%
	read.table(col.names = c("ID","gene.num")) %>%
	left_join(taxdata, by = "ID") %>%
	group_by(Taxon) %>%
	mutate(taxon_gene.num = sum(gene.num))
pfamfiles <- Sys.glob("background/*.pfam")
bkg.pfam <- lapply(pfamfiles, read.table, col.names = pfam.cols) %>%
	setNames(pfamfiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = "ID", regex = "/(.+)\\.pfam") %>%
	mutate(clan = ifelse(clan == "No_clan", NA, clan)) %>%
	mutate(label = hmm.acc, label.type = "pfam") %>%
	distinct(ID, gene, label = hmm.acc, clan, label.type = "pfam")

bkg.clan <- filter(bkg.pfam, !is.na(clan)) %>%
	distinct(ID, gene, label = clan, label.type = "clan")

bkg.orthogroup <- pipe("parallel --tag sort -uk1,1 ::: background/*.blast") %>%
	read.table(sep = "\t", col.names = c("fname", blast.cols)) %>%
    extract(fname, into = "ID", regex = "/(.+)\\.blast") %>%
	left_join(proteinortho, by = c(sseqid = "gene_repr")) %>%
	left_join(proteinortho.min.scores, by = "orthogroup") %>%
	filter(evalue < 1e-10) %>%
	mutate(gene = qseqid, label = orthogroup, label.type = "orthogroup")

bkg.labels <- bind_rows(bkg.pfam, bkg.clan, bkg.orthogroup) %>%
	select(ID, gene, label, label.type) %>%
	group_by(ID, label, label.type) %>%
	summarize(bkg.num = n()) %>%
	group_by(ID) %>%
	mutate(gene.num = n()) %>%
	ungroup

contigs <- pipe("find hmmsearch-1 -name '*.match' -size +1c | parallel --tagstring {/.} cat") %>%
	read.table(header = F, col.names = c("ID", "contig", "desc"), sep = "\t", row.names = NULL, fill = T) %>%
	mutate(contig = sub("_[0-9]+$", "", contig)) %>%
	distinct(ID, contig) %>%
	extract(ID, into = "assembly", regex = "GC[AF]_([0-9]+)", remove = F) %>%
	extract(ID, into = "wgs", regex = "^([A-Z]+)[0-9.]+", remove = F) %>%
	left_join(taxdata, by = "ID") %>%
	mutate(id = ifelse(is.na(wgs), assembly, wgs)) %>%
	mutate(Taxon = recode(Taxon, !!!taxa.synonyms))

pfamfiles <- Sys.glob("targets/*/*.pfam")
pfam <- lapply(pfamfiles, read.table, col.names = pfam.cols) %>%
	setNames(pfamfiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = "target", regex = "/(.+)\\.pfam") %>%
	mutate(clan = recode(clan, No_clan = NA_character_))

groups <- read.table("targets/phylogeny/heliorhodopsin_groups.txt", col.names = c("Representative", "Group"))
cdhit <- read.cdhit.clstr("targets/phylogeny/targets_Heliorhodopsin_uniprot.cdhit.clstr") %>%
	left_join(groups, by = "Representative") %>%
	select(Group, target.gene = Seq.Name)

clusters <- read.cdhit.clstr("targets/targets.8.clstr") %>%
	select(cluster = Cluster, Seq.Len, gene = Seq.Name, gene_repr = Representative) %>%
	left_join(proteinortho, by = "gene_repr") %>%
	mutate(orthogroup = ifelse(is.na(orthogroup), gene_repr, orthogroup))

gfffiles <- Sys.glob("targets/*/*.gff")
gff <- lapply(gfffiles, read.table, sep = "\t", col.names = c("contig", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")) %>%
	lapply(select, contig, start, end, strand, attributes) %>%
	setNames(gfffiles) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = c("hmm", "target.gene"), regex = "targets/(.+)/(.+)\\.gff") %>%
	extract(attributes, into = "gene", regex = "ID=([A-Za-z0-9_.-]+)") %>%
	left_join(contigs,  by = "contig") %>%
	left_join(clusters, by = "gene") %>%
	left_join(cdhit,    by = "target.gene") %>%
	replace_na(list(Group = "")) %>%
	group_by(target.gene) %>%
	mutate(target.pos = which(gene == target.gene), target.strand = strand[target.pos],
		to.target = ifelse(target.strand == "-", 1, -1) * (target.pos - 1:n()), to.target.abs = abs(to.target), to.target.sign = sign(to.target)) %>%
	mutate(to.min = min(to.target), to.max = max(to.target), to.total = to.max - to.min) %>%
	filter(to.min < 0 | to.max > 0) %>%
	arrange(target.gene, to.target) %>%
	mutate(swap = strand != ifelse(to.target > 0, lag(strand), lead(strand))) %>%
	mutate(swap.num = sapply(min(to.target):max(to.target), function(x) which(to.target.abs <= abs(x) & to.target.sign == sign(x)) %>% `[`(swap,.) %>% sum)) %>%
	mutate(co.express = swap.num == 0 | swap.num == 1 & to.target < 0) %>%
	mutate(cluster1 = cluster[which(to.target == 0)], cluster2 = ifelse(to.min < 0, cluster[which(to.target == -1)], NA), cluster3 = ifelse(to.max > 0, cluster[which(to.target == 1)], NA)) %>%
	arrange(-to.total) %>%
	{ write.table(., "gff.tsv", sep = "\t", row.names = F); (.) } %>%
	group_by(cluster1, cluster2) %>%
	fill(cluster3, .direction = "downup") %>%
	group_by(cluster1, cluster3) %>%
	fill(cluster2, .direction = "downup") %>%
	ungroup %>%
	distinct(cluster1, cluster2, cluster3, to.target, .keep_all = T) %>%
	filter(target.gene %in% gene)

labels <- bind_rows(
	mutate(pfam, label = hmm.acc,    label.type = "pfam") %>% left_join(gff, ., by = "gene"),
	mutate(pfam, label = clan,       label.type = "clan") %>% filter(!is.na(clan)) %>% left_join(gff, ., by = "gene") %>% distinct(gene, label, .keep_all = T),
	mutate(gff,  label = orthogroup, label.type = "orthogroup")
) %>%
	filter(to.target != 0, !is.na(label)) %>%
	ungroup %>%
	distinct(gene, label, .keep_all = T)

data <- select(labels, Taxon, hmm, Group, ID, target.gene, label, label.type, to.target, to.target.abs, to.target.sign, to.total, co.express) %>%
	filter(Group != "" | hmm == "PR_XR", label.type == "pfam")

data.bkg <- left_join(bkg.labels, distinct(data, ID, hmm, Taxon, Group), by = "ID") %>%
	group_by(hmm, Taxon, Group, label, label.type) %>%
	summarize(bkg.num = sum(bkg.num), gene.num = sum(gene.num), bkg.freq = bkg.num / gene.num, .groups = "drop")
labels.group <- arrange(data, to.target.abs) %>%
	group_by(label, label.type) %>%
	filter(n_distinct(target.gene) > 1) %>%
	ungroup %>%
	complete(nesting(hmm, Group), nesting(ID, Taxon), nesting(label, label.type)) %>%
	left_join(bkg.labels, by = c("ID", "label", "label.type")) %>%
	complete(nesting(hmm, Group), Taxon, nesting(label, label.type)) %>%
	mutate(has.bkg = !is.na(target.gene) & !is.na(bkg.num)) %>%
	group_by(label, label.type, target.gene) %>%
	mutate(to.target.max = max(to.target), to.target.min = min(to.target), to.target.span = abs(to.target.max) + ifelse(to.target.max != to.target.min, abs(to.target.min), 0)) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(target.num = n_distinct(target.gene, na.rm = T)) %>%
	group_by(hmm, Taxon, Group, label, label.type) %>%
	summarize(
		target.num = first(target.num),
		target.num.label = n_distinct(target.gene, na.rm = T), label.num.all = sum(!is.na(target.gene)), co.express = sum(co.express, na.rm = T),
		to.target.abs.med = median(to.target.abs, na.rm = T), to.target.abs.avg = mean(to.target.abs, na.rm = T), to.target.abs.sd = sd(to.target.abs, na.rm = T),
		label.num = sum(has.bkg), label.gene.num = sum(to.target.span[!duplicated(target.gene) & has.bkg], na.rm = T),
		bkg.num = sum(bkg.num[!duplicated(ID)], na.rm = T), bkg.gene.num = sum(gene.num[!duplicated(ID)], na.rm = T),
		.groups = "drop"
	)
	#left_join(data.bkg, by = c("hmm", "Taxon", "Group", "label", "label.type"))
write.table(labels.group,   file = "labels.group.txt",   row.names = F)

labels.overlap <- mutate(labels, label = paste(label, label.type, clan, sep = ":")) %>%
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
	separate(label1, into = c("label1", "label.type_1", "clan1"), sep = ":", fill = "left") %>%
	separate(label2, into = c("label2", "label.type_2", "clan2"), sep = ":", fill = "left") %>%
	mutate(clan1 = ifelse(clan1 == "NA", NA, clan1), clan2 = ifelse(clan2 == "NA", NA, clan2))

write.table(labels.overlap, file = "labels.overlap.txt", row.names = F)
