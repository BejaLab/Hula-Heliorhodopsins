
suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(bioformatr)
library(tools)

input <- snakemake@input
output <- snakemake@output

synonyms_file <- unlist(input["synonyms_file"])
taxa_txt <- unlist(input["taxa_txt"])
bak_pfam_files <- unique(unlist(input["bak_pfam"]))
pfam_files <- unique(unlist(input["pfam"]))
presearch <- unlist(input["presearch"])
helio_groups <- unlist(input["groups"])
clstr <- paste0(unlist(input["cdhit"]), ".clstr")
gff_files <- unique(unlist(input["gff"]))
match_files <- unique(unlist(input["match"]))

labels_group_file <- unlist(output["labels_group"])
gff_tsv_file <- unlist(output["gff_tsv"])

filter_empty_dataframes <- function(x) dim(x)[1] > 0

pfam.cols <- c("gene", "aln.start", "aln.end", "env.start", "env.end", "hmm.acc", "hmm.name", "type", "hmm.start", "hmm.end", "hmm.length", "bit.score", "e.value", "significance", "clan")
blast.cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

taxa.synonyms <- read.table(synonyms_file, col.names = c("key", "value")) %>%
	with(setNames(value, key))

taxdata <- lapply(taxa_txt, readLines) %>%
	lapply(data.frame) %>%
	setNames(taxa_txt) %>%
	bind_rows(.id = "fname") %>%
	rename(ID = 2) %>%
	extract(fname, into = c("TYPE", "Taxon"), regex = "(.+)-(.+)\\.txt")

contigs <- read.table(presearch, col.names = c("ID", "hmm", "contig", "num_genes")) %>%
	distinct(ID, contig) %>%
	left_join(taxdata, by = "ID") %>%
	mutate(Taxon = recode(Taxon, !!!taxa.synonyms))

bkg.pfam <- lapply(bak_pfam_files, read.table, col.names = pfam.cols) %>%
	setNames(basename(bak_pfam_files)) %>%
	bind_rows(.id = "fname") %>%
	extract(fname, into = "ID", regex = "(.+)\\.pfam\\.txt") %>%
	mutate(clan = ifelse(clan == "No_clan", NA, clan)) %>%
	distinct(ID, gene, hmm.acc, clan)

bkg.labels <- select(bkg.pfam, ID, gene, hmm.acc) %>%
	group_by(ID, hmm.acc) %>%
	summarize(bkg.num = n(), .groups = "drop") %>%
	group_by(ID) %>%
	mutate(gene.num = n()) %>%
	ungroup

pfam <- lapply(pfam_files, read.table, col.names = pfam.cols) %>%
	Filter(filter_empty_dataframes, .) %>%
	bind_rows %>%
	mutate(clan = recode(clan, No_clan = NA_character_))

groups <- read.table(helio_groups, col.names = c("Representative", "Group"))
cdhit <- read.cdhit.clstr("phylogeny/targets_Heliorhodopsin_uniprot.cdhit.clstr") %>%
	left_join(groups, by = "Representative") %>%
	select(Group, target.gene = Seq.Name)

clusters <- read.cdhit.clstr(clstr) %>%
	select(cluster = Cluster, Seq.Len, gene = Seq.Name, gene_repr = Representative)

matches <- lapply(match_files, readLines) %>%
	lapply(data.frame) %>%
	setNames(file_path_sans_ext(basename(match_files))) %>%
	bind_rows(.id = "hmm") %>%
	rename(target.gene = 2)

gff <- lapply(gff_files, read.table, sep = "\t", col.names = c("contig", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")) %>%
	lapply(select, contig, start, end, strand, attributes) %>%
	Filter(filter_empty_dataframes, .) %>%
	bind_rows %>%
	extract(attributes, into = c("gene", "target.gene"), regex = "ID=([\\w_.-]+);Target=([\\w_.-]+)") %>%
	left_join(contigs,  by = "contig") %>%
	left_join(clusters, by = "gene") %>%
	left_join(cdhit,    by = "target.gene") %>%
	left_join(matches,  by = "target.gene") %>%
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
	{ write.table(., gff_tsv_file, sep = "\t", row.names = F); (.) } %>%
	group_by(cluster1, cluster2) %>%
	fill(cluster3, .direction = "downup") %>%
	group_by(cluster1, cluster3) %>%
	fill(cluster2, .direction = "downup") %>%
	ungroup %>%
	distinct(cluster1, cluster2, cluster3, to.target, .keep_all = T) %>%
	filter(target.gene %in% gene) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(target.num = n_distinct(target.gene, na.rm = T)) %>%
	ungroup

labels <- left_join(gff, pfam, by = "gene") %>%
	filter(to.target != 0, !is.na(hmm.acc)) %>%
	ungroup %>%
	distinct(gene, hmm.acc, .keep_all = T)

data <- select(labels, Taxon, hmm, Group, ID, target.num, target.gene, hmm.acc, to.target, to.target.abs, to.target.sign, to.total, co.express) %>%
	filter(Group != "" | hmm == "PR_XR")
# write.table(data, "output/all_data.tsv")

data.bkg <- left_join(bkg.labels, distinct(data, ID, hmm, Taxon, Group), by = "ID") %>%
	group_by(hmm, Taxon, Group, hmm.acc) %>%
	summarize(bkg.num = sum(bkg.num), gene.num = sum(gene.num), bkg.freq = bkg.num / gene.num, .groups = "drop")
labels.group <- arrange(data, to.target.abs) %>%
	group_by(hmm.acc) %>%
	filter(n_distinct(target.gene) > 1) %>%
	ungroup %>%
	complete(nesting(hmm, Group), nesting(ID, Taxon), hmm.acc) %>%
	left_join(bkg.labels, by = c("ID", "hmm.acc")) %>%
	complete(nesting(hmm, Group), Taxon, nesting(hmm.acc)) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(target.num = first(na.omit(target.num))) %>%
	mutate(has.bkg = !is.na(target.gene) & !is.na(bkg.num)) %>%
	group_by(hmm.acc, target.gene) %>%
	mutate(to.target.max = max(to.target), to.target.min = min(to.target), to.target.span = abs(to.target.max) + ifelse(to.target.max != to.target.min, abs(to.target.min), 0)) %>%
	group_by(hmm, Taxon, Group, hmm.acc) %>%
	summarize(
		target.num = first(target.num),
		target.num.label = n_distinct(target.gene, na.rm = T), label.num.all = sum(!is.na(target.gene)), co.express = sum(co.express, na.rm = T),
		to.target.abs.med = median(to.target.abs, na.rm = T), to.target.abs.avg = mean(to.target.abs, na.rm = T), to.target.abs.sd = sd(to.target.abs, na.rm = T),
		label.num = sum(has.bkg), label.gene.num = sum(to.target.span[!duplicated(target.gene) & has.bkg], na.rm = T),
		bkg.num = sum(bkg.num[!duplicated(ID)], na.rm = T), bkg.gene.num = sum(gene.num[!duplicated(ID)], na.rm = T),
		.groups = "drop"
	)
	#left_join(data.bkg, by = c("hmm", "Taxon", "Group", "label"))
write.table(labels.group, file = labels_group_file, row.names = F)

#labels.overlap <- mutate(labels, label = paste(hmm.acc, clan, sep = ":")) %>%
#	distinct(gene, label) %>%
#	group_by(gene) %>%
#	summarize(label1 = paste(label, collapse = ","), num = 1, .groups = "drop") %>%
#	group_by(label1) %>%
#	summarize(num = sum(num), label2 = label1, .groups = "drop") %>%
#	filter(num > 1) %>%
#	separate_rows(label1, sep = ",") %>%
#	separate_rows(label2, sep = ",") %>%
#	filter(label1 != label2) %>%
#	group_by(label1, label2) %>%
#	summarize(num = sum(num), .groups = "drop") %>%
#	separate(label1, into = c("hmm.acc1", "clan1"), sep = ":", fill = "left") %>%
#	separate(label2, into = c("hmm.acc2", "clan2"), sep = ":", fill = "left")

#write.table(labels.overlap, file = "output/labels.overlap.txt", row.names = F)
