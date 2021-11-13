
library(dplyr)
library(tidyr)
library(ggplot2)

input <- snakemake@input
output <- snakemake@output

gff_file <- unlist(input)
filtered_file <- unlist(output["filtered"])
all_file <- unlist(output["all"])

gff <- read.table(gff_file, header = T) %>%
	mutate(Group = ifelse(Group %in% c("A","B"), Group, "_")) %>%
	filter(target.gene %in% gene) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(targets_all = n_distinct(target.gene, na.rm = T), genomes_all = n_distinct(ID, na.rm = T)) %>%
	distinct(target.unique, .keep_all = T) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(targets = n_distinct(target.gene, na.rm = T), genomes = n_distinct(ID, na.rm = T)) %>%
	group_by(hmm, Taxon, Group) %>%
	summarize(targets_all = first(targets_all), genomes_all = first(genomes_all), targets = n_distinct(target.gene, na.rm = T), genomes = n_distinct(ID, na.rm = T)) %>%
	mutate(label = ifelse(targets != genomes, sprintf("%s\n(%s)", targets, genomes), as.character(targets))) %>%
	mutate(label_all = ifelse(targets_all != genomes_all, sprintf("%s\n(%s)", targets_all, genomes_all), as.character(targets_all))) %>%
	group_by(Taxon) %>%
	mutate(targets_total = sum(targets), targets_total_all = sum(targets_all), genomes_total = sum(genomes), genomes_total_all = sum(genomes_all)) %>%
	mutate(targets_pct = 100 * targets / targets_total, targets_pct_all = 100 * targets_all / targets_total_all) %>%
	mutate(n_pos = 100 - cumsum(targets_pct) + targets_pct/2, n_pos_all = 100 - cumsum(targets_pct_all) + targets_pct_all/2) %>%
	mutate(total_label = ifelse(targets_total != genomes_total, sprintf("%s\n(%s)", targets_total, genomes_total), as.character(targets_total))) %>%
	mutate(total_label_all = ifelse(targets_total_all != genomes_total_all, sprintf("%s\n(%s)", targets_total_all, genomes_total_all), as.character(targets_total_all))) %>%
	mutate(total_label = ifelse(row_number() == 1, total_label, ""), total_label_all = ifelse(row_number() == 1, total_label_all, "")) 

p1 <- ggplot(gff, aes(fill = paste(hmm, Group), y = targets_pct, x = 2)) +
	geom_bar(stat = "identity") +
	coord_polar("y", start = 0) +
	geom_text(aes(label = label, y = n_pos), size = 3) +
	geom_text(aes(label = total_label), x = 0.5, y = 0, size = 4) +
	xlim(c(0.5, 2.5)) +
	facet_grid(Taxon ~ .) +
	theme_void()
p2 <- ggplot(gff, aes(fill = paste(hmm, Group), y = targets_pct_all, x = 2)) +
	geom_bar(stat = "identity") +
	coord_polar("y", start = 0) +
	geom_text(aes(label = label_all, y = n_pos), size = 3) +
	geom_text(aes(label = total_label_all), x = 0.5, y = 0, size = 4) +
	xlim(c(0.5, 2.5)) +
	facet_grid(Taxon ~ .) +
	theme_void()

ggsave(filtered_file, p1, height = 10)
ggsave(all_file, p2, height = 10)

