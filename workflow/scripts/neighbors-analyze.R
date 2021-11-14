
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(ggrepel)
library(igraph)
library(ggnewscale)

input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

pfam2go_file <- unlist(input["pfam2go"])
clans_file   <- unlist(input["clans"])
labels_groups_file <- unlist(input["labels_group"])
classifications_file <- unlist(input["classifications"])

freq.threshold <- unlist(params["freq_threshold"])
gene.num.threshold <- unlist(params["gene_num_threshold"])

bubble_file <- unlist(output["bubble"])
volcano_file <- unlist(output["volcano"])
pooled_tests_file <- unlist(output["pooled_tests"])
all_hmms_file <- unlist(output["all_hmms"])
bubble_table_file <- unlist(output["bubble_table"])

hmm.names <- list(DTE = "DTE proton pumps", Heliorhodopsin = "Heliorhodopsins")

pfam2go <- read.fwf(pfam2go_file, width = c(5, 7, 10000), col.names = c("prefix", "hmm.acc", "go"), comment.char = "!") %>%
	extract(go, into = c("name.short", "go.desc", "go.id"), regex = "(.+) > GO:(.+) ; GO:(\\d+)") %>%
	select(hmm.acc, go.desc) %>%
	group_by(hmm.acc) %>%
	summarize(go.desc = paste(go.desc, collapse = "; "))
classifications <- read.table(classifications_file, col.names = c("hmm.acc", "hmm.category"))

label.desc <- read.table(clans_file, sep = "\t", quote = "", col.names = c("hmm.acc", "clan", "clan.name", "hmm.name", "hmm.desc")) %>%
	select(hmm.acc, hmm.desc, hmm.name, clan.name) %>%
	left_join(pfam2go, by = "hmm.acc")

labels.group <- read.table(labels_groups_file, header = T) %>%
	filter(target.num > 10) %>%
	mutate(hmm.acc = sub("[.][0-9]+$", "", hmm.acc), co.rate = co.express / label.num.all) %>%
	left_join(label.desc, by = "hmm.acc") %>%
	group_by(hmm.acc)

Shared.Taxa <- ungroup(labels.group) %>%
	distinct(hmm, Taxon) %>%
	group_by(Taxon) %>%
	filter(n() == 2) %>%
	distinct(Taxon) %>%
	pull

bkg.tests <- filter(labels.group, target.num.label > 1, target.num > 20) %>%
	mutate(In.Present = label.num, In.Absent = label.gene.num - label.num, Out.Present = bkg.num - label.num, Out.Absent = bkg.gene.num - label.gene.num - Out.Present) %>%
	ungroup %>%
	group_by(row_number()) %>%
	group_split %>%
	lapply(function(x) {
		with(x, paste(hmm.acc, Taxon, Group)) %>% print
		with(x, data.frame(row.names = c("In", "Out"), Present = c(In.Present, Out.Present), Absent = c(In.Absent, Out.Absent))) %>%
			fisher.test %>%
			with(data.frame(row.names = NULL, p.value.bkg = p.value, odds.ratio = estimate, conf.int.lo = conf.int[1], conf.int.hi = conf.int[2])) %>%
			bind_cols(select(x, hmm.acc, hmm, Taxon, Group, In.Present, In.Absent, Out.Present, Out.Absent))
	}) %>% bind_rows %>%
	mutate(q.value.bkg = p.adjust(p.value.bkg, method = "fdr"))

pooled.tests <- mutate(labels.group, Present = target.num.label, Absent = target.num - target.num.label) %>%
	group_by(hmm.acc, hmm) %>%
	summarize(AB = n_distinct(Group[Present > 1]) > 1, N_groups = sum(Present > 1), Present = sum(Present), Absent = sum(Absent), .groups = "drop_last") %>%
	filter(sum(Present > 1) > 0) %>%
	group_split %>%
	lapply(function(x) {
		print(first(x$hmm.acc))
		helio <- filter(x, hmm == "Heliorhodopsin")
		dte <- filter(x, hmm == "DTE")
		data.frame(row.names = c("Heliorhodopsin", "DTE"), Present = c(helio$Present, dte$Present), Absent = c(helio$Absent, dte$Absent)) %>%
			fisher.test %>%
			with(data.frame(row.names = NULL,
				select(helio, Absent, Present, N_groups, AB),
				select(dte, dte.Absent = Absent, dte.Present = Present, dte.N_groups = N_groups),
				p.value.types.pooled = p.value, odds.ratio.types.pooled = estimate, conf.int.lo.types.pooled = conf.int[1], conf.int.hi.types.pooled = conf.int[2],
				hmm.acc = first(x$hmm.acc)
			))
	}) %>% bind_rows %>%
	mutate(q.value.types.pooled = p.adjust(p.value.types.pooled, method = "fdr")) %>%
	mutate(Ratio = ifelse(odds.ratio.types.pooled > 1, Present / (Present + Absent), dte.Present / (dte.Present + dte.Absent))) %>%
	mutate(N_groups = ifelse(odds.ratio.types.pooled > 1, N_groups, dte.N_groups)) %>%
	mutate(Type = ifelse(odds.ratio.types.pooled > 1, "HeR > DTE", "HeR < DTE"))

left_join(pooled.tests, label.desc, by = "hmm.acc") %>%
	write.table(file = pooled_tests_file, row.names = F, sep = "\t")

all.tests <- mutate(labels.group, Present = target.num.label, Absent = target.num - target.num.label) %>%
	group_by(hmm.acc, Taxon) %>%
	filter(n_distinct(hmm) == 2) %>%
	group_by(hmm.acc) %>%
	group_split %>%
	lapply(function(x) {
		print(first(x$hmm.acc))
		helios <- filter(x, hmm == "Heliorhodopsin") %>%
			group_split(1:n())
		dtes <- filter(x, hmm == "DTE") %>%
			group_split(1:n())
		lapply(helios, function(helio) {
			lapply(dtes, function(dte) {
				if (helio$Taxon == dte$Taxon & (helio$Present > 1 | dte$Present > 1))
					data.frame(row.names = c("Heliorhodopsin", "DTE"), Present = c(helio$Present, dte$Present), Absent = c(helio$Absent, dte$Absent)) %>%
						fisher.test %>%
						with(data.frame(row.names = NULL,
							select(helio, Taxon, Group, Absent, Present),
							select(dte, dte.Absent = Absent, dte.Present = Present, dte.Taxon = Taxon),
							p.value.types = p.value, odds.ratio.types = estimate, conf.int.lo.types = conf.int[1], conf.int.hi.types = conf.int[2],
							hmm.acc = first(x$hmm.acc)
						))
			})
		}) %>% bind_rows
	}) %>% bind_rows %>%
	mutate(q.value.types = p.adjust(p.value.types, method = "fdr"))

pooled.best <- arrange(pooled.tests, q.value.types.pooled) %>%
	group_by(odds.ratio.types.pooled < 1) %>%
	left_join(label.desc, by = "hmm.acc") %>%
	filter(1:n() < 11)
helio.top <- filter(pooled.best, odds.ratio.types.pooled > 1) %>%
	pull(hmm.acc) %>%
	first
dte.top <- filter(pooled.best, odds.ratio.types.pooled < 1) %>%
	pull(hmm.acc) %>%
	first

p <- ggplot(pooled.tests, aes(x = log2(odds.ratio.types.pooled), y = -log10(q.value.types.pooled))) +
	geom_point() +
	geom_text_repel(aes(label = paste(hmm.acc, hmm.name, sep = " - ")), size = 3, pooled.best, point.size = 5) +
	xlab("log2 of estimated odds ratio HeR/DTE") +
	ylab("-log10 FDR") +
	facet_wrap(~ Type, scales = "free") +
	theme_bw()
ggsave(volcano_file, p)

labels.all <- labels.group %>%
	left_join(all.tests, by = c("hmm.acc", "Taxon", "Group")) %>%
	left_join(bkg.tests, by = c("hmm.acc", "Taxon", "hmm", "Group")) %>%
	mutate(freq = target.num.label / target.num)

left_join(labels.all, label.desc, by = "hmm.acc") %>%
	filter(freq > 0) %>%
	write.table(file = all_hmms_file, sep = "\t", row.names = F)

labels.plot <- filter(labels.all, target.num.label > 0) %>%
	group_by(hmm.acc, hmm) %>%
	mutate(odds.ratio.types.abs = ifelse(hmm == "Heliorhodopsin", odds.ratio.types, 1/odds.ratio.types)) %>%
	mutate(helio.over = hmm.acc %in% helio.top, dte.over = hmm.acc %in% dte.top) %>%
	group_by(hmm.acc) %>%
	mutate(bkg.freq.num = sum(hmm == "Heliorhodopsin" & freq > freq.threshold, na.rm = T), bkg.over = bkg.freq.num >= gene.num.threshold) %>%
	filter(any(bkg.over) | any(dte.over) | any(helio.over)) %>%
	mutate(category = case_when(any(helio.over) ~ "HeR", any(dte.over) ~ "DTE", any(bkg.over) ~ "Common")) %>%
	mutate(hmm.label = recode(hmm, !!!hmm.names)) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(Taxon.label = sprintf("%s (%s)", Taxon, unique(target.num))) %>%
	mutate(Label.label = paste(hmm.acc, hmm.desc, sep = " - ")) %>%
	ungroup %>%
	left_join(classifications, by = "hmm.acc") %>%
	arrange(hmm.category, hmm.acc) %>%
	mutate(row = row_number())

write.table(labels.plot, file = bubble_table_file, sep = "\t", row.names = F)

p <- ggplot(labels.plot, aes(x = reorder(Label.label, -row), y = Taxon.label)) +
	geom_point(aes(size = freq, color = co.rate)) +
	scale_colour_gradient2(low = "#56B4E9", mid = "gray", high = "#D55E00", midpoint = 0.5) +
	facet_grid(category ~ hmm.label + Group, scales = "free", space = "free", switch = "both") +
	coord_flip() +
	scale_x_discrete(position = "top") +
	theme_bw() +
	theme(
		axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
		strip.text.y = element_text(angle = 0),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()
	)

ggsave(bubble_file, p, width = 8.6, height = 6.2)
