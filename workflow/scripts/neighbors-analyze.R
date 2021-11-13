
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
rate.threshold <- unlist(params["rate_threshold"])
to.target.abs.threshold <- unlist(params["to_target_abs_threshold"])
co.rate.threshold <- unlist(params["co_rate_threshold"])

freq.threashold <- unlist(params["freq_threashold"])
gene.num.threshold <- unlist(params["gene_num_threshold"])

bubble_file <- unlist(output["bubble"])
volcano_file <- unlist(output["volcano"])
table_file <- unlist(output["table"])

hmm.names <- list(PR_XR = "DTE proton pumps", Heliorhodopsin = "Heliorhodopsins")

pfam2go <- read.fwf(pfam2go_file, width = c(5, 7, 10000), col.names = c("prefix", "hmm.acc", "go"), comment.char = "!") %>%
	extract(go, into = c("name.short", "go.desc", "go.id"), regex = "(.+) > GO:(.+) ; GO:(\\d+)") %>%
	select(hmm.acc, go.desc) %>%
	group_by(hmm.acc) %>%
	summarize(go.desc = paste(go.desc, collapse = "; "))

label.desc <- read.table(clans_file, sep = "\t", quote = "", col.names = c("hmm.acc", "clan", "clan.name", "hmm.name", "hmm.desc")) %>%
	select(hmm.acc, hmm.desc, hmm.name, clan.name) %>%
	left_join(pfam2go, by = "hmm.acc")

labels.group <- read.table(labels_groups_file, header = T) %>%
	filter(target.num > 10) %>%
	mutate(hmm.acc = sub("[.][0-9]+$", "", hmm.acc), co.rate = co.express / label.num.all) %>%
	left_join(label.desc, by = "hmm.acc") %>%
	group_by(hmm.acc)
	# filter(any(freq > rate.threshold), any(to.target.abs.avg < to.target.abs.threshold), any(co.rate > co.rate.threshold))

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
		pr_xr <- filter(x, hmm == "PR_XR")
		data.frame(row.names = c("Heliorhodopsin", "PR_XR"), Present = c(helio$Present, pr_xr$Present), Absent = c(helio$Absent, pr_xr$Absent)) %>%
			fisher.test %>%
			with(data.frame(row.names = NULL,
				select(helio, Absent, Present, N_groups, AB),
				select(pr_xr, pr_xr.Absent = Absent, pr_xr.Present = Present, pr_xr.N_groups = N_groups),
				p.value.types.pooled = p.value, odds.ratio.types.pooled = estimate, conf.int.lo.types.pooled = conf.int[1], conf.int.hi.types.pooled = conf.int[2],
				hmm.acc = first(x$hmm.acc)
			))
	}) %>% bind_rows %>%
	mutate(q.value.types.pooled = p.adjust(p.value.types.pooled, method = "fdr")) %>%
	mutate(Ratio = ifelse(odds.ratio.types.pooled > 1, Present / (Present + Absent), pr_xr.Present / (pr_xr.Present + pr_xr.Absent))) %>%
	mutate(N_groups = ifelse(odds.ratio.types.pooled > 1, N_groups, pr_xr.N_groups)) %>%
	mutate(Type = ifelse(odds.ratio.types.pooled > 1, "HeR > DTE", "HeR < DTE"))

all.tests <- mutate(labels.group, Present = target.num.label, Absent = target.num - target.num.label) %>%
	group_by(hmm.acc, Taxon) %>%
	filter(n_distinct(hmm) == 2) %>%
	group_by(hmm.acc) %>%
	group_split %>%
	lapply(function(x) {
		print(first(x$hmm.acc))
		helios <- filter(x, hmm == "Heliorhodopsin") %>%
			group_split(1:n())
		pr_xrs <- filter(x, hmm == "PR_XR") %>%
			group_split(1:n())
		lapply(helios, function(helio) {
			lapply(pr_xrs, function(pr_xr) {
				if (helio$Taxon == pr_xr$Taxon & (helio$Present > 1 | pr_xr$Present > 1))
					data.frame(row.names = c("Heliorhodopsin", "PR_XR"), Present = c(helio$Present, pr_xr$Present), Absent = c(helio$Absent, pr_xr$Absent)) %>%
						fisher.test %>%
						with(data.frame(row.names = NULL,
							select(helio, Taxon, Group, Absent, Present),
							select(pr_xr, pr_xr.Absent = Absent, pr_xr.Present = Present, pr_xr.Taxon = Taxon),
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
pr_xr.top <- filter(pooled.best, odds.ratio.types.pooled < 1) %>%
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

labels.plot <- labels.group %>%
	left_join(all.tests, by = c("hmm.acc", "Taxon", "Group")) %>%
	left_join(bkg.tests, by = c("hmm.acc", "Taxon", "hmm", "Group")) %>%
	mutate(freq = target.num.label / target.num) %>%
	filter(target.num.label > 0) %>%
	group_by(hmm.acc, hmm) %>%
	mutate(odds.ratio.types.abs = ifelse(hmm == "Heliorhodopsin", odds.ratio.types, 1/odds.ratio.types)) %>%
	group_by(hmm.acc, hmm) %>%
	{ write.table(., table_file, row.names = T); (.) } %>%
	mutate(helio.over = hmm.acc %in% helio.top, pr_xr.over = hmm.acc %in% pr_xr.top) %>%
	group_by(hmm.acc) %>%
	mutate(bkg.freq.num = sum(hmm == "Heliorhodopsin" & freq > freq.threashold, na.rm = T), bkg.over = bkg.freq.num >= gene.num.threshold) %>%
	filter(any(bkg.over) | any(pr_xr.over) | any(helio.over)) %>%
	mutate(category = case_when(any(helio.over) ~ "HeR", any(pr_xr.over) ~ "DTE", any(bkg.over) ~ "Common")) %>%
	mutate(hmm.label = recode(hmm, !!!hmm.names)) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(Taxon.label = sprintf("%s (%s)", Taxon, unique(target.num))) %>%
	mutate(Label.label = paste(hmm.acc, hmm.desc, sep = " - ")) %>%
	ungroup %>%
	arrange(bkg.freq.num) %>%
	mutate(row = row_number())

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

ggsave(bubble_file, p, width = 8.5, height = 5.4)
