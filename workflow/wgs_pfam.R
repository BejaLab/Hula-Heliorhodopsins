
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(ggrepel)
library(igraph)
library(ggnewscale)

label.types <- "pfam"

hmm.names <- list(PR_XR = "DTE proton pumps", Heliorhodopsin = "Heliorhodopsins")

pfam2go <- read.fwf("pfam2go.txt", width = c(5, 7, 10000), col.names = c("prefix", "label", "go"), comment.char = "!") %>%
	extract(go, into = c("name.short", "go.desc", "go.id"), regex = "(.+) > GO:(.+) ; GO:(\\d+)") %>%
	select(label, go.desc) %>%
	group_by(label) %>%
	summarize(go.desc = paste(go.desc, collapse = "; "))
label.desc <- read.table("Pfam-A.clans.tsv", sep = "\t", quote = "", col.names = c("hmm.acc", "clan", "clan.name", "hmm.name", "hmm.desc")) %>%
	{bind_rows(select(., label = hmm.acc, desc = hmm.desc, clan.name), distinct(., label = clan, desc = clan.name, clan.name))} %>%
	left_join(pfam2go, by = "label")

rate.threshold <- 0.05
to.target.abs.threshold <- 4
co.rate.threshold <- 0.1
labels.group <- read.table("labels.group.txt", header = T) %>%
	filter(label.type %in% label.types, target.num > 10) %>%
	mutate(label = sub("^(PF[0-9]+)[.][0-9]+$", "\\1", label), co.rate = co.express / label.num.all) %>%
	left_join(label.desc, by = "label") %>%
	filter(label.type == "pfam") %>%
	group_by(label) %>%
	filter(any(target.num.label / target.num > rate.threshold), any(to.target.abs.avg < to.target.abs.threshold), any(co.rate > co.rate.threshold))

all.tests <- labels.group %>%
	mutate(Present = target.num.label, Absent = target.num - target.num.label) %>%
	group_split %>%
	lapply(function(x) {
		print(first(x$label))
		helios <- filter(x, hmm == "Heliorhodopsin") %>%
			bind_rows(summarize(., Taxon = "All", Group = "", Present = sum(Present), Absent = sum(Absent))) %>%
			group_split(1:n())
		pr_xrs <- filter(x, hmm == "PR_XR") %>%
			bind_rows(summarize(., Taxon = "All", Present = sum(Present), Absent = sum(Absent))) %>%
			mutate(Present.rank = rank(Present, tie = "first"), Ratio.rank = rank(Present/Absent, tie = "first")) %>%
			filter(Present.rank == max(Present.rank) | Ratio.rank == max(Ratio.rank) | Taxon == "All") %>%
			group_split(1:n())
		lapply(helios, function(helio) {
			lapply(pr_xrs, function(pr_xr) {
				data.frame(row.names = c("Heliorhodopsin", "PR_XR"), Present = c(helio$Present, pr_xr$Present), Absent = c(helio$Absent, pr_xr$Absent)) %>%
					fisher.test %>%
					with(data.frame(row.names = NULL,
						select(helio, Taxon, Group, Absent, Present),
						select(pr_xr, pr_xr.Absent = Absent, pr_xr.Present = Present, pr_xr.Taxon = Taxon),
						p.value, odds.ratio = estimate, conf.int.lo = conf.int[1], conf.int.hi = conf.int[2],
						label = first(x$label) #, desc = first(x$desc), clan.name = first(x$clan.name)
					))
			}) %>% bind_rows
		}) %>% bind_rows
	}) %>% bind_rows

bkg.tests <- filter(labels.group, label.num > 0, bkg.num > 0) %>%
	mutate(In.Present = label.num, In.Absent = label.gene.num - label.num, Out.Present = bkg.num - label.num, Out.Absent = bkg.gene.num - label.gene.num - Out.Present) %>%
	ungroup %>%
	group_by(row_number()) %>%
	group_split %>%
	lapply(function(x) {
		with(x, paste(label, Taxon, Group)) %>% print
		with(x, data.frame(row.names = c("In", "Out"), Present = c(In.Present, Out.Present), Absent = c(In.Absent, Out.Absent))) %>%
			fisher.test %>%
			with(data.frame(row.names = NULL, p.value, odds.ratio = estimate, conf.int.lo = conf.int[1], conf.int.hi = conf.int[2])) %>%
			bind_cols(select(x, label, hmm, Taxon, Group, In.Present, In.Absent, Out.Present, Out.Absent))
	}) %>% bind_rows

helio.tests <- filter(all.tests, Taxon != "All") %>%
	arrange(-p.value) %>%
	distinct(Taxon, Group, label, .keep_all = T) %>%
	mutate(hmm = "Heliorhodopsin")

pr_xr.tests <- filter(all.tests, Taxon == "All", pr_xr.Taxon != "All") %>%
	arrange(-p.value) %>%
	distinct(pr_xr.Taxon, label, .keep_all = T) %>%
	mutate(hmm = "PR_XR", Group = "")

freq.threashold1 <- 0.05
freq.threashold2 <- 0.02
q.value.threashold <- 0.01

labels.plot <- left_join(labels.group, bind_rows(helio.tests, pr_xr.tests), by = c("label", "hmm", "Group"), suffix = c("",".types")) %>%
	filter(hmm == "PR_XR" | Taxon == Taxon.types) %>%
	left_join(bkg.tests, by = c("label", "Taxon", "hmm", "Group"), suffix = c(".types", ".bkg")) %>%
	filter(target.num.label > 0) %>%
	mutate(q.value.types = p.adjust(p.value.types), q.value.bkg = p.adjust(p.value.bkg)) %>%
	group_by(label, hmm) %>%
	mutate(freq = target.num.label / target.num) %>%
	mutate(odds.ratio.types.abs = ifelse(hmm == "Heliorhodopsin", odds.ratio.types, 1/odds.ratio.types)) %>%
	mutate(helio.over =
		hmm == "Heliorhodopsin" &
		sum(freq > freq.threashold1, na.rm = T) >= 2 &
		sum(freq > freq.threashold2, na.rm = T) >= 4 &
		sum(q.value.types < q.value.threashold, na.rm = T) >= 3 &
		sum(co.rate > 0.5, na.rm = T) >= 3
	) %>%
	mutate(pr_xr.over =
		hmm == "PR_XR" &
		sum(freq > freq.threashold1, na.rm = T) >= 1 &
		sum(freq > freq.threashold2, na.rm = T) >= 3 &
		sum(q.value.types < q.value.threashold, na.rm = T) >= 2 &
		sum(co.rate > 0.5, na.rm = T) >= 2
	) %>%
	group_by(label) %>%
	mutate(bkg.over =
		sum(freq > freq.threashold1, na.rm = T) >= 3 &
		sum(freq > freq.threashold2, na.rm = T) >= 5 &
		sum(q.value.bkg < q.value.threashold, na.rm = T) >= 4 &
		sum(co.rate > 0.5, na.rm = T) >= 4
	) %>%
	{ write.table(., "neighbors.txt", row.names = T); (.) } %>%
	filter(any(bkg.over) | any(pr_xr.over) | any(helio.over)) %>%
	mutate(significant.types = sum(q.value.types < 0.01, na.rm = T), significant.bkg = sum(q.value.bkg < 0.01, na.rm = T)) %>%
	mutate(category = case_when(any(helio.over) ~ "HeR>DTE", any(pr_xr.over) ~ "HeR<DTE", any(bkg.over) ~ "Common")) %>%
	mutate(hmm.label = recode(hmm, !!!hmm.names)) %>%
	group_by(hmm, Taxon, Group) %>%
	mutate(Taxon.label = sprintf("%s (%s)", Taxon, unique(target.num))) %>%
	mutate(Label.label = paste(label, desc, sep = " - ")) %>%
	ungroup %>%
	arrange(significant.types, significant.bkg) %>%
	mutate(row = row_number())

p <- ggplot(labels.plot, aes(x = reorder(Label.label, row), y = Taxon.label)) +
	geom_point(
		size = 7,
		data = filter(labels.plot, q.value.bkg < 0.01),
		position = position_nudge(0, 0), shape = 1, color = "blue"
	) +
	geom_point(
		size = 9,
		data = filter(labels.plot, q.value.types < 0.01 & (helio.over | pr_xr.over)),
		position = position_nudge(0, 0), shape = 1, color = "red"
	) +
	geom_point(aes(size = freq, color = co.rate)) +
	scale_color_continuous(type = "viridis") +
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
ggsave("neighbors.pdf", p, width = 8.5, height = 6.5)
