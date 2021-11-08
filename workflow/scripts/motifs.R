
if (exists("snakemake")) {
	hmmsearch  <- as.character(snakemake@input["hmmsearch"])
	motif_file <- as.character(snakemake@input["motif_file"])
	output     <- as.character(snakemake@output)
} else {
	args <- commandArgs(T)
	hmmsearch  <- args[1]
	motif_file <- args[2]
	output     <- args[3]
}

info <- file.info(hmmsearch)
if (info$size == 0) {
	cat("", file = output)
} else {
	library(bioformatr)
	suppressPackageStartupMessages(library(tidyr))
	suppressPackageStartupMessages(library(dplyr))

	motif <- as.integer(readLines(motif_file))

	sequences <- read.hmmer(hmmsearch) %>%
		filter(inclusion == "!") %>%
		group_by(Sequence, domain) %>%
		separate_rows(M.Seq, T.Seq, sep = "") %>%
		filter(! M.Seq %in% c(".","")) %>%
		mutate(pos = hmmfrom - 1 + 1:n()) %>%
		filter(pos %in% motif, M.Seq == T.Seq) %>%
		group_by(Sequence) %>%
		filter(n() == length(motif)) %>%
		distinct(Sequence) %>%
		pull
	if (length(sequences) > 0) {
		cat(sequences, file = output, sep = "\n")
	} else {
		cat("", file = output)
	}
}
