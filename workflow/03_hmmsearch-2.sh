
# Second round of searches for rhodopsins

mkdir -p hmmsearch-2/{Heliorhodopsin,PR_XR} targets/{Heliorhodopsin,PR_XR}

search() {
	local faa=$1
	local contig=$(basename "$faa" .faa)
	local hmm=$(basename $(dirname "$faa"))
	local out=hmmsearch-2/$hmm/$contig
	hmmsearch "hmm/$hmm.hmm" "$faa" | tr \\n \\r | grep -v 'No hits detected' | tr \\r \\n > "$out.out"
	[[ -s "$out.out" && ! -s "$out.match" ]] && Rscript motifs.R "$out.out" | grep -wvf duplications.txt > "$out.match"
}
select_target() {
	local match=$1
	local pre=$(basename "$match" .match)
	local contig=${pre%-*}
	local hmm=$(basename $(dirname "$match"))
	local target=$(cut -f1 "$match" | head -n1)
	grep -wC10 "$target" "contigs/$hmm/$contig.gff" > "targets/$hmm/$target.gff"
	seqkit fx2tab "contigs/$hmm/$contig.faa" | grep -wC10 "$target" | seqkit tab2fx > "targets/$hmm/$target.faa"
}
export -f search select_target

find contigs     -name '*.faa'             | parallel search
find hmmsearch-2 -name '*.match' -size +1c | parallel select_target
