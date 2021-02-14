
# Filter hmm matches from the first round

mkdir -p contigs/{Heliorhodopsin,PR_XR}

predict_genes() {
	local wgs=$1
	local contig=$2
	mkdir -p "genemarks/$contig"
	samtools faidx "wgs/$wgs.fna" "$contig" > "genemarks/$contig/contig.fna"
	(cd "genemarks/$contig"; gmsn.pl --format GFF --prok --output /dev/stdout contig.fna | awk -F\\t '{sub(/gene_id=/,sprintf("ID=%s_gmsn_",$1),$9)}1' OFS=\\t > "../$contig.gff")
}
extract_contigs() {
	local match=$1
	local wgs=$(basename "$match" .match)
	local hmm=$(basename $(dirname "$match"))
	local contig pre
	cut -f1 "$match" | sed -E 's/_[0-9]+$//' | sort -u | while read contig; do
		pre=contigs/$hmm/$contig
		[[ -s "wgs/$wgs.gff" || -s "genemarks/$contig.gff" ]] || predict_genes "$wgs" "$contig"
		if [[ -s "genemarks/$contig.gff" ]]; then
			gffread -H -y "$pre.faa" -g "wgs/$wgs.fna" "genemarks/$contig.gff"
			awk -vc="$contig" '$1==c' "genemarks/$contig.gff" | sort -k4,4n > "$pre.gff"
		else
			awk -vc="$contig" '$1==c' "wgs/$wgs.gff" | sort -k4,4n | tee "$pre.gff" | cut -f9 | grep -oP 'ID=\K([a-zA-Z0-9_.-]+)\b' | xargs samtools faidx "wgs/$wgs.faa" > "$pre.faa"
		fi
	done
}
export -f extract_contigs predict_genes

find hmmsearch-1 -name '*.match' -size +1c | parallel extract_contigs
