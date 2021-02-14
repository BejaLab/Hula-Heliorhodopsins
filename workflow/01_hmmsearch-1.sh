
# First round of hmmsearches

mkdir -p hmmsearch-1/{Heliorhodopsin,PR_XR}

# Search for the two rhodopsin types among ORFs
search() {
	hmmsearch "hmm/$1.hmm" - | tr \\n \\r | grep -v 'No hits detected' | tr \\r \\n
}
# output re-formatted protein sequences from a faa file
faa_orfs() {
	seqkit fx2tab "$1" | grep -wvf blacklist.txt | awk '{ printf ">%s_%d\n%s\n", $2, NR, $NF }'
}
# generate orfs from fna file
fna_orfs() {
	getorf -filter < "$1"
}
export -f search fna_orfs faa_orfs

# Search for rhodopsins
cat {wgs,gbk}-*.txt | parallel --tag [[ ! -s wgs/{}.fna '||' -f hmmsearch-1/PR_XR/{}.out          ]] '||' fna_orfs wgs/{}.fna \| search PR_XR          \> hmmsearch-1/PR_XR/{}.out
cat {wgs,gbk}-*.txt | parallel --tag [[ ! -s wgs/{}.fna '||' -f hmmsearch-1/Heliorhodopsin/{}.out ]] '||' fna_orfs wgs/{}.fna \| search Heliorhodopsin \> hmmsearch-1/Heliorhodopsin/{}.out
cat faa-*.txt | parallel --tag [ -f hmmsearch-1/Heliorhodopsin/{}.out ] '||' faa_orfs wgs/{}.faa \| search Heliorhodopsin \> hmmsearch-1/Heliorhodopsin/{}.out

# Extract motifs
find hmmsearch-1 -name '*.out' -size +1c | parallel [ -f {.}.match ] '||' Rscript motifs.R {} {//} \> {.}.match
