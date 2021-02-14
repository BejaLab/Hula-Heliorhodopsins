#!/bin/bash

# This script generates data for calulcation of background
# frequencies of matches

if [ -z "$PFAMDB" ]; then
	echo "Environmental variable PFAMDB not set" >&2
	exit 1
fi

mkdir -p background

# Select assemblies of sufficient length (> 999 genes)
find hmmsearch-1 -name '*.match' -size +1c | parallel [ -s wgs/{/.}.faa ] '&&' seqkit stat -T wgs/{/.}.faa \| tail -1 | parallel --colsep \\t [ {4} -gt 999 ] '&&' cp {1} background/

# Subset pfam matches - take those that have appear at least 3 times
mkdir targets/Pfam
awk -vRS=// -F\\n '{$1=$1}1' "$PFAMDB.dat" | \
	grep -wf <(awk '!/^#/{print$6,$15}' targets/*/*.pfam | sort | uniq -c | awk '$1>2&&$2{print$2,$3}' OFS=\\n) | \
	grep -Eo PF[0-9]+[.][0-9]+ | \
	hmmfetch -f "$PFAMDB" - > targets/Pfam/Pfam-A.hmm
hmmpress targets/Pfam/Pfam-A.hmm
ln -fs "$PFAMDB.dat" targets/Pfam/Pfam-A.hmm.dat

# Search the selected assemblies with the Pfam subset
parallel [ -s {.}.pfam ] '||' pfam_scan.pl -cpu 1 -dir targets/Pfam -fasta {} -outfile {.}.pfam ::: background/*.faa

# Search the selected assemblies for de novo orthogroups
parallel [ -s {.}.blast ] '||' usearch -ublast {} -db targets/targets.8.ublast -threads 1 -evalue 1e-05 -blast6out {.}.blast ::: background/*.faa
