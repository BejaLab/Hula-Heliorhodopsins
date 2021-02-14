
# Find orthogroups among rhodopsin-containing contigs
seqkit rmdup targets/*/*.faa > targets/targets.faa 
cdhit -i targets/targets.faa -o targets/targets.cdhit -c 0.8 -d 0 -T 0
(cd targets; proteinortho -sim=0.1 -selfblast -project=targets -p=blastp+ -keep -conn=0.05 targets.cdhit)

# Search for Pfam matches
parallel pfam_scan.pl -cpu 1 -dir /Data/seqdbases/Pfam -fasta {} -outfile {.}.pfam ::: targets/*/*.faa
