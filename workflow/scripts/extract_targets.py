
import re
from collections import defaultdict
from Bio import SeqIO

gff = str(snakemake.input["gff"])
faa = str(snakemake.input["faa"])
ignore = str(snakemake.input["ignore"])
matches = snakemake.input["matches"]
output_gff = str(snakemake.output["gff"])
output_faa = str(snakemake.output["faa"])

targets = []

for match in matches:
    with open(match) as fh:
        for line in fh:
            target, *rest = line.split()
            targets.append(target)

ignore_records = []
with open(ignore) as fh:
    for line in fh:
        ignore_records.append(line.strip())

contigs = {}
records = defaultdict(list)

fasta = {}
with open(faa) as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        fasta[record.id] = str(record.seq)

with open(gff) as fh:
    for line in fh:
        line = line.rstrip()
        if line and not line.startswith('#'):
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            match = re.search('ID=([\w_.-]+)', attributes)
            assert match, "ID attribute not found"
            ID = match.group(1)
            if ID not in ignore_records and ID in fasta:
                records[contig].append((contig, source, feature, start, end, score, strand, phase, attributes, ID))
                if ID in targets:
                    contigs[ID] = contig, len(records[contig]) - 1

fasta_written = []
with open(output_gff, 'w') as gff_fh:
    with open(output_faa, 'w') as faa_fh:
        for target, (contig, ix) in contigs.items():
            start = ix - 10 if ix >= 11 else 0
            stop = ix + 11
            for contig, source, feature, start, end, score, strand, phase, attributes, ID in records[contig][start:stop]:
                attributes += ';Target=%s' % target
                line = '\t'.join([contig, source, feature, start, end, score, strand, phase, attributes]) + '\n'
                gff_fh.write(line)
                if ID not in fasta_written:
                    faa_fh.write('>%s\n%s\n' % (ID, fasta[ID]))
                    fasta_written.append(ID)
