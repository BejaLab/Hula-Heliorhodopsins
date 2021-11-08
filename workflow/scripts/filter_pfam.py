
pfam_txt = snakemake.input["pfam_txt"]
dat = str(snakemake.input["dat"])
output = str(snakemake.output)
pfam_num_min = int(snakemake.params["pfam_num_min"])

clans = {}
names = {}

for fname in pfam_txt:
	with open(fname) as fh:
		for line in fh:
			line = line.rstrip()
			if line and not line.startswith('#'):
				seq_id, aln_start, aln_end, env_start, env_end, hmm_acc, hmm_name, hmm_type, hmm_start, hmm_end, hmm_length, bit_score, E_value, significance, clan = line.split()
				if hmm_name in names:
					names[hmm_name] += 1
				else:
					names[hmm_name] = 1
				if clan in clans:
					clans[clan] += 1
				else:
					clans[clan] = 1

collect = False

with open(output, 'w') as out:
	with open(dat) as fh:
		for line in fh:
			line = line.rstrip()
			if line.startswith('#=GF'):
				prefix, tag, value = line.split(maxsplit = 2)
				if tag == "AC":
					AC = value
				elif tag == "ID":
					collect |= value in names and names[value] >= pfam_num_min
				elif tag == "CL":
					collect |= value in clans and clans[value] >= pfam_num_min
			elif line == "//":
				if collect:
					out.write(AC + '\n')
					collect = False
