# Hula-Heliorhodopsins

This is the gene neighbors workflow, it includes two parts:

* initial screen for heliorhodopsins and DTE rhodopsins: `snakemake -c[cores] -s workdlow/1-Search.snake`
* actual analysis: `snakemake -c[cores] -s workdlow/2-Analysis.snake`

Not included in the repository are:

* directories `wgs` and `database` and the derived data stored in `contigs` and `background`
* directory `Pfam` - it contains the Pfam database from `http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/` (the released used originally is 33.1)
