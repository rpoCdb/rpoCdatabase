### rpoCdb Database Build 

```bash
conda create --name 2023-rpoCdb

conda activate 2023-rpoCdb
conda install -c bioconda prokka
conda install -c bioconda agat

cd /home/allie/rpoCdb/gtdb_version1.0 # on hillary

# download required gtdb files
wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz &
wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
# pull ncbi refseq and genbank files for comparison as well
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
# unzip files
gzip -d bac120_metadata.tsv.gz
tar xzf gtdb_genomes_reps.tar.gz

# pull all genomes into single directory
find . -name "*fna.gz" | parallel -j30 'mv {} gtdb_genomes_reps_r214'
# clean up
cd gtdb_genomes_reps_r214
find -type d -empty -print | xargs rm -r
rm -r database

# unzip files for prokka
find -name "*fna.gz" | parallel 'gzip -d {}'

# run prokka on each genome
for f in gtdb_genomes_reps_r214/*fna;
	do; 
		prokka --cpus 55 \
		--force \
		--outdir $f.prokka \
		--addgenes \
		--norrna \
		--notrna \
		$f
	;done

# if prokka is interrupted, re run with force flag disabled to prevent overwriting previously generated output directories
for f in gtdb_genomes_reps_r214/*fna;
	do; 
		prokka --cpus 55 \
		--outdir $f.prokka \
		--addgenes \
		--norrna \
		--notrna \
		$f
	;done 1>prokka_restart.out 2>prokka_restart.err

##### DOMHAIN RNA Seq Database build #####

# make global GTDB gff and fna file (for RNA seq mapping -- domhain project)
cd /home/allie/rpoCdb/gtdb_version1.0
mkdir gtdb_prokka_gff && cd gtdb_genomes_reps_r214
find -type f -regex ".*\.gff" -exec cat {} + > /home/allie/rpoCdb/gtdb_version1.0/gtdb_prokka_all/all.gff
find -type f -regex ".*\.fna" -exec cat {} + > /home/allie/rpoCdb/gtdb_version1.0/gtdb_prokka_gff/all.fna
# tarball for transfer
tar czvf all.fna.tar.gz all.fna &
tar czvf all.gff.tar.gz all.gff
# transfer to palmetto
rsync -azP all*tar.gz amann3@hpcdtn01.rcd.clemson.edu:/scratch1/amann3/gtdb/

# pull all annotations for rpoc from gff files
cd /home/allie/rpoCdb/gtdb_version1.0
mkdir rpoC_annotated
cd rpoC_annotated

for f in ../gtdb_genomes_reps_r214/*prokka/*gff;
	do;
		grep "eC_number=2.7.7.6;Name=rpoC" $f
	;done > rpoC.gff &
```





