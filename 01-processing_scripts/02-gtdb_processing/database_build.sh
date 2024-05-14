conda create --name 2023-rpoCdb

conda activate 2023-rpoCdb
conda install -c bioconda prokka
conda install -c bioconda seqtk

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

## Pull global GTDB files into a new folder for ease of access
cd /home/allie/rpoCdb/gtdb_version1.0
mkdir gtdb_prokka_concat && cd gtdb_genomes_reps_r214
# pull genbank file for all samples into single genbank
find -type f -regex ".*\.gff" -exec cat {} + > /home/allie/rpoCdb/gtdb_version1.0/gtdb_prokka_concat/all.gff &
find -type f -regex ".*\.fna" -exec cat {} + > /home/allie/rpoCdb/gtdb_version1.0/gtdb_prokka_concat/all.fna &
find -type f -regex ".*\.gbf" -exec cat {} + > /home/allie/rpoCdb/gtdb_version1.0/gtdb_prokka_concat/all.gbf &

# # next need to fix the genbank file as it has an issue with the locus ids colliding with the sequence length
# # find pattern that is messing up our search parameters
# cat all.gbf | grep "LOCUS" > loci_lines
# grep "\.[0-9][0-9]* bp" loci_lines | awk '{print $2}' > problemIDs
# cat problemIDs | sed 's/\.1/\.1 /' > fixIDs

# # ok this is an easy ish fix, need to add a space after the locus ID and base pair number
# # IN FUTURE RELEASES, ADD THIS TO PROKKA CODE -- NEEDS TO HAVE CONTIG IDS LESS THAN 20 CHARACTERS
# paste problemIDs fixIDs | \
# 	sed "s/^/'0,\//" | \
# 	sed 's/\./\\./g' | \
# 	sed 's/\t/\/s\/\//' | \
# 	sed "s/$/\/'/" | \
# 	sed 's/^/sed -i /' | \
# 	sed 's/$/ all.gbk/'> tofix.sh
# bash tofix.sh # run corrections
# # see if this was fixed
# grep "LOCUS" all.gbk > loci_lines.after
# grep "\.[0-9][0-9]* bp" loci_lines.after

# running fragant script to pull single and multicopy reads from genbank file
time python3 fragant.py 1>fragant.out 2>fragant.err

# pull fasta from coordinates generated from fragant
awk -F"," '{print $2 "\t" $4 "\t" $5}' multicopy_list.csv > multi.query
awk -F"," '{print $2 "\t" $4 "\t" $5}' singlecopy_list.csv > single.query

# subseq
seqtk subseq all.fna multi.query > rpoC_multi.fa &
seqtk subseq all.fna single.query > rpoC_single.fa &

# concatenate together to dereplicate
cat rpoC_single.fa rpoC_multi.fa > rpoC_all.fa
conda install -c bioconda vsearch=2.28.1
vsearch --threads 25 --sizeout --derep_fulllength rpoC_all.fa --output rpoC_all.uniq.fa

# pull taxonomy information from each accession id
cd .. && mkdir rpoCdb_files && cd rpoCdb_files
mv ../gtdb_prokka_concat/rpoC_all.uniq.fa .
grep ">" rpoC_all.uniq.fa | awk -F":" '{print $1}' | sed 's/>//'  > accessions

# get metadata for each entry (metadata file from suzanne)
cat accessions | while read line; do grep $line -m1 metadata.tsv; done > metadata.filt.tsv &

# merge with our list of sequences to get full metadata for each sequence
# in ipython
import pandas as pd
access = pd.read_table("accessions", header=None, names=["accession"])
meta = pd.read_table("metadata.filt.tsv", header=None, names=["contig", "gca", "assembely_level", "cord1", "cord2", "length", "kingdom", "phylum", "class", "order", "family", "genus", "species", "rpoC_gca", "contig_number", "rpoC_num"])
merge = pd.merge(access, meta, left_on="accession", right_on="contig")
merge.to_csv("metadata.merge.txt", sep="\t", index=False, na_rep="NA")
# get list of only bacterial genes to filter
grep -v "Archaea" metadata.merge.txt > metadata.bac.txt
# filter out archaea
awk '{print $1}' metadata.bac.txt | grep -v "accession" > bac.ids
sed -i 's/:/ /g' rpoC_all.uniq.fa
seqtk subseq rpoC_all.uniq.fa bac.ids > rpoC_all.bac.fa
# now get length distribution
cat rpoC_all.bac.fa | awk '$0 ~ ">" {if (NR > 1) \
	{print c;} c=0;printf substr($0,2,100) "\t"; } \
	$0 !~ ">" {c+=length($0);} END { print c; }' > \
	rpoC_all.bac.lengths
awk -F" " '{print $3}' rpoC_all.bac.lengths > len
# quick analysis of length distribution in R to inform cutoffs
# below in R
dat <- read.table("len", header=F)
summary(dat)
 #       V1
 # Min.   :  126
 # 1st Qu.: 3516
 # Median : 3939
 # Mean   : 3540
 # 3rd Qu.: 4224
 # Max.   :49951
# quick visualization
pdf("bac.length.density.pdf")
d <- density(dat$V1)
plot(d)
dev.off()
sd(dat$V1)
# [1] 1140.27
# ok, so length filter is going to be the 1st quartile - sd, 3rd quartile + sd
vsearch --fastx_filter rpoC_all.bac.fa \
	--fastq_maxlen 5364 \
	--fastq_minlen 2376 \
	--lengthout \
	--fastaout rpoC_all.lenfilt.fa

# finally, filter our metadata to only include those sequences we kept
grep ">" rpoC_all.lenfilt.fa | awk '{print $1}' | sed 's/>//' > filtered.ids
cat filtered.ids| while read line; do grep $line metadata.merge.txt; done > metadata.lenfilt.txt &









