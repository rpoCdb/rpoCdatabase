# pull genbank file for all samples into single genbank
```sh
find . -maxdepth 2 -type f -name "*.gff" -exec cat {} + > ~/rpoCdb/01-processing_scripts/rpoCdbv2.0/gtdb_prokka_concat/all.gff &

find -type f -regex ".*\.fna" -exec cat {} + > ~/rpoCdb/01-processing_scripts/rpoCdbv2.0/gtdb_prokka_concat/all.fna &
find -type f -regex ".*\.gbf" -exec cat {} + > ~/rpoCdb/01-processing_scripts/rpoCdbv2.0/gtdb_prokka_concat/all.gbf &
```
# Get rpoC
```sh
grep "rpoC" all.gff | grep -v "gene=rpoC2" | grep -v "gamma" > rpoC.gff # 210844 non-rpoC2 genes, gamma subunits are found in cholorplasts
cd .. && mkdir rpoCdb_files && cd rpoCdb_files
#get number greater than a certain amount
grep "gene=rpoC_" ../gtdb_prokka_concat/rpoC.gff | sed 's/.*gene=rpoC_//' | sed 's/;.*//' | sort -n | uniq > number_annots.txt
```
## Get columns with accesion number and coordinates from gff file
```sh
awk '{print $1}' ../gtdb_prokka_concat/rpoC.gff | sort | uniq > rpoC_seqid #get just seqid
#get cords from genomes with just one rpoC
rm cords
parallel -a rpoC_seqid -j 50 -k "grep -Fw '{}' ../gtdb_prokka_concat/rpoC.gff"> cords #get cords of all rpoC genes
wc -l ../gtdb_prokka_concat/rpoC.gff # 210844
wc -l cords # 210844
awk '{print $1, $4, $5}' cords | sed 's/ /\t/g'> rpoC_cords
#get seqid to accescion number
# scp -r ./cords suzanne@hillary.clemson.edu:/home/suzanne/rpoCdb/cords
find /home/allie/rpoCdb/gtdb_version2.0/gtdb_genomes_reps_r220/ -type f -name "*.ffn" | wc -l # 113104
# find /home/allie/rpoCdb/gtdb_version2.0/gtdb_genomes_reps_r220 -type f  -name "*.gff" -exec grep -HF rpoC {} + >rpoC_files 
find /home/allie/rpoCdb/gtdb_version2.0/gtdb_genomes_reps_r220 -type f -name "*.gff" | xargs -P 60 -I {} grep -H -F 'rpoC' {} > rpoC_files
grep allie rpoC_files > rpoC_files2
sed 's/_genomic.*gff:/\t/' rpoC_files2 | sed 's/.*r220\///' | sed 's/\tprokka.*//' | sed 's/\tProdigal.*//' | sort | uniq  > acces2seqid
# scp -r ./acces2seqid suzanne@stella.clemson.edu:~/rpoCdb/rpoC_annotated
wc -l acces2seqid #98858
#get cordinates and seqid
python3 cord_flip.py #run python script that edits rpoC_cords. flip_cords.tsv output
sed '1d' flip_cords.tsv | awk '{print $2, $3, $4, $5}' | sed 's/ /\t/g'> temp
mv temp flip_cords.tsv

#join cordinate file to accescion number
python3 seqid2acces.py #outputs cords.tsv
#check for any missing accesion numbers

awk '{print $1}' cords.tsv | sed '1d '> annot_cords
awk '{print $1}' cords | sort | uniq >seqids_needed
python3 miss_ids.py #output get_seqs.tsv
awk '{print $1}' get_seqs.tsv | wc -l # should be empty
```
Get the full taxonomy using asscesnion number
```sh
#get archea taxonomy
wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz
gzip -d ar53_metadata.tsv.gz
awk -F "\t" '{print $1, $49, $20}' ar53_metadata.tsv | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' > accs2taxa
#get bacteria taxonomy
awk -F "\t" '{print $1, $49, $20}' bac120_metadata.tsv | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa #should be missing but had 0
cat ar53_metadata.tsv bac120_metadata.tsv > full_metadata.tsv
parallel -a missing_taxa -j 50 -k "grep -wm 1 '{}' full_metadata.tsv >> missed_taxa"
awk -F "\t" '{print $55, $49, $20}' missed_taxa | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa second round
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa # 0
awk -F "\t" '{print $1, $49, $20}' full_metadata.tsv | sed '1d' | sed 's/[GR][BS]_G/G/' | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa #should be 0
sort <(sed '1d' taxa.tsv) | uniq > taxa_uniq.tsv
cat <(head -n 1 taxa.tsv) taxa_uniq.tsv > temp
mv temp taxa.tsv
```
## Get metadata file ready
Prepare file for R
```sh
awk '{print $1, $5, $6, $2, $3, $4, $7}' taxa.tsv | sed '1d' | sed 's/ /\t/g' | sort | uniq > length_taxa
grep -c ";s__" length_taxa #numbers should match
grep -c "d__" length_taxa
sed 's/d__//' length_taxa | sed 's/;[pcofgs]__/\t/g' > length_taxa2

#get number of time genome appears
parallel -a <(awk '{print $2}' length_taxa2 | sort | uniq) -j 50 -k "grep -c '{}' length_taxa2"> rpoC_count #get cords of all rpoC genes
#get number of contigs
for i in $(awk '{print $2}' length_taxa2 | sort | uniq); do grep $i length_taxa2 | awk '{print $1}' | sort | uniq | wc -l; done > rpoC_contigs
awk '{print $2}' length_taxa2 | sort | uniq > uniq_GCAs
paste -d '\t' uniq_GCAs rpoC_count rpoC_contigs > counts

#number of rpoC on that contig
for i in $(awk '{print $1}' length_taxa2 | sort | uniq); do grep -c $i length_taxa2; done > contig_rpoC
awk '{print $1}' length_taxa2 | sort | uniq > uniq_contigs
paste -d '\t' uniq_contigs contig_rpoC > temp
mv temp contig_rpoC

#combine using python script
python3 taxawcontigs.py #output metadata.tsv
parallel -a <(awk '{print $1}' metadata.tsv) -j 50 -k "grep -wm 1 '{}' ../gtdb_prokka_concat/rpoC.gff " > strand_sense
awk -F '\t' '{print $7}' strand_sense > temp
sed -i '1s/^/strand\n/' temp
mv temp strand_sense
paste -d "\t" metadata.tsv strand_sense > temp
mv temp metadata.tsv
```
Get exact number of how many genomes (both archea and bacteria) had a certain number of rpoC annotations
```sh
#get actual number of genomes with annoteted rpoC (does not include rpoC2)
awk '{print $2}' metadata.tsv | sed '1d' | sort | uniq > gca_rpoc 
#should get 96615 uniq genomes
rm gca_number
for i in $(cat gca_rpoc); do grep -c $i metadata.tsv >> gca_number; done
#combine asscension number with the amount of rpoC annotations 
paste gca_rpoc gca_number > temp
mv temp gca_number
#make file with number or rpoC genes annotated and how many genomes were annotated
echo "num_of_genomes" > num_rpoC # header
wc -l gca_rpoc | awk '{print $1}' >> num_rpoC # total number of genomes with any amount of rpoc annotated
awk '{print $14}' metadata.tsv | sed '1d' | sort -n | uniq > number_annots.txt
for i in $(cat number_annots.txt); do sed 's/\./_/' gca_number| grep -wc $i >> num_rpoC; done #get number of genomes with x amount of rpoc annotated 
#get rownames
echo "num_rpoC_genes" > row_name
echo "any" >> row_name 
for i in $(cat number_annots.txt); do echo "$i" >> row_name; done 
#paste rownames and numbers together
paste row_name num_rpoC > temp
mv temp num_rpoC
cat num_rpoC
num_rpoC_genes	num_of_genomes
# any	96615
# 1	90084
# 2	5532
# 3	424
# 4	169
# 5	180
# 6	149
# 7	14
# 8	41
# 9	11
# 10	7
# 12	2
# 13	1
# 15	1
```
Get exact number of how many genomes (only bacteria) had a certain number of rpoC annotations
```sh
#get actual number of genomes with annoteted rpoC (does not include rpoC2)
grep Bacteria metadata.tsv  | awk '{print $2}' | sort | uniq > gca_rpoc 
#should get 69633 uniq genomes
rm gca_number
for i in $(cat gca_rpoc); do grep -c $i metadata.tsv >> gca_number; done
#combine asscension number with the amount of rpoC annotations 
paste gca_rpoc gca_number > temp
mv temp gca_number
#make file with number or rpoC genes annotated and how many genomes were annotated
echo "num_of_genomes" > num_rpoC # header
wc -l gca_rpoc | awk '{print $1}' >> num_rpoC # total number of genomes with any amount of rpoc annotated
grep Bacteria metadata.tsv | awk '{print $14}' | sort -n | uniq > number_annots.txt
for i in $(cat number_annots.txt); do sed 's/\./_/' gca_number| grep -wc $i >> num_rpoC; done #get number of genomes with x amount of rpoc annotated 
#get rownames
echo "num_rpoC_genes" > row_name
echo "any" >> row_name 
for i in $(cat number_annots.txt); do echo "$i" >> row_name; done 
#paste rownames and numbers together
paste row_name num_rpoC > temp
mv temp num_rpoC
cat num_rpoC
# num_rpoC_genes	num_of_genomes
# any	91304
# 1	87787
# 2	2599
# 3	356
# 4	157
# 5	180
# 6	149
# 7	14
# 8	41
# 9	11
# 10	7
# 12	2
# 13	1
```
## Length distro for single copy
```sh
awk -F '\t' '$14 == 1' metadata.tsv > single_rpoC #search for which genomes have 1 rpoC annotations
```
Make length distro plots in R
```R
library(ggplot2)
distro <- read.table("./single_rpoC",row.names=1)
colnames(distro)[2] <- "genome" #rename column
colnames(distro)[5] <- "length" #rename column
colnames(distro)[6] <- "kingdom" #rename column
colnames(distro)[7] <- "phylum" #rename column
colnames(distro)[8] <- "class" #rename column
colnames(distro)[9] <- "order" #rename column
colnames(distro)[10] <- "family" #rename column
colnames(distro)[11] <- "genus" #rename column
colnames(distro)[12] <- "species" #rename column
colnames(distro)[13] <- "rpoC_gca" #rename column
colnames(distro)[14] <- "contig_number" #rename column
colnames(distro)[15] <- "rpoC_num" #rename column
colnames(distro)[15] <- "strand" #rename column

median(distro$length)
mean(distro$length)

pdf("length_distro.kingdom.pdf", width = 20, height =20)
ggplot(distro, aes(x=length, fill =kingdom)) + 
  geom_histogram()+
  theme_minimal()
dev.off()

pdf("length_distro.genome.pdf", width = 20, height =20)
ggplot(distro, aes(x=length, fill =genome)) + 
  geom_histogram()+
  theme_minimal()
dev.off()

pdf("length_distro.phyla.pdf", width = 20, height =20)
ggplot(distro, aes(x=length, fill =phylum)) + 
  geom_histogram()+
  facet_wrap(~kingdom)+
  theme_minimal()
dev.off()
#only bacteria 
bact.distro <- distro[is.element(distro$kingdom, c("Bacteria")), ]
median(bact.distro$length)
# 4124
sd(bact.distro$length)
# 893.6991
```