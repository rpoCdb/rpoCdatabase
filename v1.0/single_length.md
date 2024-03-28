# Get length distrubtion of single copy rpoC genes
```sh
grep -v "gene=rpoC2" rpoC.gff | grep -v "gamma" > non_rpoC2.gff #80363 non-rpoC2 genes, gamma subunits are found in cholorplasts
#get number greater than a certain amount
grep "gene=rpoC_" non_rpoC2.gff | sed 's/.*gene=rpoC_//' | sed 's/;.*//' | sort -n | uniq > number_annots.txt #get annotation number for each contig. highest number is 13
```
## Get columns with accesion number and coordinates from gff file
```sh
#grep "gene=rpoC;" non_rpoC2.gff | sed 's/\t.*//' > single_acces
awk '{print $1}' non_rpoC2.gff | sort | uniq > rpoC_seqid #get just seqid
#get cords from genomes with just one rpoC
rm cords
for i in $(cat rpoC_seqid); do grep -Fw $i non_rpoC2.gff >> cords; done #get cords of all rpoC genes
wc -l non_rpoC2.gff #80363
wc -l cords #80363
awk '{print $1, $4, $5}' cords | sed 's/ /\t/g'> rpoC_cords
#get seqid to accescion number
scp -r ./cords suzanne@hillary.clemson.edu:/home/suzanne/rpoCdb/cords
find . -type f -name "*.ffn" | wc -l #85099
find /home/allie/rpoCdb/gtdb_version1.0/gtdb_genomes_reps_r214 -type f  -name "*.gff" -exec grep -HF rpoC {} + >rpoC_files 
sed 's/_genomic.*gff:/\t/' rpoC_files | sed 's/.*r214\///' | sed 's/\tprokka.*//' | sed 's/\tProdigal.*//' | sort | uniq  > acces2seqid
scp -r ./acces2seqid suzanne@stella.clemson.edu:~/rpoCdb/rpoC_annotated
wc -l acces2seqid #75259
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
awk -F "\t" '{print $1, $46, $17}' ar53_metadata.tsv | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' > accs2taxa
#get bacteria taxonomy
awk -F "\t" '{print $1, $46, $17}' bac120_metadata.tsv | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa #should be missing but had 24793
cat ar53_metadata.tsv bac120_metadata.tsv > full_metadata.tsv
n=0
maxjobs=5 #can increase this number if you are feeling brave
rm missed_taxa
for i in $(cat missing_taxa); do
    grep -wm 1 $i full_metadata.tsv >> missed_taxa &
    # limit jobs
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
awk -F "\t" '{print $55, $46, $17}' missed_taxa | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa second round
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa
awk -F "\t" '{print $1, $46, $17}' full_metadata.tsv | sed '1d' | sed 's/[GR][BS]_G/G/' | sed 's/Complete Genome/Complete_genome/' | sed 's/ /\t/' | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa #should be 0

#for i in $(awk '{print $4}' cords.tsv | sed '1d'); do grep -w $i bac120_metadata.tsv >> temp; done
#for i in $(awk '{print $1}' cords); do grep -m 1 $i bac120_metadata.tsv || echo $i "no taxonomy" ; done > bac_rpoc #use the || pipe
#grep "no taxonomy" bac_rpoc # should return empty
#use at your own peril
#for i in $(awk '{print $1}' ~/rpoCdb/cords); do find /home/allie/rpoCdb/gtdb_version1.0/gtdb_genomes_reps_r214 -type f  -name "*.gff" -exec grep -HFm 1 $i {} + >> ~/rpoCdb/accescion.txt; done
#wait
#ps -ef | grep "grep" | awk '{print $2}' | head -n 40 > kil2
#cat kil2 | while read line; do kill $line; done
```
## Get metadata file ready
Prepare file for R
```sh
awk '{print $1, $5, $6, $2, $3, $4, $7}' taxa.tsv | sed '1d' | sed 's/ /\t/g' | sort | uniq > length_taxa
grep -c ";s__" length_taxa #numbers should match
grep -c "d__" length_taxa
sed 's/d__//' length_taxa | sed 's/;[pcofgs]__/\t/g' > length_taxa2

#get number of time genome appears
for i in $(awk '{print $2}' length_taxa2 | sort | uniq); do grep -c $i length_taxa2; done > rpoC_count
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
```
Get exact number of how many genomes (both archea and bacteria) had a certain number of rpoC annotations
```sh
#get actual number of genomes with annoteted rpoC (does not include rpoC2)
awk '{print $2}' metadata.tsv | sed '1d' | sort | uniq > gca_rpoc 
#should get 73656 uniq genomes
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
# num_rpoC_genes  num_of_genomes
# any 73656
# 1 68628
# 2 4288
# 3 331
# 4 115
# 5 130
# 6 114
# 7 38
# 8 7
# 9 3
# 10  1
# 13  1
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
# any	69633
# 1	66982
# 2	1978
# 3	273
# 4	107
# 5	130
# 6	114
# 7	38
# 8	7
# 9	3
# 10	1
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
#4127
sd(bact.distro$length)
#857.8345
```
