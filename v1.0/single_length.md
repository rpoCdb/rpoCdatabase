# Get length distrubtion of single copy rpoC genes
## Find how many prokka annotated
```sh
grep -v "gene=rpoC2" rpoC.gff | grep -v "gamma" > non_rpoC2.gff #80363 non-rpoC2 genes, gamma subunits are found in cholorplasts
#get number greater than a certain amount
grep "gene=rpoC_" non_rpoC2.gff | sed 's/.*gene=rpoC_//' | sed 's/;.*//' | sort -n | uniq > number_annots.txt #get annotation number. highest number is 13
```
Get how many had x or above annotated
```sh
grep -c "gene=rpoC" non_rpoC2.gff > rpoc_annots.txt #get all rpoc annotations
grep -c "gene=rpoC;" non_rpoC2.gff >> rpoc_annots.txt #get single copy rpoc annotations
for i in $(cat number_annots.txt); do grep -c "gene=rpoC_$i" non_rpoC2.gff >> rpoc_annots.txt; done #get number of duplicate rpoc annotions
#get rownmames
echo "total_rpoc_annots" > row_name
echo "single_copy_annots" >> row_name
for i in $(cat number_annots.txt); do echo "$i≥rpoc" >> row_name; done

paste row_name rpoc_annots.txt > temp
mv temp rpoc_annots.txt
```
Get exact number of how many genomes had a certain number of rpoC annotations
```sh
#get actual number of genomes with annoteted rpoC (does not include rpoC2)
awk '{print $1}' non_rpoC2.gff | sort | uniq > acces_rpoc 
#should get 74828 uniq seqids numbers
for i in $(cat acces_rpoc); do grep -c $i non_rpoC2.gff >> acces_number; done
#combine asscension number with the amount of rpoC annotations 
paste acces_rpoc acces_number > temp
mv temp acces_number
#make file with number or rpoC genes annotated and how many genomes were annotated
echo "num_of_genomes" > num_rpoC # header
wc -l acces_rpoc | awk '{print $1}' >> num_rpoC # total number of genomes with any amount of rpoc annotated
for i in $(cat number_annots.txt); do sed 's/\./_/' acces_number| grep -wc $i >> num_rpoC; done #get number of genomes with x amount of rpoc annotated 
#get rownames
echo "num_rpoC_genes" > row_name
echo "any" >> row_name 
for i in $(cat number_annots.txt); do echo "$i" >> row_name; done 
#paste rownames and numbers together
paste row_name num_rpoC > temp
mv temp num_rpoC
cat num_rpoC
# num_rpoC_genes	num_of_genomes
# any	74828
# 1	70732
# 2	3516
# 3	208
# 4	96
# 5	121
# 6	110
# 7	37
# 8	5
# 9	3
# 10	0
# 11	0
# 12	0
# 13	0
```
## Get columns with accesion number and coordinates from gff file
```sh
#grep "gene=rpoC;" non_rpoC2.gff | sed 's/\t.*//' > single_acces
awk -F '\t' '$2 == 1' acces_number | awk '{print $1}' > single_acces #search for which accesions have 1 rpoC annotations
#get cords from genomes with just one rpoC
rm cords
for i in $(cat single_acces); do grep -Fw $i non_rpoC2.gff >> cords; done
grep "gene=rpoC;" cords| awk {'print $1, $4, $5'}
wc -l cords #out put should match the next ones output
head -n 2 rpoc_annots.txt | awk '{print $2}' | tail -n 1 #make sure out put matches from cords
awk '{print $1, $4, $5}' cords | sed 's/ /\t/g'> rpoC_cords
#get seqid to accescion number
scp -r ./cords suzanne@hillary.clemson.edu:/home/suzanne/rpoCdb/cords
find /home/allie/rpoCdb/gtdb_version1.0/gtdb_genomes_reps_r214 -type f  -name "*.gff" -exec grep -HF rpoC {} + >rpoC_files 
sed 's/_genomic.*gff:/\t/' rpoC_files | sed 's/.*r214\///' | sed 's/\tprokka.*//' | sed 's/\tProdigal.*//' | sort | uniq  > acces2seqid
scp -r ./acces2seqid suzanne@stella.clemson.edu:~/rpoCdb/rpoC_annotated

#get cordinates and seqid
python3 cord_flip.py #run python script that edits rpoC_cords
sed '1d' flip_cords.tsv | awk '{print $2, $3, $4, $5}' | sed 's/ /\t/g'> temp
mv temp flip_cords.tsv

#join cordinate file to accescion number
python3 seqid2acces.py #outputs cords.tsv
#check for any missing accesion numbers

awk '{print $1}' cords.tsv | sed '1d '> annot_cords
awk '{print $1}' cords >seqids_needed
python3 miss_ids.py
awk '{print $1}' get_seqs.tsv | wc -l # should be empty
```
Get the full taxonomy using asscesnion number
```sh
#get archea taxonomy
wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz
gzip -d ar53_metadata.tsv.gz
awk -F "\t" '{print $1, $17}' ar53_metadata.tsv | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' > accs2taxa
#get bacteria taxonomy
awk -F "\t" '{print $1, $17}' bac120_metadata.tsv | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/.*G/G/' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa
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
awk -F "\t" '{print $55, $17}' missed_taxa | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
python3 chords2taxa.py #taxa.tsv is the output

#get missing taxa second round
grep left_only taxa.tsv | awk '{print $5}' > missing_taxa
wc -l missing_taxa
awk -F "\t" '{print $1, $17}' full_metadata.tsv | sed '1d' | sed 's/[GR][BS]_G/G/' | sed 's/ /\t/' | sed 's/ /_/g' >> accs2taxa
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
## Get length distro for rpoC2
Prepare file for R
```sh
awk '{print $1, $5, $2, $3, $4, $6}' taxa.tsv | sed '1d' | sort | uniq > length_taxa
grep -c ";s__" length_taxa #numbers should match
grep -c "d__" length_taxa
sed 's/d__//' length_taxa | sed 's/;[pcofgs]__/\t/g' > length_taxa2
```
Make length distro plots in R
```R
library(ggplot2)
distro <- read.table("./length_taxa2",row.names=1)
colnames(distro)[4] <- "length" #rename column
colnames(distro)[5] <- "kingdom" #rename column
colnames(distro)[6] <- "phylum" #rename column
colnames(distro)[7] <- "class" #rename column
colnames(distro)[8] <- "order" #rename column
colnames(distro)[9] <- "family" #rename column
colnames(distro)[10] <- "genus" #rename column
colnames(distro)[11] <- "species" #rename column
median(distro$length)
mean(distro$length)

pdf("length_distro.kingdom.pdf", width = 20, height =20)
ggplot(distro, aes(x=length, fill =kingdom)) + 
  geom_histogram()+
  theme_minimal()
dev.off()

pdf("length_distro.phyla.pdf", width = 20, height =20)
ggplot(distro, aes(x=length, fill =phylum)) + 
  geom_histogram()+
  facet_wrap(~kingdom)+
  theme_minimal()
dev.off()
```