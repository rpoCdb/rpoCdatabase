#create list of unique named species
python3 unique_species.py #output unique_species.tsv
# subset unique named species by genus list
python3 sub_genomes.py #output filtered_species.tsv

#how many genomes with that genera present
parallel -a <(cat pom_genera.txt) -j 50 -k "grep '{}' taxa.tsv | wc -l"> genera_count 
paste pom_genera.txt genera_count
# Acidovorax	35
# Actinobaculum	1
# Actinomyces	39
# Aerococcus	22
# Alloscardovia	4
# Anaerococcus	28
# Aquabacterium	47
# Atopobium	2
# Bacteroides	121
# Bifidobacterium	121
# Blautia	116
# Campylobacter	10
# Cloacibacterium	5
# Clostridiales	172
# Coriobacteriaceae	802
# Corynebacterium	157
# Dialister	33
# Diaphorobacter	8
# Facklamia	4
# Faecalibacterium	38
# Finegoldia	6
# Gardnerella	0
# Lachnospiraceae	3094
# Lactobacillaceae	431
# Lactobacillales	1164
# Lactobacillus	55
# Megasphaera	26
# Mobiluncus	6
# Oscillibacter	5
# Parvimonas	5
# Peptoniphilus	17
# Prevotella	582
# Propionibacterium	3
# Propionimicrobium	1
# Ruminococcaceae	994
# Ruminococcus	321
# Saccharofermentans	84
# Sneathia	3
# Staphylococcus	58
# Streptococcus	437
# Subdoligranulum	0
# Varibaculum	5
# Veillonella	117
# Leptotrichia	37
# Escherichia	10

#how many unique species in each genus that fultered for

parallel -a <(cat pom_genera.txt) -j 50 -k "grep -w '{}' filtered_species.tsv | wc -l"> genera_count 
paste pom_genera.txt genera_count
# genus	1
# Acidovorax	8
# Actinobaculum	1
# Actinomyces	27
# Aerococcus	6
# Alloscardovia	4
# Anaerococcus	11
# Aquabacterium	7
# Atopobium	2
# Bacteroides	34
# Bifidobacterium	82
# Blautia	12
# Campylobacter	0
# Cloacibacterium	2
# Clostridiales	69
# Coriobacteriaceae	9
# Corynebacterium	105
# Dialister	4
# Diaphorobacter	1
# Facklamia	4
# Faecalibacterium	5
# Finegoldia	1
# Gardnerella	0
# Lachnospiraceae	226
# Lactobacillaceae	319
# Lactobacillales	608
# Lactobacillus	31
# Megasphaera	4
# Mobiluncus	3
# Oscillibacter	2
# Parvimonas	1
# Peptoniphilus	1
# Prevotella	51
# Propionibacterium	3
# Propionimicrobium	1
# Ruminococcaceae	66
# Ruminococcus	1
# Saccharofermentans	0
# Sneathia	2
# Staphylococcus	51
# Streptococcus	94
# Subdoligranulum	0
# Varibaculum	1
# Veillonella	11
# Leptotrichia	6
# Escherichia	5

#grab the rpoC sequences
# awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' ../gtdb_prokka_concat/all.fna > all.fna
awk -F "\t" '{print $1, $4, $5}' filtered_species.tsv | sed 's/ /\t/g' | sed '1d'> ids.txt
seqtk subseq ../gtdb_prokka_concat/all.fna ids.txt > pom_rpoC.fa
sed '/^.*$/N; />\n/!P; D' pom_rpoC.fa | sed 's/ .*//' > rpoC_pom.fa
awk '{print $1}' ids.txt | while read line; do grep -m 1 -A 1 $line rpoC_pom.fa; done > rpoc_pom.fa
#fixing headers
sed -i 's/:.*//' rpoc_pom.fa
paste -d ' ' <(grep '^>' rpoc_pom.fa) <(awk '{print $13}' filtered_species.tsv | sed '1d') > combined.txt
awk '{print $1"_"$2}' combined.txt > headers.txt

paste -d "\t" <(grep '^>' rpoc_pom.fa) headers.txt | sed 's/>//g '> names.txt
seqkit replace -p "(.+)" -r '{kv}|$1' -k names.txt rpoc_pom.fa |sed 's/|.*//' > rpoc_species.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' rpoc_species.fa > temp
mv temp rpoc_species.fa
#reverse files that are on a negative strand
grep -w "-"  filtered_species.tsv | awk -F "\t" '{print $1}' | while read line; do grep -A 1 $line rpoc_species.fa; done > rpoc_negative.fa
seqtk seq -r rpoc_negative.fa > negative.fa
# bioawk -c fastx '{print ">"$name;print reverse($seq)}' rpoc_negative.fa > reverse.fa
grep -w "+"  filtered_species.tsv | awk -F "\t" '{print $1}' | while read line; do grep -A 1 $line rpoc_species.fa; done > postive.fa
cat postive.fa negative.fa > rpoc_pom.fa
