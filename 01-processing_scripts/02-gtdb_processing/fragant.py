#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
import os
import gzip
import pandas as pd

gb_file = "all.gbf"
geneName = "rpoC"
basename = "rpoC"

######### MULTI COPY PROCESSING
multi = []
print("Starting to process multi copy genes")
for gb_record in SeqIO.parse(gb_file, "genbank"):
	for feature in gb_record.features:
		if feature.type == "CDS" and 'gene' in feature.qualifiers:
			if geneName in feature.qualifiers['gene'][0] and "_" in feature.qualifiers['gene'][0]:
				multi += [[gb_record.name, \
				feature.qualifiers['gene'][0], \
				feature.location.start, \
				feature.location.end, \
				len(gb_record.seq[feature.location.start:feature.location.end])]]

###### convert to pandas dataframe
df = pd.DataFrame(multi, columns = ['ContigID', 'annotated_as', 'start', 'end', 'length'])
# group into length and coordinates for each contig
df1 = df.groupby("ContigID")["length"].sum().reset_index().merge(df.groupby("ContigID")["start"].min().reset_index())
df = df1.merge(df.groupby("ContigID")["end"].max().reset_index())

# get actual length from coordinates
df.insert(4, "coordinate_len", df['end'] - df['start'])
df.to_csv("multicopy_list.csv") 

print("Multi copy genes identified")

######### SINGLE COPY PROCESSING
single = []
print("Starting to process single copy genes")
# first pull all single copy genes into a single fasta file
for gb_record in SeqIO.parse(gb_file, "genbank"):
	for feature in gb_record.features:
		if feature.type == "CDS" and 'gene' in feature.qualifiers:
			if geneName in feature.qualifiers['gene'][0] and "_" not in feature.qualifiers['gene'][0]:
				single += [[gb_record.name, \
				feature.qualifiers['gene'][0], \
				feature.location.start, \
				feature.location.end, \
				len(gb_record.seq[feature.location.start:feature.location.end])]]


###### convert to pandas dataframe
df = pd.DataFrame(single, columns = ['ContigID', 'annotated_as', 'start', 'end', 'length'])
df.to_csv("singlecopy_list.csv") 

print("Single copy genes identified")