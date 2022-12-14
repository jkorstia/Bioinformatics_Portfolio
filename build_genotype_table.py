# nway processing
# Jenny Korstian
# August 25, 2020

#############################################
# This script was written during testing and validation of a new bioinformatic tool "nway" 
# briefly, the "nway" tool uses pairwise genome alignments to determine presence and absence of a TE insertion across species
# Data must be processed through the "nway" tool to produce a .csv file whie is then run through nway_processing.py before this script can be used.
# the nway_processing.py script extracts the sequence for each putative TE insertion and blasts against a library of known TE consensi.

# This script uses information from the blast results (generated by the nway_processing.py script) 
# to create a table of independently determined genotypes that can be compared to those produced by the "nway" tool to estimate error rates of the new tool


# USAGE STATEMENT:
# build_genotype_table.py -s [full path to a list of fasta files to process]

# This script creates a table of binary genotypes for each species
# iterates through each fasta file in the input file (-s option) and appends one new row to the genotype table for each TE insertion
# Also writes progress and summary to stdout

# 3 output files are created: 
# genotype_table.csv
# 	1= TE insertion present in that species, 0= TE insertion absent in that species
# 	'?' genotype assigned when element is too short, a bad match %, or the blast TE type doesnt match the TE type from the nway results
# 	see nway_processing.py for description of how those three groups were determined
# keepers.csv (list of all TE loci that survived all filtering steps)
# losers.csv (list of all TE loci that FAILED filtering steps)





import argparse
import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import csv

##############################################
################## Arguments #################
##############################################

parser = argparse.ArgumentParser(description="Will generate a summary table for each locus after processing with the nway_processing.py script", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s', '--samples', type=str, help='Full path to a list of fasta files to be processed, one file per line.', required = True)

args = parser.parse_args()
LIST_PATH = args.samples

#############################################

# Import list of samples to process
fastalist = []
with open(LIST_PATH, 'r') as f:
	fastalist=f.read().splitlines()
#print(fastalist)
#make empty dataframe to store results table in
genotypes_table=pd.DataFrame()
binary=pd.DataFrame()

########################  Build genotypes table  #######################

# work through each fasta in the user provided list of samples to process & create genotypes table locus by locus
for fasta in fastalist:
	#strip ".fa" from the end of the fasta filename
	locus=fasta.split(".")
	locus=locus[0]
	#create an empty dataframe to store the genotypes for a single locus during the loop.  will be rewritten/emptied at the beginning of each loop.
	g=pd.DataFrame()
	s=pd.DataFrame()
	#create empty dataframe to store all genotype results in
	for record in SeqIO.parse(fasta, "fasta"): 
	#get header from each sequence in the fasta file
		header=record.id 
		 #reformat so that * can be used to separate out chunks
		header = re.sub('_Nway', '*Nway', header) 
		header = re.sub('_Blast', '*Blast', header) 
		 #split header on *'s -> this is a list now
		header = header.split("*")
		g['locus']=[locus]
		s['locus']=[locus]
		#get genotypes and load into into genotypes_table dataframe
		if len(header) <= 1: 
			#then it is the element consensus sequence
			# store the element name as a list
			element=header[0]
			#create a column name as a list
			col=["element"]
			#create a dataframe from with the element (.T transposes rows and columns)
			g['element']=(element)
		else: 
			# use components of split list to build column headers
			col=[header[0]+"_Nway", header[0]+"_Blast"]
			build=header[1:]
			#create a column for the nway genotype & populate with value
			g[col[0]]=(build[0])
			#create a column for the blast genotype & populate with value
			g[col[1]]=(build[1])
			#determine simplified genotype
			blast=build[1].split("_")[0]
			detail=build[1].split("_")[1]
			simple=header[0]+"_simple"
			#print(simple)
			# create a column for a simplified genotype	
			if blast == "Blast0":
				s[simple]=["-"]
			elif blast == "Blast1":
				s[simple]=["+"]
			elif blast == "BlastN":
				if detail == "wrongTE":
					s[simple]=["?"]
				if detail == "BadMatch":
					s[simple]=["?"]
				if detail == "TooShort":
					s[simple]=["?"]
			else:
				s[simple]=["error"]
			#print(g)
	#add the g (all of the entries for that locus) to the genotypes table
	genotypes_table=genotypes_table.append([g])
	binary=binary.append([s])
number_of_species=len(s.columns)-1
bad=int(0.2*number_of_species)
binary=binary.set_index(keys='locus')
pivot=binary.apply(pd.Series.value_counts, axis=1).fillna(0)
#make empty lists to store the keepers and loci to toss.
keep=[]
toss=[]
#questionable=[]
verdict=pd.DataFrame()

for index, row in pivot.iterrows():
	why=pd.DataFrame()
	if row['+'] == number_of_species:
		#print("monomorphic")
		toss.append(index)
		why['reason']=["insertion present in all species"]
	elif row['+'] == [0]:
		#print("not present in at least one species")
		toss.append(index)
		why['reason']=["not present in at least one species"]
	elif row['-'] == [0]:
		#print("not absent in at least one species")
		toss.append(index)
		why['reason']=["not absent in at least one species"]
	elif row['?'] != [0]:
		#print("missing too many loci.  NO missing/undetermined genotypes permitted.")
		toss.append(index)
		why['reason']=["Questionable. NO missing/undetermined genotypes permitted."]
	else:
		#print("survived filtering")
		keep.append(index)
		why['reason']=["survived filtering"]
	verdict=verdict.append([why])
df=verdict['reason'].value_counts()
		
#convert keep & toss lists to pandas dataframes because I dont wanna figure out another way to write them out to a csv file.
keepers=pd.DataFrame(keep)
losers=pd.DataFrame(toss)
#question=pd.DataFrame(questionable)

#write out complete genotypes table to a file
genotypes_table.to_csv("genotype_table.csv", index=False)

#write out keepers list to a file
keepers.to_csv("keepers.csv", index=False, header=False)
#write out losers list to a file
losers.to_csv("losers.csv", index=False, header=False)
#write out questionables list to a file
#question.to_csv("questionable_needs_validation.csv", index=False, header=False)


print("-----------------------")
print("Summary of final decisions below:")
print(df)
print("\nThis run allowed 0 species with undetermined genotypes.  \n    This run classified BlastN_wrongTE genotypes as ?'s, BadMatches as ?'s and shorties as ?'s")

print("-----------------------")
print("finished.")



