# ISC_final
---
title: "ISC Final"
author: "Tiffany Marchell"
date: "March 19, 2016"
output: html_document
---
# Overview 
For my final project I wanted to analyze a raw file of B cell receptor sequences starting with the raw results from the sequencing core. The raw sequence file contains ~200,000 sequences of the B cell receptors from a single patient, collected by Patrick Wilson's lab. These sequences were compiled from patient blood samples at day 7, 14, and 90 post influenze vaccination. 


The day 7 sequences are representative of the receptor repertoire of the plasmablast B cell population--these are the newly expanded B cell clones in response to the challenge with influenza, and should all be relatively influenza-specific.
The day 14 sequences are representative of the early memory B cell repertoire. These cells should be a mixture of memory B cells in response to influenze vaccination, as well as to the patient's history of immunologic challenges. 
The day 90 sequences are representative of the long lived memory B cell repertoire. These cells are a stable population that comprise B cells that once responded to prior immunologic challenges, at any point within this patient's life. 

In the files provided, the file "rawseq.txt" contains all of the raw sequence results. 

The variable domain is also referred to as the V region and is the most important region for binding to antigens. To be specific, variable loops of β-strands, three each on the light (VL) and heavy (VH) chains are responsible for binding to the antigen. I am interested in seeing if there are large differences between useage of specific variable genes in the different populations of B cells. My hypothesis is that the plasmablast (influenza specific) B cell populations will have a bias towards the use of a single V gene or a group of a few V genes that are useful in generating flu-antigen specificity. I would expect that the distribution of V gene usage for the plasmablast population to be different than the 14 or 90 day memory B cell population. 

# Plan of attack: 
###1. Remove duplicate sequences. 
Because of the way receptors were PCR'ed and sequenced, the same B cell receptor may have multiple cDNA copies made. When the cDNA "library" was made by putting Illumina sequencing identifiers on the ends, the same cDNA may have ended up with 2 unique library identifiers. When we then get the sequencing results, this could lead to a false positive overrepresentation of certain V genes so its important to get rid of duplicates. 

###2. Separate individual populations from combined raw sequences file. 

###3. To determine what V genes were being used, run sequences through online IMGT/V-Quest IgBLAST.
Unfortunately, there was no good way to automate this. The immunoglobulin sequence database is not available to call via the NCBI BLAST search function we learned in class. Also, due to the high throughput nature, I had to sign up for a "High Throughput Analysis" account, submit my sequences and recive the alignment results via csv file download via email a few days later. Here is the online basic version: http://www.imgt.org/IMGT_vquest/vquest?livret=0&Option=humanIg 
  
###4. Extract just V genes segment of summary from results of IgBLAST analysis.

###5. Import the Vgenes data for each of the 3 populations, count the useage of speicific V genes across the 3, and print results. 


## Remove Duplicates 

A single python program file "prepseq_for_submit.py" imports the next two functions: 
remove_duplicates()
separatebytype()

The first function is remove duplicates. It is specifically set up so that the rawseq.txt file is opened and it will run and save the output file in the same folder. 

It can be called by simply entering remove_duplicates(). 

```{python}

from Bio import SeqIO

def remove_duplicates():
	
	sequences = {}
	min_length = 5 
	#set min length at 5 just so if there happen to be any "empty" sequences they wont be chosen. any number >0 would work
	with open ('rawseq.txt', "r") as f:
		for seq_record in SeqIO.parse(f, "fasta"):
			sequence = str(seq_record.seq)

		#make sure sequence is larger than minimum length
			if (len(sequence)>= min_length ): 
				if sequence not in sequences:
					sequences[sequence] = seq_record.id
					#if its already there then just continue
				else: 
					continue

		#save the unique sequences to a file w/ same name as 
		#input file, just add nodup_ to the name of the OG file 
		output_file = open("nodup_rawseq.txt", "w")
		for sequence in sequences:
			output_file.write( ">" + sequences[sequence] + "\n" + sequence + "\n")
		output_file.close()

	print ("done.")

```

The output of this function is a new file with "nodup_" added to the title. 
This file "nodup_rawseq.txt" then becomes the input file for the next function. 

## Separate by Type

To separate plasmablast, 14 day memory, and 90 day memory sequences into 3 separate fasta files for submission to the IMGT Vquest IgBLAST. 

Just run function separatebytype() in command line. 

```{python}

import re

def separatebytype():
	input_file = open('nodup_rawseq.txt', 'r')
#these are the files to write the specificied sequences to
	output_handle1 = open("PB_seq.fasta", "w")
	output_handle2 = open("14dMemory_seq.fasta", "w")
	output_handle3 = open("90dMemory_seq.fasta", "w")

	for record in SeqIO.parse(input_file, "fasta"):
	#find and separate all plasmablast, day 14 memory, and day 90 memory seq
		if re.search("PB-HT", record.description):
			SeqIO.write(record, output_handle1, "fasta")

		if re.search("14low-|14hi", record.description):
			SeqIO.write(record, output_handle2, "fasta")

		if re.search("90Hi-|90low-", record.description):
			SeqIO.write(record, output_handle3, "fasta")

	output_handle1.close() 
	output_handle2.close()
	output_handle3.close()
	input_file.close()

	print ("done")
```

The output of this function is 3 files: "PB_seq.fasta", "14dMemory_seq.fasta", and "90dMemory_seq.fasta". 
These were then submitted online to the IMGT/V-QUEST. 


## IgBLAST on IMGT/V-QUEST Database
Because I needed to do sequence alignment on 200,000 sequences to find out what V genes each was using, I needed to sign up for a high throughput version of the tool (http://www.imgt.org/HighV-QUEST/login.action).

I was able to upload all of my data at once as FASTA files for it will run the analysis on their servers. Several days later when I was notified the alignments were complete I then could download my results as a .csv file. 

For some reason my 90d_Memory sequences submission is, unfortunately, still queued. I submitted these 3/16/16 but the website does say it may take several days for the analysis to be finished. Strictly for this 90d_Memory file only (for submitting purposes) I ran the max (50 sequences) submission for the 90 day sequences on the basic web version of V-quest and will use those results for the remainder of the analysis. 

These files are: 
PB_BLAST.csv
14d_BLAST.csv
90d_BLAST.csv

## Extract V gene Info

Each of these files are then run through extract_Vgenescolumn() function to separate out just the column containing the V gene annotation information. 

Import file extract_vgenes.py and run extract_Vgenescolumn() function. 

This function won't need an imput, the files to open are already specified within the code itself. 

```{python}
import csv

def extract_Vgenescolumn():
	
	genes = []
	output_file = open("Vgenes_PB_BLAST.csv", "w")
	values = csv.reader(open("PB_BLAST.csv" , 'r'), delimiter = ',')
	#want just the row w/ the Vgenes
	for row in values: 
		genes.append(row[3])

		#write them to file
	for gene in genes: 
		output_file.write("%s\n" % gene)
	output_file.close()
	
	print ("Plasmablast V genes extracted...")

# repeat for 14 day memory sequences 
	
	genes2 = []
	output_file = open("Vgenes_14d_BLAST.csv", "w")
	values = csv.reader(open("14dBLAST.csv" , 'r'), delimiter = ',')
	#want just the row w/ the Vgenes
	for row in values: 
		genes2.append(row[3])

		#write them to file
	for gene2 in genes2: 
		output_file.write("%s\n" % gene2)
	output_file.close() 
	
	print ("14 day Memory V genes extracted...")

#repeat for 90 day memory sequences 
	genes3 = []
	output_file = open("Vgenes_90d_BLAST.csv", "w")
	values = csv.reader(open("90dBLAST.csv" , 'r'), delimiter = ',')
	#want just the row w/ the Vgenes
	for row in values: 
		genes3.append(row[3])

		#write them to file
	for gene3 in genes3: 
		output_file.write("%s\n" % gene3)
	output_file.close() 
	
	print ("90 Day memory V genes extracted...")
	print ("done.")
```

The output of this function is 3 files:
Vgenes_PB_BLAST.csv
Vgenes_14d_BLAST.csv
Vgenes_90d_BLAST.csv

## Analyze V gene Usage Across Populations

The last function Vgene_count( ), takes the 3 Vgenes_population_BLAST.csv files as input, and analyzes the frequency of use of each V gene.

Import Vgene_counts.py and run function Vgene_count('input_file_1', 'input_file_2', 'input_file_3').

You will need to type in the 3 names of the input files for this one to work correctly. 

```{python}

import re
import csv


def Vgene_count(input_file_1, input_file_2, input_file_3):
	
	with open(input_file_1) as f: 
		my_file1 = f.read()
		vgenes1 = re.findall(r'\nHomsap\sIGHV\d?', my_file1)
		stats1 = {}
		for gene1 in vgenes1:
			if gene1 in stats1:
				stats1[gene1] += 1

			else:
				stats1[gene1] = 1
		print ('-----------------------------------PLASMABLASTS----------------------------------------:')
		
		for gene1 in sorted(stats1, key=stats1.get):
			print ("%s : %d " % (gene1, stats1[gene1]))
			print (" %d %% " % (int((stats1[gene1]) / float(sum(stats1.values())) * int(100))))
		print ('\n')
		


	with open(input_file_2) as f: 
		my_file2 = f.read()
		vgenes2 = re.findall(r'\nHomsap\sIGHV\d?', my_file2)
		stats2 = {}
		for gene2 in vgenes2:
			if gene2 in stats2:
				stats2[gene2] += 1

			else:
				stats2[gene2] = 1
		print ('-----------------------------------14 DAY MEMORY---------------------------------------:')
		
		for gene2 in sorted(stats2, key=stats2.get):
			print ("%s : %d " % (gene2, stats2[gene2]))
			print (" %d %% " % (int((stats2[gene2]) / float(sum(stats2.values())) * int(100))))
			
		print ('\n')


	with open(input_file_3) as f: 
		my_file3 = f.read()
		vgenes3 = re.findall(r'\nHomsap\sIGHV\d?', my_file3)
		stats3 = {}
		for gene3 in vgenes3:
			if gene3 in stats3:
				stats3[gene3] += 1

			else:
				stats3[gene3] = 1
		print ('-----------------------------------90 DAY MEMORY---------------------------------------:')
		for gene3 in sorted(stats3, key=stats3.get):
			print ("%s : %d " % (gene3, stats3[gene3]))
			print (" %d %% " % (int((stats3[gene3]) / float(sum(stats3.values())) * int(100))))

		print ('\n')


#save dictionaries w/ genes and counts into csv files to analyze
	with open('PBresultsdict.csv', 'w') as r:
		w = csv.DictWriter(r, stats1.keys())
		w.writeheader()
		w.writerow(stats1)

	with open('14day_resultsdict.csv', 'w') as r: 
		w = csv.DictWriter(r, stats2.keys())
		w.writeheader()
		w.writerow(stats2)

	with open('90day_resultsdict.csv', 'w') as r:
		w = csv.DictWriter(r, stats3.keys())
		w.writeheader()
		w.writerow(stats3)

```


The printout at the end is the V genes found, how many times that V gene was used, and percentage of use for that V gene out of the total V gene sequences. 

This file also outputs the 3 dictionaries of the V gene counts to csv files for further use in analysis if wanted. 

#Conclusion
From this analysis it looks as though V3 gene is the most abundantly used V gene across all 3 B cell populations. It does look like the 90 Day memory population does use more of a variety of V genes, i.e. IGHV5, IGHV6, IGHV7. Although I only ran 50 sequences for this population, I would expect similar results from running all 90 day memory sequences. 



