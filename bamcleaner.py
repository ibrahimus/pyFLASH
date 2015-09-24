#!/usr/bin/python
#todo: work with paired-end data
from __future__ import division

__author__		= "Ibrahim Ilik"
__copyright__	= "Copyleft 2015"
__version__		= "0.0.2"
__credits__		= ["Ibrahim Ilik"]
__maintainer__	= "None"
__email__		= "ibrahimus@gmail.com"

import Levenshtein
import time
import argparse
import hashlib
import subprocess
import os
import random
import multiprocessing as mp
import glob
import sys

def DupRemover(inputbam, output, bedout, chunk, quality, multi):
	start = time.time()
	#pre-processing. Generating unique names for temporary files, in case multiple instances of the script is run at the same time.

	mm = hashlib.sha1()
	mm.update(inputbam.split("/")[-1]+str(random.random()))
	tempbed = ".temp_bed-"+mm.hexdigest()
	tempbed0 = ".0temp_bed-"+mm.hexdigest()
	tempsam = ".temp_sam-"+mm.hexdigest()
	ticker = False
	ticker2 = True
#	logfile = "cleaned/"+inputbam.split(".bam")[0]+"_cleaned.log"
	logfile = output.split(".bam")[0]+".log"

	if multi:
		with open(logfile, 'a') as ff:
			msg = "Working on file %s\n" %inputbam
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
	else:
		printStatus("Working on the file: %s" %inputbam)
	
	output = output.split(".bam")[0]+".sam"

	# check if the input file is from bbmap or bowtie

	cmd0 = "bedtools bamtobed -i %s | head > %s" %(inputbam, tempbed0)
	subprocess.check_call(cmd0, shell=True)

	xx = open(tempbed0, 'r')

	# This routine tries to understand if the input file came from bowtie2 or bbmap (because bbmap generate some funny .sam files). 
	# It may not work for other aligners. Although I guess it might also work.

	if len(xx.next().strip().split(" ")) == 1:
		cmd = "bedtools bamtobed -i %s | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(inputbam, tempbed)
	elif len(xx.next().strip().split(" ")) == 2:
		ticker2 = False # marks  as "bbmap" output
		cmd = "bedtools bamtobed -i %s | awk '{print $1,$2,$3,$4,$6,$7}' | tr [:blank:] \\\\t | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(inputbam, tempbed)
	elif len(xx.next().strip().split(" ")) == 3:
		ticker2 = False#  marks  as "bbmap" output
		cmd = "bedtools bamtobed -i %s | awk '{print $1,$3,$4,$5,$7,$8}' | tr [:blank:] \\\\t | sed 's/:/ /g' | sort -k1,1 -k2,2n -k3,3n -k11,11 -k12,12n -k13,13 > %s" %(inputbam, tempbed)
	else:
		print "Unrecognized .bam formatting."
		sys.exit(1)
	xx.close()
	os.remove(tempbed0)

	if multi:
		with open(logfile, 'a') as ff:
			msg = "Generating some temp files.\n"
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
	else:
		printStatus("Generating some temp files.")
	cmd2 = "samtools view -h %s > %s" %(inputbam, tempsam)
	if multi:
		with open(logfile, 'a') as ff:
			msg = "Running command: \n\t%s\n"%cmd
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))
	else:
		printStatus("Running command: \n\t%s"%cmd) 
	subprocess.check_call(cmd, shell=True)
	if multi:
		with open(logfile, 'a') as ff:
			msg = "Running command: \n\t%s"%cmd2
			ff.write("%s: %s\n" % (time.strftime("%X %x"), msg))	
	else:
		printStatus("Running command: \n\t%s"%cmd2) 
	subprocess.check_call(cmd2, shell=True)

	if multi:
		with open(logfile, 'a') as ff:
			msg = "Done. Starting to work on the temp files.\n"
			ff.write("%s: %s" % (time.strftime("%X %x"), msg))	
	else:
		printStatus("Done. Starting to work on the temp files.")

	aa = open(tempbed, 'r')

	seq_set = set()
	name_set = set()
	sam_set = set()
	counter = 0 # reads that are kept
	counter2 = 2 # total number of reads
	counter3 = 0 # reads that are below the quality threshold. 
	multiplier = 1

	# Initiate the whole thing with the first two lines.
	line1 = aa.next()
	line2 = aa.next()

	line1_bcode = line1.split(" ")[-1].split("\t")[0]
	line2_bcode = line2.split(" ")[-1].split("\t")[0]
	line1_qual = int(line1.split("\t")[-2])
	line1_strand = line1.strip().split("\t")[-1]
	line2_strand = line2.strip().split("\t")[-1]
	line1_begin = line1.split("\t")[0] + "," + line1.split("\t")[1]
	line1_end = line1.split("\t")[0] + "," + line1.split("\t")[2]
	line2_begin = line2.split("\t")[0] + "," + line2.split("\t")[1]
	line2_end = line2.split("\t")[0] + "," + line2.split("\t")[2]
	levent = Levenshtein.hamming(line1_bcode, line2_bcode)
	
	#comparing consecutive lines with each other. 

	if line1_strand == line2_strand:
		if line1_strand == "+":
			if line1_begin != line2_begin:
				if bedout:
					if line1_qual >= quality:
						seq_set.add(":".join(line1.split(" ")))
				if line1_qual >= quality:	
					name_set.add(":".join(line1.split(" ")).split("\t")[3])
					counter += 1
				else:
					counter3 +=1
			else:	
				if levent > 1:
					if bedout:
						if line1_qual >= quality:
							seq_set.add(":".join(line1.split(" ")))
					if line1_qual >= quality:
						name_set.add(":".join(line1.split(" ")).split("\t")[3])
						counter += 1
					else:
						counter3 += 1

		elif line1_strand =="-":
			if line1_end != line2_end:
				if line1_qual >= quality:
					if bedout:
						seq_set.add(":".join(line1.split(" ")))
				if line1_qual >= quality:
					name_set.add(":".join(line1.split(" ")).split("\t")[3])
					counter += 1
				else:
					counter3 += 1
			else:	
				if levent > 1:
					if line1_qual >= quality:
						if bedout:
							seq_set.add(":".join(line1.split(" ")))
					if line1_qual >= quality:
						name_set.add(":".join(line1.split(" ")).split("\t")[3])
						counter += 1
					else:
						counter3 += 1
		else:
			printStatus("The strand is neither + or - Something must be wrong, exiting")
			aa.close()
			os.remove(tempbed)
			os.remove(tempsam)
			sys.exit(1)
	else:		
		if line1_qual >= quality:
			if bedout:
				seq_set.add(":".join(line1.split(" ")))
		if line1_qual >= quality:
			name_set.add(":".join(line1.split(" ")).split("\t")[3])
			counter += 1
		else:
			counter3 += 1

	#switching lines, so we can iteratively work with consecutive lines
		
	line1 = line2

	while True:
		try:
			line1_bcode = line1.split(" ")[-1].split("\t")[0]
			line1_qual = int(line1.split("\t")[-2])
			line2 = aa.next()
			counter2 += 1 #this counter holds the total number of lines in the file. 
			line2_bcode = line2.split(" ")[-1].split("\t")[0]
			line1_strand = line1.strip().split("\t")[-1]
			line2_strand = line2.strip().split("\t")[-1]
			line1_begin = line1.split("\t")[0] + "," + line1.split("\t")[1]
			line1_end = line1.split("\t")[0] + "," + line1.split("\t")[2]
			line2_begin = line2.split("\t")[0] + "," + line2.split("\t")[1]
			line2_end = line2.split("\t")[0] + "," + line2.split("\t")[2]			
			levent = Levenshtein.hamming(line1_bcode, line2_bcode)
			
			if line1_strand == line2_strand:
				if line1_strand == "+":
					if line1_begin != line2_begin:
						if line1_qual >= quality:
							if bedout:
								seq_set.add(":".join(line1.split(" ")))
							name_set.add(":".join(line1.split(" ")).split("\t")[3])
							counter += 1
						else:
							counter3 += 1
					else:	
						if levent > 1:
							if line1_qual >= quality:
								if bedout:
									seq_set.add(":".join(line1.split(" ")))
								name_set.add(":".join(line1.split(" ")).split("\t")[3])
								counter += 1
							else:
								counter3 += 1

				elif line1_strand =="-":
					if line1_end != line2_end:
						if line1_qual >= quality:
							if bedout:
								seq_set.add(":".join(line1.split(" ")))
							name_set.add(":".join(line1.split(" ")).split("\t")[3])
							counter += 1
						else:
							counter3 += 1
					else:	
						if levent > 1:
							if line1_qual >= quality:
								if bedout:
									seq_set.add(":".join(line1.split(" ")))
								name_set.add(":".join(line1.split(" ")).split("\t")[3])
								counter += 1
							else:
								counter3 += 1
				else:
					printStatus("The strand is neither + or - Something must be wrong, exiting.")
					aa.close()
					os.remove(tempbed)
					os.remove(tempsam)
					sys.exit(1)
			else:		
				if line1_qual >= quality:
					if bedout:
						seq_set.add(":".join(line1.split(" ")))
				if line1_qual >= quality:
					name_set.add(":".join(line1.split(" ")).split("\t")[3])
					counter += 1
				else:
					counter3 += 1

			line1 = line2

			# Every 100000 unique reads (or multiples of 100000 as determined by the chunk variable), the set that holds unique reads are flushed. 

			if bedout:
				if counter == multiplier * chunk:
					for i in seq_set:
						with open(bedout, 'a') as ff:
							ff.write(i)
					multiplier += 1
					seq_set = set()

		except StopIteration:
			if bedout:
				for i in seq_set:
					with open(bedout, 'a') as ff:
						ff.write(i)
			if multi:
				with open(logfile, 'a') as ff:
					msg ="%s of %s reads reads were removed. %s reads were removed because they were below the set quality threshold: Q%s.\n" %(counter2-counter, counter2, counter3, quality)
					msg2 = "%%%s of the reads reads were removed. Of the reads thet were removed %%%s were below the set quality threshold\n" %(round(((counter2-counter)/counter2)*100), round((counter3/(counter2-counter))*100))
					ff.write("%s: %s" % (time.strftime("%X %x"), msg))
					ff.write("%s: %s" % (time.strftime("%X %x"), msg2))
			else:
				printStatus(" %s%% of the reads were removed, they were either PCR-duplicates or were below the quality threshold." %(round((1-(counter/(counter+counter2)))*100)))
				printStatus("%s of %s reads reads were removed. %s of the reads were removed because they were below the set quality threshold: Q%s.\n" %(counter2-counter, counter2, counter3, quality))
				printStatus("%%%s of the reads reads were removed. Of the reads thet were removed %%%s were below the set quality threshold" %(round(((counter2-counter)/counter2)*100), round((counter3/(counter2-counter))*100)))
			break
	aa.close()

	bb = open(tempsam, 'r')
	list_of_stuff = ['@SQ', '@PG', '@HD']
	counter3 = 0
	multiplier2 = 1
	if not multi:
		printStatus("Generating the .sam file now.")

	# Capturing the header from the .bam file.

	while True:
		line = bb.next()
		if line.split("\t")[0] in list_of_stuff:
			with open(output, 'a') as ff:
				ff.write(line)
		else:
			if line.split("\t")[0] in name_set:
				sam_set.add(line)
				counter3 += 1
			break


	while True:
		try:
			sam_line = bb.next()
			if ticker2: # If input was from bowtie
				if sam_line.split("\t")[0] in name_set:
					sam_set.add(sam_line)
					counter3 += 1 #this is a useless counter. Should be the same as counter
			else: #if input was from BBMap
				if sam_line.split(" ")[0] in name_set:
					sam_set.add(sam_line.split(" ")[0]+"\t"+"\t".join(sam_line.split("\t")[1:]))
					counter3 += 1 #this is a useless counter. Should be the same as counter				

			if counter3 == multiplier2*chunk:
				for i in sam_set:
					with open(output, 'a') as ff:
						ff.write(i)
				multiplier2 += 1
				sam_set = set()

		except StopIteration:
			for i in sam_set:
				with open(output, 'a') as ff:
					ff.write(i)
			sam_set = set()
			break
	bb.close()

#cleaning up

	os.remove(tempbed)
	os.remove(tempsam)

	if not multi:
		printStatus("Generating a bam file now.")
	cmd_bam = "samtools view -Sb %s | samtools sort - %s" %(output, output.split(".sam")[0])
	subprocess.check_call(cmd_bam, shell=True)
	if not multi:
		printStatus("Generating the index for that bam.")
	cmd_sam = "samtools index %s" %output.split(".sam")[0]+".bam"
	subprocess.check_call(cmd_sam, shell=True)
	os.remove(output)

	if multi:
		end = time.time()
		elapsed = round(end-start)
		with open(logfile, 'a') as ff:
				msg ="Done! Elapsed time: %s seconds.\n" %elapsed
				ff.write("%s: %s" % (time.strftime("%X %x"), msg))

def main():
	parser = argparse.ArgumentParser(description="Remove PCR-duplicates by making use of random barcodes in uvCLAP data. \
		The header of your fastq file should look like this: @39V34V1:117:H9FK3ADXX:1:2205:20357:12179:AACCCGCCAA 1:N:0:GGTAGC where AACCCGCCAA is the random barcode (of any length)")
	parser.add_argument('-i', '--input', help='The input .bam file. Doesn\'t have to be ordered.', required=True)
	parser.add_argument('-r', '--reverse-reads', help='For paired-end data. (not done yet')
	parser.add_argument('-k', '--keep-bed', help='Keep the cleaned .bed file? If so where?')
	parser.add_argument('-o', '--output', help='Output bam file.', required=True)
	parser.add_argument('-c', '--chunk', type=int, help='Chunk size. Number of reads kept in memory before written to disk. This value will be multiplied by 100000. Default = 1', default = 1)
	parser.add_argument('-m', '--min-quality', type=int, help='Minimum mapping quality. Default is 10', default=10)
	args = parser.parse_args()

	
	if args.input and args.output == 'all':
		print "Working on all bam files in this folder. You can find the the cleaned .bam, the index and the log file under the 'cleaned' folder."
		global multi
		multi = True
		if not os.path.exists("cleaned"):
				os.makedirs("cleaned")
		files = glob.glob("*.bam")
		jobs = []
		for i in files:
			p = mp.Process(target=DupRemover, args=(i, "cleaned/"+i.split(".bam")[0]+"_cleaned.bam", args.keep_bed, args.chunk*100000, args.min_quality, multi))
			jobs.append(p)
			p.start()
	elif args.input and args.output == 'all':
		print "Both -i and -o should be 'all' for this to work."
		sys.exit(1)
	else:
		multi = False
		if args.output.split(".")[-1] != 'bam':
			print "The output has to be .bam!"
			print args.output
			print args.output.split(".")[-1]
			sys.exit(1)
		DupRemover(inputbam=args.input,  output=args.output, bedout=args.keep_bed, chunk=args.chunk, quality=args.min_quality, multi=multi)

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

if __name__ == "__main__":
	multi = False
	start1 = time.time()
	main()
	if not multi:
		end = time.time()
		elapsed = round(end-start1)
		printStatus("Done! Elapsed time: %s seconds." %elapsed)