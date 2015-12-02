#!/home/semplicio/virtualenv/bin/python -u
# change the path to your python interpreter

# to do:
# mapping of spike-ins --done
# implement bbduk
# put bbmerge file in a common folder (structure changes)
# put bad barcode logs in a common folder (structure changes maybe)
# i think bbmap only is missing -implement


from __future__ import division
__author__		= "Ibrahim Ilik"
__copyright__	= "Copyright 2015"
__version__		= "0.0.1"
__credits__		= ["Ibrahim Ilik"]
__maintainer__	= "None"
__email__		= "ibrahimus@gmail.com"
__status__		= "Production"

# Some parameters for BBMap and BBMerge. usejni=t only after compiling the C code!
BBMAP_PATH = "/data/akhtar/group/Ibrahim_data_2014/__BBMap__/bbmap"
BBMAP_OPTIONS = "usejni=f vslow=t maxindel=100k sam=1.3"

import re
import string
import time
import sys
import time
import gzip
import Levenshtein
import argparse
import os
import subprocess
import datetime
import glob
import multiprocessing as mp
import itertools
#import numpy as np
#import matplotlib.pyplot as plt
try:
	import bamcleaner as sbp
except ImportError:
	print "bamcleaner.py not found. Stopping."
	sys.exit(1)

current_dir = os.getcwd()+"/"

def processFastqFiles(forward, reverse, barcodes, file_type, bbmerge, ignore, chunk, mapping, gzipp, reps_merged, script_process, gzipafter, user_barcode, hamming_distance, processor, aligner, spikes_mapper):
	# creating folders, extracting and defining barcodes ##
			
	folders = ["barcode_sorted", "final_files", "final_files/bams", "final_files/bam_coverage", "final_files/bam_coverage/bedGraphs", "final_files/bam_coverage/bigwigs", "final_files/xlinks",
				"final_files/xlinks/bed", "final_files/xlinks/bedGraph_coverage", "final_files/xlinks/bigwig_coverage", "barcode_sorted/bowtie2", "barcode_sorted/bowtie2/log",
				"barcode_sorted/bowtie2/cleaned", "barcode_sorted/bowtie2/unmapped", "barcode_sorted/hisat", "barcode_sorted/hisat/log",
				"barcode_sorted/hisat/unmapped", "barcode_sorted/hisat/cleaned", "barcode_sorted/bbmap", "barcode_sorted/bbmap/log", "barcode_sorted/bbmap/cleaned", 
				"final_files/bams/replicates_merged", "final_files/bam_coverage/bedGraphs/replicates_merged", "final_files/bam_coverage/bigwigs/replicates_merged",
				"final_files/xlinks/replicates_merged", "final_files/xlinks/replicates_merged/bed", "final_files/xlinks/replicates_merged/bedGraph", "final_files/xlinks/replicates_merged/bigwig",
				"spikes", "spikes/cleaned"]

	for i in folders:
		if not os.path.exists(i):
			os.makedirs(i)

	index_barcode_names, index_barcodes = extract_barcodes(barcodes)

	index_uniq = {x for x in index_barcodes}

	if len(index_uniq) != len(index_barcodes):
		printStatus("Something is wrong with your barcodes! Are they all uniqe?")
		sys.exit(1)


	# Fill these in as reqired. Should also change the available options down at the main() function, --mapping option.

	indices = {'genome':['path/to/bowtie2_index','bbmap_index_designator_as_build_number','path/to/chrom.sizes', 'path/to/hisat_index'],
				'hg19':['/data/akhtar/group/ibonomics/hg19', 1, '/data/akhtar/group/ibonomics/hg19.chrom.sizes'],
				'dm3':['/data/akhtar/group/ibonomics/dm3', 2, '/data/akhtar/group/ibonomics/dm3.chrom.sizes', '/data/akhtar/group/ibonomics/hisat_index/hisat2/dm3/dm3'], 
				'dm6':['/data/akhtar/group/ibonomics/dm6/dm6', 6, '/data/akhtar/group/ibonomics/dm6/dm6.chrom.sizes', '/data/repository/organisms/dm6_ensembl/HISAT2Index/genome'],
				'mm9':['/data/akhtar/group/ibonomics/mm9/mm9', 3, '/data/akhtar/group/ibonomics/mm9/mm9.chrom.sizes'],
				'mm10':['/data/akhtar/group/ibonomics/mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome', 4, '/data/akhtar/group/ibonomics/mm10.chrom.sizes'],
				'hg19_simple':['/data/akhtar/group/ibonomics/hg19_simple/hg19_simple', 5, '/data/akhtar/group/ibonomics/hg19_simple/hg19_simple.chrom.sizes'],
				'hg38':['/data/akhtar/group/ibonomics/gencodegenes/GRCh38.p3/fasta/GRCh38.p3.genome', 7, '/data/akhtar/group/ibonomics/gencodegenes/GRCh38.p3/fasta/GRCh38.p3.genome.chrom.sizes']}
	
	spike_index = {'spike':'path/to/spike_index', 
					'dm3':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/drosophila/spikes',
					'dm6':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/drosophila/spikes',
					'hg19':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/human/spikes',
					'hg38':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/human/spikes',
					'mm9':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/mouse/spikes',
					'mm10':'/data/akhtar/group/Giuseppe/INDEXES/FETISH_spikes/mouse/spikes'}

	### iterating over the input files and splitting reads according to their barcodes and read lengths ###
	
	merged_statistics = []
	unmerged_statistics = []

	#Figure out structure of the internal barcode, and start the splitter accordingly
	
	if user_barcode.lower() == 'iclip':
		left_user_barcode, right_user_barcode = 'NNNXXXXNN', None
	elif user_barcode.lower() == 'uvclap':
		left_user_barcode, right_user_barcode = 'NNNXXXXXNN', 'YRRYN'
	elif user_barcode.lower() == 'flash':
		left_user_barcode, right_user_barcode = None, 'NNRRNXXXXXXNN'
	else:
		left_user_barcode = None
		right_user_barcode = None
		user_barcode_list = user_barcode.lower().split(',')
		for i in user_barcode_list:
			if re.search('l=', i):
				left_user_barcode = i.split('l=')[1]
			elif re.search('r=', i):
				right_user_barcode = i.split('r=')[1]
	
	left_params = ([0,0], [0,0], [0,0], [0,0], [0,0])
	right_params = ([0,0], [0,0], [0,0], [0,0], [0,0])
	
	if left_user_barcode and right_user_barcode:
		left_params = GenerateBinary(left_user_barcode)
		right_params = GenerateBinary(right_user_barcode)
	elif left_user_barcode:
		left_params = GenerateBinary(left_user_barcode)
	elif right_user_barcode:
		right_params = GenerateBinary(right_user_barcode)
	else:
		printStatus("No internal barcodes/indices, or your -ub input is incorrect, or perhaps this is not the right tool for you?")
		sys.exit(1)

	print "****************************************"	
	print "Barcoding scheme:", user_barcode
	print "The barcode structure of forward reads: ", left_user_barcode
	print "The barcode structure of reverse reads: ", right_user_barcode
	print "The 'left' parameters", left_params
	print "The 'right' parameters", right_params
	print "****************************************"

	if not ignore: #ignore means "ignore unmerged reads."
		if script_process['BarcodeSplittingUnmerged'] == 'Finished':
			printStatus("Splitting of unmerged R1 and R2 reads was completed in a previous run. Skipping.")
			pass		
		elif script_process['BarcodeSplittingUnmerged'] == 'Started':
			printStatus("Barcode splitting has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob("barcode_sorted/*fastq*")
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Started:Again@ %s\n"%time.strftime("%X %x"))
			
			unmerged_statistics = UnmergedSplitter(forward, reverse, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance)
		
		else:
			script_process['BarcodeSplittingUnmerged'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Started: %s\n"%time.strftime("%X %x"))
			
			unmerged_statistics = UnmergedSplitter(forward, reverse, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance)

		printStatus("Done with the R1 and R2 reads.")	

		
	if bbmerge:
		if script_process['BarcodeSplittingMerged'] == 'Finished':
			printStatus("Merged reads were split successfully in a previous run. Skipping.")
			pass

		elif script_process['BarcodeSplittingMerged'] == 'Started':
			printStatus("Barcode splitting has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob("barcode_sorted/*_merged_*")
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingMerged:Started:Again@ %s\n"%time.strftime("%X %x"))

			merged_statistics = MergedSplitter(gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance)

		else:
			script_process['BarcodeSplittingMerged'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingMerged:Started: %s\n"%time.strftime("%X %x"))
			merged_statistics = MergedSplitter(gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance)
	
	if merged_statistics and unmerged_statistics:
		total_reads = merged_statistics[0] + unmerged_statistics[0]
		total_alien = merged_statistics[1] + unmerged_statistics[1]
		total_bad = merged_statistics[2] + unmerged_statistics[2]
		if total_alien > 0:
			print "\n\n############################################################\n"
			print "Here's some barcode splitting statistics:\n"
			print "Total number of reads:\t %s" %(total_reads)
			print "Total merged reads:\t %s" %merged_statistics[0]
			print "Total unmerged reads:\t %s" %unmerged_statistics[0]
			print "Number of reads that did not conform to the A/B binary barcode scheme:\t %s (%%%s of all reads)" %(total_alien, round((total_alien/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R1 reads" %(merged_statistics[1], unmerged_statistics[1])
			print "Number of reads that contained an unidentifiable barcode:\t %s (%%%s of all reads)" %(total_bad, round((total_bad/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R1 reads\n" %(merged_statistics[2], unmerged_statistics[2])
			for i in range(len(index_barcode_names)):
				print "Reads from replicate A of %s: %s" %(index_barcode_names[i], merged_statistics[3][i] + unmerged_statistics[3][i])
				print "Reads from replicate B of %s: %s\n" %(index_barcode_names[i], merged_statistics[4][i] + unmerged_statistics[4][i])
			print "\n############################################################\n"
		else:
			print "\n\n############################################################\n"
			print "Here's some barcode splitting statistics:\n"
			print "Total number of reads:\t %s" %(total_reads)
			print "Total merged reads:\t %s" %merged_statistics[0]
			print "Total unmerged reads:\t %s" %unmerged_statistics[0]
			print "Number of reads that contained an unidentifiable barcode:\t %s (%%%s of all reads)" %(total_bad, round((total_bad/total_reads)*100))
			print "\t%s of these came from merged reads\n\t%s of these came from R1 reads\n" %(merged_statistics[2], unmerged_statistics[2])
			for i in range(len(index_barcode_names)):
				print "Reads from sample %s: %s" %(index_barcode_names[i], merged_statistics[3][i] + unmerged_statistics[3][i])
			print "\n############################################################\n"			

		cmd_bad_barcodes = "sort bad_barcodes_merged.txt | uniq -c | sort -k1,1n > bad_barcodes_merged_with_counts.txt"
		subprocess.check_call(cmd_bad_barcodes, shell=True)
		cmd_bad_barcodes_R1 = "sort bad_barcodes_R1.txt | uniq -c | sort -k1,1n > bad_barcodes_R1_with_counts.txt"
		subprocess.check_call(cmd_bad_barcodes_R1, shell=True)

	if mapping and aligner == 'bt2-bbmap':
		merged_files = glob.glob("barcode_sorted/*merged*")

		if script_process['bowtie2'] == 'Finished':
			printStatus("Mapping with bowtie2 was successfully completed in a previous run. Skipping.")
			pass
		elif script_process['bowtie2'] == 'Started':
			aa = glob.glob("barcode_sorted/bowtie2/*.bam")
			bb = glob.glob("barcode_sorted/bowtie2/unmapped/*fastq*")
			cc = glob.glob("barcode_sorted/bowtie2/log/*.log")
			files_from_before = aa + bb + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** bowtie2:Started:Again@ %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore)

		else:
			script_process['bowtie2'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** bowtie2:Started: %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore)
		
		unmapped_reads = glob.glob("barcode_sorted/bowtie2/unmapped/*_unmapped.fastq.gz")

		if script_process['bbmap'] == 'Finished':
			printStatus("BBMap has finished running before. Skipping.")
			pass
		elif script_process['bbmap'] == 'Started':
			aa = glob.glob("barcode_sorted/bbmap/*.bam")
			bb = glob.glob("barcode_sorted/bbmap/log/*.log")
			files_from_before = aa + bb
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** bbmap:Started:Again@ %s\n"%time.strftime("%X %x"))

			runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor)

		else:
			script_process['bbmap'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** bbmap:Started: %s\n"%time.strftime("%X %x"))

			runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor)

		bowtie_bams = glob.glob("barcode_sorted/bowtie2/*.bam")
		bbmap_bams = glob.glob("barcode_sorted/bbmap/*.bam")

		if script_process['sambedcleaner'] == 'Finished':
			printStatus("BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass			
		elif script_process['sambedcleaner'] == 'Started':
			aa = glob.glob("barcode_sorted/bbmap/cleaned/*")
			bb = glob.glob("barcode_sorted/bowtie2/cleaned/*")
			cc = glob.glob("*temp*")
			files_from_before = aa + bb + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		else:
			script_process['sambedcleaner'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		# Merge the bam files from bowtie2 and bbmap
		cleaned_files = glob.glob("barcode_sorted/bowtie2/cleaned/*.bam")
		
		if script_process['MergeBowtie2BBMap'] == 'Finished':
			printStatus("BBMap and bowtie2 outputs were merged before. Skipping.")
			pass
		elif script_process['MergeBowtie2BBMap'] == 'Started':
			files_from_before = glob.glob("final_files/bams/*bam*")
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** MergeBowtie2BBMap:Started:Again@ %s\n"%time.strftime("%X %x"))
			BamMerger(cleaned_files, script_process)
		else:
			script_process['MergeBowtie2BBMap'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** MergeBowtie2BBMap:Started: %s\n"%time.strftime("%X %x"))
			BamMerger(cleaned_files, script_process)

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass
		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob("final_files/bam_coverage/bedGraphs/*.bedGraph")
			bb = glob.glob("final_files/bam_coverage/bigwigs/*.bw")
			cc = glob.glob("final_files/xlinks/replicates_merged/bed/*_xlinks.bed")
			dd = glob.glob("final_files/xlinks/bedGraph_coverage/*_coverage.bedGraph")
			ee = glob.glob("final_files/xlinks/bigwig_coverage/*bw")

			files_from_before = aa + bb + cc + dd + ee

			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)
		else:
			script_process['GenerateCov'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)
	
	elif mapping and aligner == 'hisat':
		merged_files = glob.glob("barcode_sorted/*merged*")

		if script_process['hisat'] == 'Finished':
			printStatus("Mapping with HISAT was successfully completed in a previous run. Skipping.")
			pass
		elif script_process['hisat'] == 'Started':
			aa = glob.glob("barcode_sorted/hisat/*.bam")
			cc = glob.glob("barcode_sorted/hisat/log/*.log")
			files_from_before = aa + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** hisat:Started:Again@ %s\n"%time.strftime("%X %x"))
			runHiSat(merged_files, script_process, indices, mapping, processor)

		else:
			script_process['hisat'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** hisat:Started: %s\n"%time.strftime("%X %x"))
			runHiSat(merged_files, script_process, indices, mapping, processor)
	
		hisat_bams = glob.glob("barcode_sorted/hisat/*.bam")
		bbmap_bams = None
		unmapped_reads = glob.glob("barcode_sorted/hisat/unmapped/*_unmapped.fastq.gz")

		if script_process['sambedcleaner'] == 'Finished':
			printStatus("BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass			
		elif script_process['sambedcleaner'] == 'Started':
			aa = glob.glob("barcode_sorted/hisat/cleaned/*")
			cc = glob.glob("*temp*")
			files_from_before = aa + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(hisat_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		else:
			script_process['sambedcleaner'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(hisat_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		# Copy bams to final folder
		cleaned_files = glob.glob("barcode_sorted/hisat/cleaned/*.bam")
		
		if script_process['MoveBams'] == 'Finished':
			printStatus("HISAT .bams moved and indexed. Skipping.")
			pass
		elif script_process['MoveBams'] == 'Started':
			files_from_before = glob.glob("final_files/bams/*bam*")
			for i in files_from_before:
				os.remove(i)
			
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Started:Again@ %s\n"%time.strftime("%X %x"))
			
			for i in cleaned_files:
				cmd = "cp %s final_files/bams/%s" %(i, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)
			
			moved_bams = glob.glob("final_files/bams/*.bam")

			for i in moved_bams:
				cmd = "samtools index %s" %i
				subprocess.check_call(cmd, shell=True)

			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))				
		
		else:
			script_process['MoveBams'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Started: %s\n"%time.strftime("%X %x"))
			
			for i in cleaned_files:
				cmd = "cp %s final_files/bams/%s" %(i, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)
			
			moved_bams = glob.glob("final_files/bams/*.bam")

			for i in moved_bams:
				cmd = "samtools index %s" %i
				subprocess.check_call(cmd, shell=True)
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))		

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass
		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob("final_files/bam_coverage/bedGraphs/*.bedGraph")
			bb = glob.glob("final_files/bam_coverage/bigwigs/*.bw")
			cc = glob.glob("final_files/xlinks/replicates_merged/bed/*_xlinks.bed")
			dd = glob.glob("final_files/xlinks/bedGraph_coverage/*_coverage.bedGraph")
			ee = glob.glob("final_files/xlinks/bigwig_coverage/*bw")

			files_from_before = aa + bb + cc + dd + ee

			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)
		else:
			script_process['GenerateCov'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)
	
	if mapping and aligner == 'bt2':
		merged_files = glob.glob("barcode_sorted/*merged*")

		if script_process['bowtie2-only'] == 'Finished':
			printStatus("Mapping with Bowtie2 was successfully completed in a previous run. Skipping")
			pass
		elif script_process['bowtie2-only'] == 'Started':
			aa = glob.glob("barcode_sorted/bowtie2/*.bam")
			cc = glob.glob("barcode_sorted/bowtie2/log/*.log")
			files_from_before = aa + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** bowtie2-only:Started:Again@ %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore)

		else:
			script_process['bowtie2-only'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** bowtie2-only:Started: %s\n"%time.strftime("%X %x"))
			runBowtie2(merged_files, script_process, indices, mapping, processor, ignore)

		bowtie_bams = glob.glob("barcode_sorted/bowtie2/*.bam")
		bbmap_bams = None
		unmapped_reads = glob.glob("barcode_sorted/bowtie2/unmapped/*_unmapped.fastq.gz")

		if script_process['sambedcleaner'] == 'Finished':
			printStatus("BAM files were de-duplicated and cleaned up with sambedcleaner. Skipping.")
			pass			
		elif script_process['sambedcleaner'] == 'Started':
			aa = glob.glob("barcode_sorted/bowtie2/cleaned/*")
			cc = glob.glob("*temp*")
			files_from_before = aa + cc
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started:Again@ %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(bowtie_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		else:
			script_process['sambedcleaner'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** sambedcleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(hisat_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files)

		# Copy bams to final folder
		cleaned_files = glob.glob("barcode_sorted/bowtie2/cleaned/*.bam")
		
		if script_process['MoveBams'] == 'Finished':
			printStatus("bowtie2-only.bams moved and indexed. Skipping.")
			pass
		elif script_process['MoveBams'] == 'Started':
			files_from_before = glob.glob("final_files/bams/*bam*")
			for i in files_from_before:
				os.remove(i)
			
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Started:Again@ %s\n"%time.strftime("%X %x"))
			
			for i in cleaned_files:
				cmd = "cp %s final_files/bams/%s" %(i, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)
			
			moved_bams = glob.glob("final_files/bams/*.bam")

			for i in moved_bams:
				cmd = "samtools index %s" %i
				subprocess.check_call(cmd, shell=True)

			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))				
		
		else:
			script_process['MoveBams'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Started: %s\n"%time.strftime("%X %x"))
			
			for i in cleaned_files:
				cmd = "cp %s final_files/bams/%s" %(i, i.split('/')[-1])
				subprocess.check_call(cmd, shell=True)
			
			moved_bams = glob.glob("final_files/bams/*.bam")

			for i in moved_bams:
				cmd = "samtools index %s" %i
				subprocess.check_call(cmd, shell=True)
			with open("logfile.log", 'a') as ff:
				ff.write("** MoveBams:Finished: %s\n"%time.strftime("%X %x"))		

		# Generate coverage files for individual replicates
		if script_process['GenerateCov'] == 'Finished':
			printStatus("Coverages were generated before. Skipping")
			pass
		elif script_process['GenerateCov'] == 'Started':
			aa = glob.glob("final_files/bam_coverage/bedGraphs/*.bedGraph")
			bb = glob.glob("final_files/bam_coverage/bigwigs/*.bw")
			cc = glob.glob("final_files/xlinks/replicates_merged/bed/*_xlinks.bed")
			dd = glob.glob("final_files/xlinks/bedGraph_coverage/*_coverage.bedGraph")
			ee = glob.glob("final_files/xlinks/bigwig_coverage/*bw")

			files_from_before = aa + bb + cc + dd + ee

			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started:Again@ %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)
		else:
			script_process['GenerateCov'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** GenerateCov:Started: %s\n"%time.strftime("%X %x"))
			GenerateCoverage(script_process, indices, mapping)

	if spikes_mapper:
		if script_process['spikes'] == 'Finished':
			printStatus("Spikes were already mapped successfully in a previous run. Skipping.")
			pass
		elif script_process['spikes'] == 'Started':
			printStatus("Spikes mapping has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob("barcode_sorted/*")
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** Spikes mapping:Started:Again@ %s\n"%time.strftime("%X %x"))

			spikeMapper(unmapped_reads, script_process, spike_index, mapping, processor)	

		else:
			script_process['spikes'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** Spikes mapping:Started: %s\n"%time.strftime("%X %x"))
		
			spikeMapper(unmapped_reads, script_process, spike_index, mapping, processor)	

		spike_bams = glob.glob("spikes/*.bam")	

		if script_process['spikesBamCleaner'] == 'Finished':
			printStatus("Spikes were de-duplicated and cleaned up with BamCleaner. Skipping.")
			pass
		elif script_process['spikesBamCleaner'] == 'Started':
			printStatus("BamCleaner on Spikes has been initiated before but not completed. Deleting previous files and restarting.")
			files_from_before = glob.glob("spikes/cleaned/*")
			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** Spikes BamCleaner:Started:Again@ %s\n"%time.strftime("%X %x"))
			
			runSamBedCleaner(spike_bams, script_process, chunk, index_barcode_names,  )	
		else:
			script_process['spikesBamCleaner'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** Spikes BamCleaner:Started: %s\n"%time.strftime("%X %x"))

			runSamBedCleaner(spike_bams, script_process, chunk, index_barcode_names,  )	

	# Merge replicates 
	if reps_merged:
		if script_process['MergeReps'] == 'Finished':
			printStatus("Replicates were merged before. Skipping.")
			pass
		elif script_process['MergeReps'] == 'Started':
			aa = glob.glob("final_files/bams/replicates_merged/*bam*")
			bb = glob.glob("final_files/bam_coverage/bedGraphs/replicates_merged/*.bedGraph")
			cc = glob.glob("final_files/bam_coverage/bigwigs/replicates_merged/*.bw")
			dd = glob.glob("final_files/xlinks/replicates_merged/bed/*_xlinks.bed")
			ee = glob.glob("final_files/xlinks/replicates_merged/bedGraph/*.bedGraph")
			gg = glob.glob("final_files/xlinks/replicates_merged/bigwig/*.bw")

			files_from_before = aa + bb + cc + dd + ee + gg

			for i in files_from_before:
				os.remove(i)
			with open("logfile.log", 'a') as ff:
				ff.write("** MergeReps:Started:Again@ %s\n"%time.strftime("%X %x"))
		
			RepMerger(index_barcode_names, script_process, indices, mapping, aligner)
		
		else:
			script_process['MergeReps'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** MergeReps:Started: %s\n"%time.strftime("%X %x"))
		
			RepMerger(index_barcode_names, script_process, indices, mapping, aligner)	

	if gzipafter:
		fastq_files = glob.glob("barcode_sorted/*.fastq")
		if len(fastq_files) > 0:
			cmd = "gzip barcode_sorted/*.fastq"
			printStatus("Gzipping all the fastq files generated after barcode splitting.")
			subprocess.check_call(cmd, shell=True)
			printStatus("Gzipping done.")


def main():
	parser = argparse.ArgumentParser(description='Demultiplexer for FLASH data. Takes in fastq files, list of barcodes, throws out the reads separated for their barcode. Can be used also for iCLIP, uvCLAP or PAR-iCLIP. Just set your internal barcoding scheme with the -ub option.')
	parser.add_argument('-f', '--forward', metavar='path/to/left_reads_R1.fastq.gz', help='The fastq file with the forward reads', required=True)
	parser.add_argument('-r', '--reverse', metavar='path/to/right_reads_R2.fastq.gz', help='The fastq file with the reverse reads')
	parser.add_argument('-b', '--barcodes', metavar='path/to/barcodes.fa', help='The fasta file that contains the barcodes.' ,required=True)
	parser.add_argument('-bm', '--bbmerge', action='store_true', help='Use SeqPrep (you need to have it in your PATH) to preprocess the reads. Highly recommended for uvCLAP data. Using this option will create and use the merged data from SeqPrep (saved as seqprep_merged.fastq.gz).')
	parser.add_argument('-i', '--ignore', action='store_true', help='Only use the merged file generated by BBMerge and ignore the non-mergable ones, assumption being that there\'s something wrong with them. Use only with the -bm option, otherwise the program won\'t do anything.')
	parser.add_argument('-c', '--chunk', type=int, help='Determine the chunk size. The formula is chunk * 100000 per file. Default is 2.', default=2)
	parser.add_argument('-gz', '--gzip', action='store_true', help='Gzip the fastq files on the fly to save space. Slower. I just use -gza instead of this, but this is also an option.')
	parser.add_argument('-mr', '--reps_merged', action='store_true', help='Merge the biological replicates A and B in addition to the normal stuff. Slower')
	parser.add_argument('-m', '--mapping', choices=['hg19','hg38', 'dm3', 'dm6', 'mm10', 'mm9', 'hg19_simple'])
	parser.add_argument('-gza', '--gzipafter', action='store_true', help='Runs gzip on the split files, AFTER everything else is done. Much faster than -gz option.')
	parser.add_argument('-ub', '--user_barcode', type=str, help='User defined internal barcode. For iCLIP it would be L=NNNXXXXNN; for uvCLAP it is L=NNNXXXXXNN,R=YRRYN; for FLASH it is R=NNYYNXXXXXXNN. Default is FLASH', default='FLASH')
	parser.add_argument('-hm', '--hamming_distance', type=int, help="Maximum allowed mutations in the internal index. Default is 1 mutation allowed. For shorter barcodes, can 0 (i.e. perfect match). For longer barcodes can use larger numbers. Experiment with this number for best results", default=1)
	parser.add_argument('-p', '--processor', type=int, help="Number of processors to use. This is only passed to bowtie2 BBMap, and bbduk. Default = 8 ", default=8)
	parser.add_argument('-mi', '--minimum_insert', type=int, help='Minimun insert size for BBMerge. Default is 30 (with all the internal barcodes, so subtract that to get the actualy minimum insert size you are setting.)', default=30)
	parser.add_argument('-al', '--aligner', type=str, choices=['bt2', 'bbmap', 'bt2-bbmap', 'hisat'], default='bt2-bbmap', help='Aligner to use. Default is bowtie2 first, and BBMap on the reads that couldn\'t be aligned with bowtie2')
	parser.add_argument('-tr', '--trim', type=str, default='TACACTCTTTCCCTACACGACGCTCTTCCGATCT,AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help='Remove adapters with BBDuk. Default is Illumina adapters')
	parser.add_argument('-s', '--spikes_mapper', action='store_true', help='Use Bowtie2 to map the spike-ins.')
	args = parser.parse_args()
	current_dir = os.getcwd()+"/"

	script_process = {"BBMerge":None, "BarcodeSplittingUnmerged":None, "BarcodeSplittingMerged":None, 
	"bowtie2":None, "bbmap":None, "sambedcleaner":None, "MergeBowtie2BBMap":None, "MergedSort":None, 
	"MergeReps":None, "GenerateCov":None, "hisat":None, "MoveBams":None, "bowtie2-only":None, "bbmap-only":None,
	"spikes":None, "spikesBamCleaner":None}

	try:
		with open("logfile.log", 'r') as ff:
			printStatus("Reading the logfile.")
			while True:
				try:
					line = ff.next().strip()
					if line.split(" ")[0] == "**":
						process, result = line.split(" ")[1].split(":")[0], line.split(" ")[1].split(":")[1]
						script_process[process] = result
				except StopIteration:
					break
	except IOError:
		printStatus("Fresh start.")
		with open("logfile.log", 'w') as ff:
			date_time = time.strftime("%X %x")
			ff.write("Starting at: %s\n\nExecuting script with the command:\n\n%s\n\n" %(date_time, " ".join(sys.argv[:])))
	
	if re.search('fastq', args.forward, re.IGNORECASE):
		if re.search('fastq.gz', args.forward, re.IGNORECASE):
			filetype = 'fastq.gz'
		else:
			filetype = 'fastq'
	else:
		print("I don't understand your file type. Your files should end either with .fastq or .fastq.gz")
		sys.exit(1)

	name_forward = 'bbmerged_' + os.path.splitext(args.forward)[0].split('/')[-1].split(".fastq")[0] + '.fastq.gz'
	name_reverse = 'bbmerged_' + os.path.splitext(args.reverse)[0].split('/')[-1].split(".fastq")[0] + '.fastq.gz'

	if args.bbmerge:
		if script_process['BBMerge'] == 'Finished':
			printStatus("Using the output of a previous BBmerge run. If you don't want this, stop the script, delete bbmap_merged.fastq.gz and other BBmerge files and start again.")
			processFastqFiles(name_forward, name_reverse, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000, 
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance+1, 
				args.processor, args.aligner, args.spikes_mapper)
		elif script_process['BBMerge'] == 'Started':
			printStatus("BBMerge was initiated but not completed, deleting partially processed files and restarting.")
			os.remove(name_forward)
			os.remove(name_reverse)
			os.remove('bmap_merged.fastq.gz')
			RunBBMerge(script_process, args, current_dir, name_forward, name_reverse, args.minimum_insert)
			print "\n"
			printStatus("Now running the splitter.")	
			processFastqFiles(name_forward, name_reverse, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000, 
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance+1, 
				args.processor, args.aligner, args.spikes_mapper)
		else:
			script_process['BBMerge'] == 'Started'
			with open("logfile.log", 'a') as ff:
				ff.write("** BBMerge:Started: %s\n"%time.strftime("%X %x"))
			RunBBMerge(script_process, args, current_dir, name_forward, name_reverse, args.minimum_insert)
			printStatus("Now running the splitter.")	
			processFastqFiles(name_forward, name_reverse, args.barcodes, filetype, args.bbmerge, args.ignore, args.chunk*100000, 
				args.mapping, args.gzip, args.reps_merged, script_process, args.gzipafter, args.user_barcode, args.hamming_distance+1, 
				args.processor, args.aligner, args.spikes_mapper)
	else:
		processFastqFiles(args.forward, args.reverse, args.barcodes, filetype, False, False, args.chunk*100000, 
			False, args.gzip, False, script_process, args.gzipafter, args.user_barcode, args.hamming_distance+1, 
			args.processor, args.spikes_mapper)

####################################################
# Some useful functions here 
####################################################

def extract_barcodes(barcodes):
	"""This script reads the barcodes fasta file and returns two lists, 
	one with the names of the samples, the other with the barcodes """
	barcode_names = []
	barcode_seqs =[]
	with open(barcodes, 'r') as ff:
		while True:
			try:
				barcode_names.append(ff.next().strip().split('>')[1])
				barcode_seqs.append(ff.next().strip())
			except StopIteration:
				break
		return barcode_names, barcode_seqs

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

#def WriteLog(msg):
#	with open("logfile.log", 'a') as ff:
#		ff.write("** %s %s\n"%(msg,time.strftime("%X %x")))

def RunBBMerge(script_process, args, current_dir, name_forward, name_reverse, minimum_insert):
	""" This just runs BBMerge on the input files. Doesn't return anything (probably returns 1 or 0?)
	The output is named bbmap_merged.fastq.gz. This cannot be changed right now, but I guess should be changed at some point."""

	printStatus("Running BBMerge now.\n")
	BBMERGE_COMMAND = "mininsert=%s in1=%s in2=%s out=%s outu1=%s outu2=%s" %(minimum_insert, os.path.abspath(args.forward), os.path.abspath(args.reverse), 
		current_dir+"bbmap_merged.fastq.gz", current_dir+name_forward, current_dir+name_reverse)
	cmd = "bbmerge.sh %s" %(BBMERGE_COMMAND)
	subprocess.check_call(cmd, shell=True)
	printStatus("Finished BBmerge.")
	script_process['BBMerge'] = 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** BBMerge:Finished: %s\n"%time.strftime("%X %x"))

def RunBBDuk():
	"""Remove adapters on unmergable reads or just unmerged reads"""

def GenerateBinary(user_barcode):
	"""Takes in a string that is composed of n, x, y or r letters
	such as 'nnnxxxyrrynn' (lower case), returns two sets that contain all possible combinations
	of the binary part (y/r) in DNA language, index of the binary barcode as a list, where the first elemenent
	is the beginning position and the second the end, 
	index of the internal index as a list (same as above), and finall a list that contains
	the positions of random nucleotides as its elemenents."""

	barcode = user_barcode.lower()
	code_a = ''.join(''.join(barcode.split('n')).split('x'))
	code_a_index = [ barcode.index(code_a), barcode.index(code_a) + len(code_a)]
	internal_index = ''.join(''.join(''.join(barcode.split('r')).split('y')).split('n'))
	internal_index_index = [barcode.index(internal_index), barcode.index(internal_index) + len(internal_index)]
	all_indices = [x for x in range(len(barcode))]
	
	#remove binary and internal index from the user barcode, what is left is the positions of the random barcode
	#e.g. read = 'ABNHJKIOLKJHDS' 
	#rand = ''.join([read[x] for x in all_indices]) will give the random part as a string

	for i in range(code_a_index[0],code_a_index[1]):
		all_indices.remove(i)
	
	for i in range(internal_index_index[0],internal_index_index[1]):
		all_indices.remove(i)

	convert = {'y':'r', 'r':'y'}
	code_b = "".join([convert[x] for x in code_a])

	letters = {'r':['A', 'G'], 'y':['C', 'T']}
	
	indices = list(itertools.product([0, 1], repeat=len(code_a)))
	clist_a = [x for x in code_a]
	clist_b = [x for x in code_b]
	bcode_a = set()
	bcode_b = set()

	for i in indices:
		new_list = []
		for x in range(len(code_a)):
			a = i[x]
			new_list.append(letters[clist_a[x]][a])
		new_code = "".join(new_list)
		bcode_a.add(new_code)

	for i in indices:
		new_list = []
		for x in range(len(code_b)):
			b = i[x]
			new_list.append(letters[clist_b[x]][b])
		new_code = "".join(new_list)
		bcode_b.add(new_code)
	
	return bcode_a, bcode_b, code_a_index, internal_index_index, all_indices

def UnmergedSplitter(forward, reverse, gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, file_type, left_user_barcode, right_user_barcode, hamming_distance):
	""" This is the first splitter function. This one works on unmerged R1 and R2 files. 
	It basically does 3 things at the same time:
	1) Finds to what sample the reads belongs to. Looks for an index that mathces with at most 1 mismatch
	2) See if the reads has an A or B binary barcode, and split accordingly
	3) Collect splitting statistics to quickly check how many reads comes from which sample in the end """

	printStatus("Will process R1 and R2 files now.")
	if file_type == 'fastq.gz':
		left_reads = gzip.open(forward, 'rb')
		right_reads = gzip.open(reverse, 'rb')
	elif file_type == 'fastq':
		left_reads = open(forward, 'rb')
		right_reads = open(reverse, 'rb')
	multiplier = 1
	count1 = 0 #total number of reads
	alien_count = 0 # not A or B.
	bad_barcode = 0 # A or B but unknown barcode.

	single_reads_f_dict = {}
	single_reads_r_dict = {}
	bad_set = []
	a_barcodes_unmerged = [0 for x in range(len(index_barcode_names))]
	b_barcodes_unmerged = [0 for x in range(len(index_barcode_names))]

	while True:
		try:
			name_f, name_r	= left_reads.next().strip(), right_reads.next().strip()
			seq_f,	seq_r	= left_reads.next().strip(), right_reads.next().strip()
			plus_f, plus_r	= left_reads.next().strip(), right_reads.next().strip()
			qual_f, qual_r	= left_reads.next().strip(), right_reads.next().strip()
			count1 += 1
			# collect the random barcode from both ends, wherever it is defined

			if sum(left_params[4]) != 0 and sum(right_params[4]) != 0:
				ran_bcode = ''.join([seq_f[x] for x in left_params[4]]) + ''.join([seq_r[x] for x in right_params[4]])
				left_trim_l = len(left_user_barcode)
				right_trim_l = -len(right_user_barcode)
				left_trim_r = -len(left_user_barcode)
				right_trim_r = len(right_user_barcode)

			elif sum(left_params[4]):
				ran_bcode = ''.join([seq_f[x] for x in left_params[4]])
				left_trim_l = len(left_user_barcode)
				right_trim_l = None
				left_trim_r = -len(left_user_barcode)
				right_trim_r = None

			elif sum(right_params[4]):
				ran_bcode = ''.join([seq_r[x] for x in right_params[4]])
				left_trim_l = None
				right_trim_l = -len(right_user_barcode)
				left_trim_r = None
				right_trim_r = len(right_user_barcode)

			# also add the randomness from the binary part, if it exists
			if sum(left_params[2]):
				ran_bcode += seq_f[left_params[2][0]:left_params[2][1]]
			elif sum(right_params[2]):
				ran_bcode += seq_r[right_params[2][0]:right_params[2][1]]

			#Get the internal index from one side
			if sum(left_params[3]):
				index_bcode = seq_f[left_params[3][0]:left_params[3][1]]
			elif sum(right_params[3]):
				index_bcode = seq_r[right_params[3][0]:right_params[3][1]]
			
			#Get the binary barcode from one side, if it exists
			if sum(left_params[2]):
				binary_bcode = seq_f[left_params[2][0]:left_params[2][1]]
				list_of_A_barcodes = left_params[0]
				list_of_B_barcodes = left_params[1]
			elif sum(right_params[2]):
				binary_bcode = seq_r[right_params[2][0]:right_params[2][1]]
				list_of_A_barcodes = right_params[0]
				list_of_B_barcodes = right_params[1]
			else:
				binary_bcode = None


			name_f1 = name_f.split(' ')
			name_f1.insert(1, ":" + ran_bcode + " ")
			name_f = ''.join(name_f1)

			name_r1 = name_r.split(' ')
			name_r1.insert(1, ":" + ran_bcode + " ")
			name_r = ''.join(name_r1)

			scores = []

			for i in index_barcodes:
				scores.append(Levenshtein.hamming(i, index_bcode))

			# If binary barcode scheme has been used, split accordingly, assuming that A or B refers to biological replicates or just replicates of some sort
			if binary_bcode:
				if binary_bcode in list_of_A_barcodes: # For YY sequences
					if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
						a_barcodes_unmerged[scores.index(min(scores))] += 1
						try:
							single_reads_f_dict["barcode_sorted/A_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
							single_reads_r_dict["barcode_sorted/A_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
						except KeyError:
							single_reads_f_dict["barcode_sorted/A_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
							single_reads_r_dict["barcode_sorted/A_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
					else:
						bad_barcode += 1
						bad_set.append(index_bcode)

				elif binary_bcode in list_of_B_barcodes: # For RR sequences
					if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
						b_barcodes_unmerged[scores.index(min(scores))] +=1 
						try:
							single_reads_f_dict["barcode_sorted/B_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
							single_reads_r_dict["barcode_sorted/B_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
						except KeyError:
							single_reads_f_dict["barcode_sorted/B_" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
							single_reads_r_dict["barcode_sorted/B_" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
					else:
						bad_set.append(index_bcode)
						bad_barcode += 1

				else:
					alien_count += 1
			else:
				if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 and unique
					a_barcodes_unmerged[scores.index(min(scores))] += 1
					try:
						single_reads_f_dict["barcode_sorted/" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"].append(["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])])
						single_reads_r_dict["barcode_sorted/" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"].append(["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])])
					except KeyError:
						single_reads_f_dict["barcode_sorted/" + index_barcode_names[scores.index(min(scores))] + "_R1.fastq"] = [["%s\n%s\n+\n%s\n" % (name_f, seq_f[left_trim_l:right_trim_l], qual_f[left_trim_l:right_trim_l])]]
						single_reads_r_dict["barcode_sorted/" + index_barcode_names[scores.index(min(scores))] + "_R2.fastq"] = [["%s\n%s\n+\n%s\n" % (name_r, seq_r[right_trim_r:left_trim_r], qual_r[right_trim_r:left_trim_r])]]
				else:
					bad_barcode += 1
					bad_set.append(index_bcode)


			# This is the routine that makes sure that the dictionary is emptied every "chunk" times
			# I guess one could write a function to find out the optimal "chunk" size by checking system recources
			# and figuring out how much one should use. I cannot do this at the moment unfortunately.

			if count1 == multiplier*chunk:
				if gzipp:
					for key, value in single_reads_f_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				
					for key, value in single_reads_r_dict.iteritems():
						with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				else:
					for key, value in single_reads_f_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				
					for key, value in single_reads_r_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				multiplier += 1
				single_reads_f_dict = {}
				single_reads_r_dict = {}

		except StopIteration:
			if gzipp:
				for key, value in single_reads_f_dict.iteritems():
					with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
		
				for key, value in single_reads_r_dict.iteritems():
					with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
			else:	
				for key, value in single_reads_f_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
		
				for key, value in single_reads_r_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
			for i in bad_set:
				with open('bad_barcodes_R1.txt', 'a') as ff:
					ff.write(i+'\n')

			single_reads_f_dict = {}
			single_reads_r_dict = {}

			script_process['BarcodeSplittingUnmerged'] = 'Finished'
			
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingUnmerged:Finished: %s\n"%time.strftime("%X %x"))
			break

	left_reads.close()
	right_reads.close()
	return count1, alien_count, bad_barcode, a_barcodes_unmerged, b_barcodes_unmerged

def MergedSplitter(gzipp, script_process,index_barcode_names, index_barcodes, chunk, left_params, right_params, left_user_barcode, right_user_barcode, hamming_distance):
	""" Same as above, just for the merged reads. 
	The barcodes are reverse complemented to make sure that A for unmerged means A for the merged reads. """
	printStatus("Will process merged reads now.")

	merged_reads_dict = {}
	alien_count = 0
	bad_barcode = 0
	a_barcodes_merged = [0 for x in range(len(index_barcode_names))]
	b_barcodes_merged = [0 for x in range(len(index_barcode_names))]
	count2 = 0
	multiplier = 1
	bad_set = []

	seqprep_merged = gzip.open("bbmap_merged.fastq.gz", 'rb')
	
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N' : 'N'}

	while True:
		try:
			merged_name = seqprep_merged.next().strip()
			merged_seq = seqprep_merged.next().strip()
			merged_plus = seqprep_merged.next().strip()
			merged_qual = seqprep_merged.next().strip()

			count2 += 1
			
			# collect the random barcode from both ends, wherever it is defined
			left_trim = None
			right_trim = None
			if sum(left_params[4]) != 0 and sum(right_params[4]) != 0: 
				ran_bcode = ''.join([merged_seq[x] for x in left_params[4]]) + ''.join(["".join(complement.get(base, base) for base in reversed(merged_seq))[x] for x in right_params[4]])
				left_trim = len(left_user_barcode)
				right_trim = -len(right_user_barcode)
			elif sum(left_params[4]):
				ran_bcode = ''.join([merged_seq[x] for x in left_params[4]])
				left_trim = len(left_user_barcode)
				right_trim = None
			elif sum(right_params[4]):
				ran_bcode = ''.join(["".join(complement.get(base, base) for base in reversed(merged_seq))[x] for x in right_params[4]])
				left_trim = None
				right_trim = -len(right_user_barcode)

			# also add the randomness from the binary part, if it exists
			if sum(left_params[2]):
				ran_bcode += merged_seq[left_params[2][0]:left_params[2][1]]
			elif sum(right_params[2]):
				ran_bcode += merged_seq[right_params[2][0]:right_params[2][1]]

			#Get the internal index from one side
			if sum(left_params[3]):
				index_bcode = merged_seq[left_params[3][0]:left_params[3][1]]
			elif sum(right_params[3]):
				index_bcode = "".join(complement.get(base, base) for base in reversed(merged_seq))[right_params[3][0]:right_params[3][1]]
			
			#Get the binary barcode from one side, if it exists
			if sum(left_params[2]):
				binary_bcode = merged_seq[left_params[2][0]:left_params[2][1]]
				list_of_A_barcodes = left_params[0]
				list_of_B_barcodes = left_params[1]
			elif sum(right_params[2]):
				binary_bcode = "".join(complement.get(base, base) for base in reversed(merged_seq))[right_params[2][0]:right_params[2][1]]
				list_of_A_barcodes = right_params[0]
				list_of_B_barcodes = right_params[1]
			else:
				binary_bcode = None

			merged_name_1 = merged_name.split(' ')
			merged_name_1.insert(1, ":" + ran_bcode + " ")
			merged_name = ''.join(merged_name_1)
			
			scores = []

			for i in index_barcodes:
				scores.append(Levenshtein.hamming(i, index_bcode))
			if binary_bcode:
				if binary_bcode in list_of_A_barcodes: # For RYYR sequences
					if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
						a_barcodes_merged[scores.index(min(scores))] += 1
						try:
							merged_reads_dict["barcode_sorted/A_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
						except KeyError:
							merged_reads_dict["barcode_sorted/A_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
					else:
						bad_barcode += 1
						bad_set.append(index_bcode)

				elif binary_bcode in list_of_B_barcodes: # For YRRY sequences
					if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
						b_barcodes_merged[scores.index(min(scores))] += 1
						try:
							merged_reads_dict["barcode_sorted/B_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
						except KeyError:
							merged_reads_dict["barcode_sorted/B_merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
					else:
						bad_barcode += 1
						bad_set.append(index_bcode)
				
				else:
					alien_count += 1
			else:
				if min(scores) < hamming_distance and scores.count(min(scores)) == 1: #The distance has to be less than 2 (0 or 1) and unique
					a_barcodes_merged[scores.index(min(scores))] += 1
					try:
						merged_reads_dict["barcode_sorted/merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"].append(['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])])
					except KeyError:
						merged_reads_dict["barcode_sorted/merged_" + index_barcode_names[scores.index(min(scores))] + ".fastq"] = [['%s\n%s\n+\n%s\n' %(merged_name, merged_seq[left_trim:right_trim], merged_qual[left_trim:right_trim])]]
				else:
					bad_barcode += 1
					bad_set.append(index_bcode)			

			if count2 == multiplier*chunk:
				if gzipp:
					for key, value in merged_reads_dict.iteritems():
							with gzip.open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
								for i in range(len(value)):
									ff.write('%s' %value[i][0])
				else:
					for key, value in merged_reads_dict.iteritems():
						with open('%s' %key, 'a+') as ff:
							for i in range(len(value)):
								ff.write('%s' %value[i][0])
				multiplier += 1
				merged_reads_dict = {}

		except StopIteration:
			if gzipp:
				for key, value in merged_reads_dict.iteritems():
					with open('%s' %key.split(".fastq")[0]+".fastq.gz", 'ab') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])
			else:
				for key, value in merged_reads_dict.iteritems():
					with open('%s' %key, 'a+') as ff:
						for i in range(len(value)):
							ff.write('%s' %value[i][0])

			for i in bad_set:
				with open('bad_barcodes_merged.txt', 'a') as ff:
					ff.write(i+'\n')

			merged_reads_dict = {}
			script_process['BarcodeSplittingMerged'] == 'Finished'
			with open("logfile.log", 'a') as ff:
				ff.write("** BarcodeSplittingMerged:Finished: %s\n"%time.strftime("%X %x"))
			break
		
	seqprep_merged.close()
	return count2, alien_count, bad_barcode, a_barcodes_merged, b_barcodes_merged

def runBowtie2(merged_files, script_process, indices, mapping, processor, ignore):
	counting_files = 1
	total_files = len(merged_files)

	for i in merged_files:
		if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				if ignore:
					both_files = i 
				else:
					both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				if ignore:
					both_files = i
				else:
					both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]
		else:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				if ignore:
					both_files = i
				else:
					both_files = i + "," + "barcode_sorted/" + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				if ignore:
					both_files = i
				else:
					both_files = i + "," + "barcode_sorted/" + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("merged_")[1].split(".fastq")[0]

		cmd = "bowtie2 -p %s -x %s --un-gz %s -U %s 2> %s | samtools view -b - > %s" %(processor, indices[mapping][0], 
			"barcode_sorted/bowtie2/unmapped/"+ new_name.split("/")[-1] +"_unmapped.fastq.gz", 
			both_files, "barcode_sorted/bowtie2/log/"+ new_name.split("/")[-1] + ".bowtie2.log", "barcode_sorted/bowtie2/"+ new_name.split("/")[-1] +".bam")
	
		printStatus("Running bowtie2 for file (%s of %s): %s " %(counting_files, total_files, i.split("/")[-1]))
		print "\t\t", cmd, "\n"
		
		counting_files += 1
		
		subprocess.check_call(cmd, shell=True)
	
	script_process['bowtie2'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** bowtie2:Finished: %s\n"%time.strftime("%X %x"))

def runBBMap(unmapped_reads, script_process, indices, mapping, merged_files, processor):
	counting_files2 = 1
	total_files = len(merged_files)
	for i in unmapped_reads:
		cmd = "cd %s && bbmap.sh %s threads=%s build=%s in=%s out=%s 2> %s" %(BBMAP_PATH,
																	BBMAP_OPTIONS, processor,
																	indices[mapping][1],
																	current_dir+i, 
																	current_dir+"barcode_sorted/bbmap/"+i.split("/")[-1].split("_unmapped.fastq")[0]+".bam", 
																	current_dir+"barcode_sorted/bbmap/log/"+i.split("/")[-1].split("_unmapped.fastq")[0]+".bbmap.log")
		
		printStatus("Running BBMap on the file (%s of %s): %s" %(counting_files2, total_files, i))
		print "\t\t", cmd, "\n"
		counting_files2 += 1
		subprocess.check_call(cmd, shell=True)

	script_process['bbmap'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** bbmap:Finished: %s\n"%time.strftime("%X %x"))

def runHiSat(merged_files, script_process, indices, mapping, processor):
	counting_files = 1
	total_files = len(merged_files)

	for i in merged_files:
		if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				both_files = i + "," + i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("_merged")[0]+i.split("_merged")[1].split(".fastq")[0]
		else:
			if merged_files[0].split(".fastq")[1] == ".gz": #If the files are gzipped to save space.
				both_files = i + "," + "barcode_sorted/" + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq.gz"
			else:
				both_files = i + "," + "barcode_sorted/" + i.split("merged_")[1].split(".fastq")[0]+"_R1.fastq"

			new_name = i.split("merged_")[1].split(".fastq")[0]

		# Options to add/change:
		# --pen-cansplice <int>              penalty for a canonical splice site (0)
		# --pen-noncansplice <int>           penalty for a non-canonical splice site (12)
		# --min-intronlen <int>              minimum intron length (20)
		# --max-intronlen <int>              maximum intron length (500000)
		# --rna-strandness <string>          Specify strand-specific information (unstranded)
		# --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
		# --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
		# --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
		# --rdg <int>,<int>  read gap open, extend penalties (5,3)
		
		
		hisat_options = "--rna-strandness F --max-intronlen 100000"

				
		cmd = "hisat2 %s -p %s -x %s -U %s --un-gz %s 2> %s | samtools view -b - > %s" %(hisat_options, processor, indices[mapping][3], 
				both_files, "barcode_sorted/hisat/unmapped/"+ new_name("/")[-1] +"_unmapped.fastq.gz",
				"barcode_sorted/hisat/log/"+ new_name.split("/")[-1] + ".hisat.log", "barcode_sorted/hisat/"+ new_name.split("/")[-1] +".bam")
		
		printStatus("Running HISAT for file (%s of %s): %s " %(counting_files, total_files, i.split("/")[-1]))
		print "\t\t", cmd, "\n"
		
		counting_files += 1
		
		subprocess.check_call(cmd, shell=True)
	
	script_process['hisat'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** hisat:Finished: %s\n"%time.strftime("%X %x"))	

def spikeMapper(unmapped_reads, script_process, spike_index, mapping, processor):
	counting_files = 1
	total_files = len(merged_files)
	for i in unmapped_reads:
		cmd = "bowtie2 -p %s -x %s -U %s 2> %s | samtools view -b > %s" %(processor, spike_index[mapping][0],
			i, "spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes.log", 
			"/spikes/"+i.split("/")[-1].split("_unmapped.fastq")[0]+"_spikes.bam")

	printStatus("Running bowtie2 on Spikes for file (%s of %s): %s " %(counting_files, total_files, i.split("/")[-1]))
	print "\t\t", cmd, "\n"

	counting_files += 1

	subprocess.check_call(cmd, shell=True)

	script_process['spikes'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** Spikes mapping:Finished: %s\n"%time.strftime("%X %x"))


def runSamBedCleaner(bowtie_bams, bbmap_bams, script_process, chunk, index_barcode_names, merged_files):
	#### sambedcleaner options:
	keep_bed = False 
	min_quality = 10 # Reads with quality less than this will be discarded. 
	multi = True # This is marker that tells the program that all the .bam files will be processed at the same time.
	
	if aligner == 'bt2-bbmap' and bowtie_bams and bbmap_bams:
		jobs1 = []
		jobs2 = []
		printStatus("Running sambedcleaner on bowtie2 generated .bam files now.")
		for i in bowtie_bams:
			p1 = mp.Process(target=sbp.DupRemover, args=(i, "barcode_sorted/bowtie2/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi))
			jobs1.append(p1)
			p1.start()
		
		printStatus("Running sambedcleaner on BBMap generated .bam files now.")
		for i in bbmap_bams:
			p2 = mp.Process(target=sbp.DupRemover, args=(i, "barcode_sorted/bbmap/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi))
			jobs2.append(p2)
			p2.start()
		
		p1.join()
		p2.join()

		# Normally, the join() thingy should work and the script should proceed only after all the files are processed.
		# For some reason this doesn't work, and the script prematurely proceed only to fail later on.
		# So this routine checks if the number of .bam files is correct, otherwise it will not proceed.
		# This might actually lead to a zombie program that is stuct at this stage if something goes wrong with the 
		# de-duplication without an exit signal. But I can live with this.

		#while True:
		#	aa = glob.glob("barcode_sorted/bowtie2/cleaned/*.bam")
		#	bb = glob.glob("barcode_sorted/bbmap/cleaned/*.bam")
		#	cleaned_files = aa + bb
		#	if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
		#		if len(cleaned_files) == len(index_barcode_names)*4:
		#			break
		#		else:
		#			time.sleep(10)
		#	else:
		#		if len(cleaned_files) == len(index_barcode_names)*2:
		#			break
		#		else:
		#			time.sleep(10)

		time.sleep(120)

	elif aligner == 'bt2' and bowtie_bams:
		jobs1 = []
		printStatus("Running sambedcleaner on bowtie2 generated .bam files now.")
		for i in bowtie_bams:
			p1 = mp.Process(target=sbp.DupRemover, args=(i, "barcode_sorted/hisat/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi))
			jobs1.append(p1)
			p1.start()

		p1.join()
		
		#while True:
		#	cleaned_files = glob.glob("barcode_sorted/hisat/cleaned/*.bam")
		#	if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
		#		if len(cleaned_files) == len(index_barcode_names)*2:
		#			break
		#		else:
		#			time.sleep(10)
		#	else:
		#		if len(cleaned_files) == len(index_barcode_names):
		#			break
		#		else:
		#			time.sleep(10)

		time.sleep(120)

	elif aligner == 'bbmap' and bbmap_bams:
		jobs1 = []
		printStatus("Running bamcleaner on bbmap generated .bam files now.")
		for i in bbmap_bams:
			p1 = mp.Process(target=sbp.DupRemover, args=(i, "barcode_sorted/bbmap/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi))
			jobs1.append(p1)
			p1.start()

		p1.join()
		
		#while True:
		#	cleaned_files = glob.glob("barcode_sorted/bbmap/cleaned/*.bam")
		#	if merged_files[0].split('/')[-1].split('merged')[0] in ['A_', 'B_']:
		#		if len(cleaned_files) == len(index_barcode_names)*2:
		#			break
		#		else:
		#			time.sleep(10)
		#	else:
		#		if len(cleaned_files) == len(index_barcode_names)*1:
		#			break
		#		else:
		#			time.sleep(10)

		time.sleep(120)

	elif aligner == 'hisat' and hisat_bams:
		jobs1 = []
		printStatus("Running bamcleaner on hisat2 gnerated .bam files now")
		p1 = mp.Process(target=sbp.DupRemover, args=(i, "barcode_sorted/bbmap/cleaned/"+i.split("/")[-1].split(".bam")[0]+"_cleaned.bam", keep_bed, chunk, min_quality, multi))
			jobs1.append(p1)
			p1.start()

		p1.join

		time.sleep(120)

	script_process['sambedcleaner'] == 'Finished'

	with open("logfile.log", 'a') as ff:
		ff.write("** sambedcleaner:Finished: %s\n"%time.strftime("%X %x"))

def BamMerger(cleaned_files, script_process):
	""" This just merges the bowtie2 and bbmap generated .bams to get the best of both worlds. """

	printStatus("Merging and sorting and indexing bowtie2 and BBMap outputs.")

	for i in cleaned_files:
		cmd = "samtools merge -h %s %s %s %s" %(i, "final_files/bams/temp-"+i.split("/")[-1].split("_cleaned")[0]+"_cleaned+merged.bam", i, "barcode_sorted/bbmap/cleaned/"+i.split("/")[-1])
		subprocess.check_call(cmd, shell=True)
	
	merged_bams = glob.glob("final_files/bams/*.bam")
	
	for i in merged_bams:
		cmd = "samtools sort %s %s" %(i, "final_files/bams/"+i.split("temp-")[-1].split(".bam")[0])
		printStatus("Sorting .bam file: %s"%i.split("/")[-1])
		subprocess.check_call(cmd, shell=True)

	for i in merged_bams:
		os.remove(i)
	
	merged_bams2 = glob.glob("final_files/bams/*.bam")

	for i in merged_bams2:
		cmd = "samtools index %s" %i
		subprocess.check_call(cmd, shell=True)
	
	script_process['MergeBowtie2BBMap'] == 'Finished'

	with open("logfile.log", 'a') as ff:
		ff.write("** MergeBowtie2BBMap:Finished: %s\n"%time.strftime("%X %x"))

def RepMerger(index_barcode_names, script_process, indices, mapping, aligner):
	""" With too many profiles at hand, it doesn't make sense to look at replicates all the time. 
	This function merges A and B replicates into one and calculates a bunch of coverage files using 
	genomeCoverageBed and converts them into bigwigs using bedGraphToBigWig."""
	
	printStatus("Merging biological replicates.")
	if aligner == 'bt2-bbmap':
		for i in index_barcode_names:
			name1 = "final_files/bams/A_%s_cleaned+merged.bam" %i
			name2 = "final_files/bams/B_%s_cleaned+merged.bam" %i
			out = "final_files/bams/replicates_merged/A+B_%s_cleaned+merged.bam" %i
			cmd = "samtools merge -h %s %s %s %s" %(name1, out, name1, name2)
			printStatus("Merging to generate:\n\t%s"%out)
			subprocess.check_call(cmd, shell=True)
	else:
		for i in index_barcode_names:
			name1 = "final_files/bams/A_%s_cleaned.bam" %i
			name2 = "final_files/bams/B_%s_cleaned.bam" %i
			out = "final_files/bams/replicates_merged/A+B_%s_cleaned.bam" %i
			cmd = "samtools merge -h %s %s %s %s" %(name1, out, name1, name2)
			printStatus("Merging to generate:\n\t%s"%out)
			subprocess.check_call(cmd, shell=True)		

	merged_reps_bams = glob.glob("final_files/bams/replicates_merged/*.bam")
	
	printStatus("Generating indices for those merged bams.")

	for i in merged_reps_bams:
		cmd = "samtools index %s" %i
		subprocess.check_call(cmd, shell=True)

	printStatus("Generating coverage files for those merged bams.")
	for i in merged_reps_bams:
		cmd1 = "genomeCoverageBed -bg -strand - -ibam %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/bam_coverage/bedGraphs/replicates_merged/"+i.split("/")[-1].split(".bam")[0]+"_minus_coverage.bedGraph")
		cmd2 = "genomeCoverageBed -bg -strand + -ibam %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/bam_coverage/bedGraphs/replicates_merged/"+i.split("/")[-1].split(".bam")[0]+"_plus_coverage.bedGraph")
		printStatus("Calculating coverage for the minus strand for the file: %s"%i)
		subprocess.check_call(cmd1, shell=True)
		printStatus("Calculating coverage for the plus strand for the file: %s"%i)
		subprocess.check_call(cmd2, shell=True)
	bedgraph_files_merged = glob.glob("final_files/bam_coverage/bedGraphs/replicates_merged/*.bedGraph")
	printStatus("Generating BigWigs.")
	
	for i in bedgraph_files_merged:
		if os.path.getsize(i) > 0:
			cmd = "bedGraphToBigWig %s %s %s" %(i, indices[mapping][2], "final_files/bam_coverage/bigwigs/replicates_merged/"+i.split("/")[-1].split(".bedGraph")[0]+'.bw')
			subprocess.check_call(cmd, shell=True)
		else:
			printStatus("Skipping file %s because it's empty." %i.split("/")[-1])

	printStatus("Generating the 'cross-linked' nucleotide file from the replicates merged bam files")

	for i in merged_reps_bams:
		cmd = "bedtools bamtobed -i %s | awk '{ if ($6 == \"+\") print $1,$2-1,$2,$4,$5,$6; else print $1,$3,$3+1,$4,$5,$6}' | awk '{if($2>0) print $0}' | tr [:blank:] \\\\t > %s" %(i, "final_files/xlinks/replicates_merged/bed/"+i.split("/")[-1].split(".bam")[0]+"_xlinks.bed")
		printStatus("Generating cross-linked nucleotide for the file: %s"%i)
		subprocess.check_call(cmd, shell=True)
	printStatus("Calculating the coverages of the replicates-merged x-linked nucleotide .bed files")
	
	xlink_bed_merged = glob.glob("final_files/xlinks/replicates_merged/bed/*_xlinks.bed")
	
	for i in xlink_bed_merged:
		cmd1 = "genomeCoverageBed -bg -strand - -i %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/xlinks/replicates_merged/bedGraph/"+i.split("/")[-1].split(".bed")[0]+"_minus_coverage.bedGraph")
		cmd2 = "genomeCoverageBed -bg -strand + -i %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/xlinks/replicates_merged/bedGraph/"+i.split("/")[-1].split(".bed")[0]+"_plus_coverage.bedGraph")
		printStatus("Calculating minus strand coverage of the file: %s"%i)
		subprocess.check_call(cmd1, shell=True)
		printStatus("Calculating plus strand coverage of the file: %s"%i)
		subprocess.check_call(cmd2, shell=True)
	
	printStatus("Generating bigwigs")
	coverage_bed_merged = glob.glob("final_files/xlinks/replicates_merged/bedGraph/*.bedGraph")
	for i in coverage_bed_merged:
		if os.path.getsize(i) > 0:
			cmd = "bedGraphToBigWig %s %s %s" %(i, indices[mapping][2], "final_files/xlinks/replicates_merged/bigwig/"+i.split("/")[-1].split(".bedGraph")[0]+".bw")
			subprocess.check_call(cmd, shell=True)
	script_process['MergeReps'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** MergeReps:Finished: %s\n"%time.strftime("%X %x"))

def GenerateCoverage(script_process, indices, mapping):
	
	merged_bams2 = glob.glob("final_files/bams/*.bam")
	
	for i in merged_bams2:
		cmd1 = "genomeCoverageBed -bg -strand - -ibam %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/bam_coverage/bedGraphs/"+i.split("/")[-1].split(".bam")[0]+"_minus_coverage.bedGraph")
		cmd2 = "genomeCoverageBed -bg -strand + -ibam %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/bam_coverage/bedGraphs/"+i.split("/")[-1].split(".bam")[0]+"_plus_coverage.bedGraph")
		printStatus("Calculating coverage for the minus strand for the file: %s"%i)
		subprocess.check_call(cmd1, shell=True)
		printStatus("Calculating coverage for the plus strand for the file: %s"%i)
		subprocess.check_call(cmd2, shell=True)

	bedgraph_files = glob.glob("final_files/bam_coverage/bedGraphs/*.bedGraph")
	
	printStatus("Generating BigWigs.")
	for i in bedgraph_files:
		if os.path.getsize(i) > 0:
			cmd = "bedGraphToBigWig %s %s %s" %(i, indices[mapping][2], "final_files/bam_coverage/bigwigs/"+i.split("/")[-1].split(".bedGraph")[0]+'.bw')
			subprocess.check_call(cmd, shell=True)
		else:
			printStatus("Skipping file %s because it's empty." %i.split("/")[-1])

	# Prepare the cross-linked nucleotide files:

	for i in merged_bams2:
		printStatus("Generating the 'cross-linked' nucleotide file from the bam file: %s" %i.split("/")[-1])
		cmd = "bedtools bamtobed -i %s | awk '{ if ($6 == \"+\") print $1,$2-1,$2,$4,$5,$6; else print $1,$3,$3+1,$4,$5,$6}' | awk '{if($2>0) print $0}' | tr [:blank:] \\\\t > %s" %(i, "final_files/xlinks/bed/"+i.split("/")[-1].split(".bam")[0]+"_xlinks.bed")
		subprocess.check_call(cmd, shell=True)

	xlink_bed = glob.glob("final_files/xlinks/bed/*_xlinks.bed")

	printStatus("Generating coverages for those .bed files")
	for i in xlink_bed:
		cmd1 = "genomeCoverageBed -bg -strand - -i %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/xlinks/bedGraph_coverage/"+i.split("/")[-1].split(".bed")[0]+"_minus_coverage.bedGraph")
		cmd2 = "genomeCoverageBed -bg -strand + -i %s -g %s | sort -k1,1 -k2,2n > %s" %(i, indices[mapping][2], "final_files/xlinks/bedGraph_coverage/"+i.split("/")[-1].split(".bed")[0]+"_plus_coverage.bedGraph")
		subprocess.check_call(cmd1, shell=True)
		subprocess.check_call(cmd2, shell=True)

	coverage_bed = glob.glob("final_files/xlinks/bedGraph_coverage/*_coverage.bedGraph")

	printStatus("Making bigwigs out of bedGraphs")
	for i in coverage_bed:
		if os.path.getsize(i) > 0:
			cmd = "bedGraphToBigWig %s %s %s" %(i, indices[mapping][2], "final_files/xlinks/bigwig_coverage/"+i.split("/")[-1].split(".bedGraph")[0]+".bw")
			subprocess.check_call(cmd, shell=True)
	script_process['GenerateCov'] == 'Finished'
	with open("logfile.log", 'a') as ff:
		ff.write("** GenerateCov:Finished: %s\n"%time.strftime("%X %x"))

if __name__ == "__main__":
	start_time = time.time()
	printStatus("Initializing the script.")
	main()
	end_time = time.time() - start_time
	printStatus("Done. Elapsed time: %s seconds." %round(end_time))