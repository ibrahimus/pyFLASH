#!/usr/bin/python

from __future__ import division
import hashlib
import os
import glob
import argparse
import time
import gzip
import sys




def FileBreaker(fastq_file, chunk, output, file_type):
	full_path = os.path.abspath(fastq_file)
	full_output_path = os.path.abspath(output)
	path = '/'.join(full_path.split('/')[:-1]) + '/'
	file_name = full_path.split('/')[-1]
	out_path = '/'.join(full_output_path.split('/')[:-1]) + '/'
	count1 = 0

#iteration 1: 
	printStatus("Creating a dictionary with the sha1 sum of read + barcode pairs with the name of the read as the 'value'")

	if file_type =='fastq.gz':
		fastq_file_handle = gzip.open(fastq_file, 'rb')
	elif file_type == 'fastq':
		fastq_file_handle = open(fastq_file, 'r')
	else:
		printStatus("Don't know this file type, exiting.")
		sys.exit(1)

	sha1_dict = {}
	quality_score = [x for x in range(42)]
	qualities_solexa = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
	qual_list = [x for x in qualities_solexa]

	while True:
		try:
			name = fastq_file_handle.next().strip()
			seq = fastq_file_handle.next().strip()
			fastq_file_handle.next()
			qual_ascii = fastq_file_handle.next().strip()
			quality = 0
			count1 += 1

			for i in qual_ascii:
				quality += quality_score[qualities_solexa.index(i)]


			mm = hashlib.sha1()
			try:
				mm.update(seq+name.split(' ')[0].split(':')[7])
			except IndexError:
				type_of_opeation = True
				mm.update(seq)
			
			digest = mm.hexdigest()

			try:
				if sha1_dict[digest][1] < quality:
					sha1_dict[digest] = [name, quality]
			except KeyError:
				sha1_dict[digest] = [name, quality]

		except StopIteration:
			break

	fastq_file_handle.close()

	if type_of_opeation:
		printStatus("WARNING!\n\tThis is probably not uvCLAP/FLASH data, falling back to normal duplicate removal. Check your fastq headers.")

	printStatus("Creating a dictionary with unique read + barcode pairs and writing them into the final file.")
#iteration 2: Writing the unique reads into a new file.
	count2 = 0
	count3 = 0
	multiplier = 1
	if file_type =='fastq.gz':
		fastq_file_handle = gzip.open(fastq_file, 'rb')
	elif file_type == 'fastq':
		fastq_file_handle = open(fastq_file, 'r')
	else:
		printStatus("Don't know this file type, exiting.")
		sys.exit(1)

	seqs_to_write = {}

	while True:
		try:
			name = fastq_file_handle.next().strip()
			seq = fastq_file_handle.next().strip()
			fastq_file_handle.next()
			quality = fastq_file_handle.next().strip()

			mm = hashlib.sha1()
			try:
				mm.update(seq+name.split(' ')[0].split(':')[7])
			except IndexError:
				mm.update(seq)
			digest = mm.hexdigest()

			if sha1_dict[digest][0] == name:
				seqs_to_write[name] = [seq, quality]
				count2 += 1
			else:
				count3 += 1

			if count2 == chunk*multiplier:
				for key, value in seqs_to_write.iteritems():
					with open(output, 'a+') as ff:
						ff.write('%s\n%s\n+\n%s\n' %(key,value[0],value[1]))
				multiplier += 1
				seqs_to_write = {}

		except StopIteration:
			for key, value in seqs_to_write.iteritems():
				with open(output, 'a+') as ff:
					ff.write('%s\n%s\n+\n%s\n' %(key,value[0],value[1]))
			seqs_to_write = {}
			break

	fastq_file_handle.close()
	
	removed_reads_percentage = round((count3/count1)*100)
	printStatus("%s reads were removed. That is %%%s of the reads." %(count3,removed_reads_percentage))

def main():
	parser = argparse.ArgumentParser(description='Cleans up duplicates in fastq files. Requires the 10nt barcode to be written to the head with pyCLAP. Use with extreme caution.')
	parser.add_argument('-f', '--file', metavar='path/to/TheFileToBeProcessed.fastq.gz', help='The gzip compressed fastq file. ', required=True)
	parser.add_argument('-o', '--output', metavar='ProcessedFile.fastq', required=True)
	parser.add_argument('-c', '--chunk', type=int, help='Determine the chunk size. The formula is chunk * 100000 per file. Default is 2.', default=1)
	parser.add_argument('-t', '--type', type=str, choices=['fastq', 'fastq.gz'])
	args = parser.parse_args()

	FileBreaker(args.file, args.chunk*100000, args.output, args.type)

def printStatus(msg):
	print "%s: %s" % (time.strftime("%X %x"), msg)

if __name__ == "__main__":
	start_time = time.time()
	printStatus("Initializing the fastq_cleaner script.")
	main()
	end_time = round(time.time() - start_time)
	printStatus("Done. Elapsed time: %s seconds." %end_time)