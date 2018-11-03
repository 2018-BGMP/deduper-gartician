#!/usr/env/bin python3 

# Title: De-Duper 

import argparse 
import re 
from collections import defaultdict 

paired = None

def get_arguments():
	parser = argparse.ArgumentParser(description="Removes single-end PCR duplicates from a regularly formatted SAM file.")
	parser.add_argument("-p", "--paired", required=False, help="OPTIONAL: Denotes file contains paired end alignment reads (SAM format).", action="store_true")
	parser.add_argument("-u", "--umi", required=True, type=str, help="REQUIRED: Select a file containing UMI's (delimited by new line).")
	parser.add_argument("-f", "--file", required=True, type=str, help="REQUIRED: Select a standard SAM file for PCR duplicate removal")
	return parser.parse_args()

args = get_arguments()
input = args.file
umi_file = args.umi
paired = args.paired

def f(paired): # shows error, quits program if user entered paire-end flag. 
	if paired == True:
		print("Error: this program only works on single-end reads. Goodbye.")
		quit()
	else: 
		pass
f(paired)

umi_set = set()
def get_umi(umi_file): 
	'''Retrieves UMI's from file (specified with -u option) and saves each UMI into a set.'''
	with open(umi_file, "r") as umi_file:
		for line in umi_file:
			umi_set.add(line.strip())
get_umi(umi_file)

# pre-compiled regex for speed optimization
FlagS = re.compile("^\d+S")
FlagS2 = re.compile("\d+S$") # captures last S of CIGAR string for flagRev
FlagI = re.compile("\d+I")
FlagRev = re.compile("\d+[MDN]")
BP = re.compile("\d+") 

def adjust_positionfw(CIGAR, pos): # make sure to add position argument. 
	'''Adjusts leftmost starting position to of FW reads to 5' priming position based on soft clipping.'''
	if 'S' in CIGAR[0:4]: 
		result1 = int(CIGAR.split('S')[0]) # given 10S70M or 70M10S or 10S70M10S ... captures first 10S flag 
		adjusted_pos = int(pos) - result1 # adjust the leftmost position 
		return(adjusted_pos) 
	else: 
		return(pos)

def adjust_positionrv(CIGAR, pos): 
	'''Adjusts leftmost starting position of RV reads to 5' priming position based on CIGAR string.'''
	'''adjusted_pos = pos + M + D + N + last S '''
	result1 = FlagRev.findall(CIGAR)
	result1.extend(FlagS2.findall(CIGAR)) # add 'S' at end of CIGAR string
	result1 = [int(x[:-1]) for x in result1] # remove '[MNDS]', convert numbers to type int
	adjusted_pos = int(pos) + sum(result1) - 1
	return(adjusted_pos) # sum the bp

Dict = defaultdict(dict)
LN = 1
output = input.split('.')[0] + "_deduped.sam"
prev_chr = 0

with open(input, "r") as input, open(output, "a") as output:
	while True:
		alignment = input.readline().strip()
		if alignment == "":
			break # end program at EOF 
		alignment = alignment.split() # aids with ending program at EOF, doesn't split empty space
		if alignment[0].startswith('@'): 
			print("\t".join(alignment), file=output) # save headers to new file
			continue 
		if 'N' in alignment[0][-8:]:
			continue # throws a read away if there is ambiguous nucleotide in the umi
		if LN % 1000 == 0:
			print("De-Duping line ", LN) # status monitor
		chr = str(alignment[2])
		position = str(alignment[3])
		direction = str(alignment[1])
		umi = str(alignment[0][-8:])
		CIGAR = alignment[5]
		if chr != prev_chr: # refresh the dictionary based on chromosome
			Dict = defaultdict(dict)
			prev_chr = chr
		read_info = "_".join([chr, position, direction, umi])
		if umi in umi_set: # if provided umi is in known set of umi 
			
			### SPLIT CONTROL FLOW BY STRAND DIRECTION
			
			if (direction == '0'): # fw reads 
				adjusted_pos = adjust_positionfw(CIGAR, position) # adjust position forst 
				if read_info not in Dict: # if encountering unique read (based on chromosome, adjusted_pos, strand, umi)
					Dict[read_info] = None # save it to a 'ban dictionary' where KEY = 'CHR_ADJUSTEDPOSITION_STRAND_UMI'; VALUE = NONE
					print("\t".join(alignment), file=output)
				else: 
					continue
			if (direction == '16'): # rv reads 
				adjusted_pos = adjust_positionrv(CIGAR, position)
				if read_info not in Dict:
					Dict[read_info] = None 
					print("\t".join(alignment), file=output)
				else: 
					continue
			
		LN += 1

print("PCR Duplicate Removal Finished")