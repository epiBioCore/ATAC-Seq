#!/usr/local/bin/python2.7


import pysam
import sys
import argparse

parser = argparse.ArgumentParser(prog = "shiftBam.py", usage = "%(prog)s --bam in.bam --out out.bam",description = "This program will shift reads aligned to positive strand + 4bps and reads aligned to negative strand - 5bp as required for ATAC -seq")
parser.add_argument("--bam","-b", help="input bam file")
parser.add_argument("--out","-o",help="name of corrected bam file")
args = parser.parse_args()

fname=args.bam
outname=args.out

bamfile=pysam.AlignmentFile(filename=fname,mode='rb')
shifted=pysam.AlignmentFile(filename=outname,mode='wb',template=bamfile)

for read in bamfile:
	if read.is_reverse:
		read.reference_start=read.reference_start-5
		if read.reference_start < 0:
			read.reference_start = 0
		if read.reference_start > bamfile.lengths[read.reference_id]:
			read.reference_start = bamfile.lengths[read.reference_id]
		shifted.write(read)
	else:
		read.reference_start = read.reference_start + 4
		if read.reference_start < 0:
                        read.reference_start = 0
                if read.reference_start > bamfile.lengths[read.reference_id]:
                        read.reference_start = bamfile.lengths[read.reference_id]
		shifted.write(read)

bamfile.close()
shifted.close()
