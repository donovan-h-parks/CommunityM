#!/usr/bin/env python

"""
Produce FASTQ file that mirrors reads in a single FASTA file.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os, sys, argparse

class MirrorSeqs(object):
  def __init__(self):
    pass

  def extractSeqsFromFastq(self, fastqFile, seqsToExtract, fout):
    lineNum = 0
    for line in open(fastqFile):
      if lineNum == 0:
        seqId = line[1:line.find(' ')].strip()
        bWrite = (seqId in seqsToExtract)

      if bWrite:
        fout.write(line)

      lineNum += 1
      lineNum = lineNum % 4

  def run(self, fastaFile, pairedEnd1, pairedEnd2):
    # extract sequences from FASTA file
    seqIds1 = set()
    seqIds2 = set()
    for line in open(fastaFile):
      if line[0] == '>':
        seqId = line[1:line.find(' ')].strip()
        if '/1' in seqId:
          seqIds1.add(seqId)
        elif '/2' in seqId:
          seqIds2.add(seqId)
        else:
          print "[Error] Invalid file format. Pairs must be specified by a '/1' and '/2'."
          sys.exit()

    # produce FASTQ file with reads in FASTA file
    fout = open(fastaFile[0:fastaFile.rfind('.')] + '.fq', 'w')
    self.extractSeqsFromFastq(pairedEnd1, seqIds1, fout)
    self.extractSeqsFromFastq(pairedEnd2, seqIds2, fout)
    fout.close()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Produce FASTQ file that mirrors reads in a single FASTA file.")
  parser.add_argument('fasta_file', help='FASTA file containing sequences to mirror.')
  parser.add_argument('pair1', help='First paired-end read file in FASTQ format.')
  parser.add_argument('pair2', help='Second paired-end read file in FASTQ format.')

  args = parser.parse_args()

  mirrorSeqs = MirrorSeqs()
  mirrorSeqs.run(args.fasta_file, args.pair1, args.pair2)
