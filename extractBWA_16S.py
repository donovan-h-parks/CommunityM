#!/usr/bin/env python

"""
Extract reads that mapped to the GreenGenes database using BWA.
"""

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2012"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import sys
import os
import argparse
import time
import ntpath

from readConfig import ReadConfig
from seqUtils import extractSeqs

import pysam
from numpy import mean, std

class ExtractBWA(object):
  def __init__(self):
    self.mappingQualityThreshold = 10
    self.minLength = 90
    pass

  def readPairedBAM(self, bamFile):
    # read compressed BAM file and report basic statistics
    bam = pysam.Samfile(bamFile, 'rb')

    print 'Reference 16S sequences: ' + str(bam.nreferences)
    print 'Mapped reads: ' + str(bam.mapped)

    # find all reads that mapped to a 16S sequence
    readsMappedTo16S_1 = set()
    readsMappedTo16S_2 = set()
    numMultipleTopHits = 0
    for read in bam.fetch(until_eof=True):
      if read.is_unmapped or read.is_secondary or (read.mapq < self.mappingQualityThreshold):
        continue

      if read.alen < self.minLength:
        continue

      for tag in read.tags:
        if tag[0] == 'X0' and int(tag[1]) > 1:
          numMultipleTopHits += 1

      if read.is_read1:
        readsMappedTo16S_1.add(read.qname)
      elif read.is_read2:
        readsMappedTo16S_2.add(read.qname)
      else:
        print '[Error] Read is neither read 1 or read 2.'
        sys.exit()

    bam.close()

    print '  Reads with multiple top hits: ' + str(numMultipleTopHits)

    return readsMappedTo16S_1, readsMappedTo16S_2

  def readSingleBAM(self, bamFile):
    # read compressed BAM file and report basic statistics
    bam = pysam.Samfile(bamFile, 'rb')

    print 'Reference 16S sequences: ' + str(bam.nreferences)
    print 'Mapped reads: ' + str(bam.mapped)

    # find all reads that mapped to a 16S sequence
    readsMappedTo16S = set()
    numMultipleTopHits = 0
    for read in bam.fetch(until_eof=True):
      if read.is_unmapped or read.is_secondary or (read.mapq <= self.mappingQualityThreshold):
        continue

      if read.alen < self.minLength:
        continue

      for tag in read.tags:
        if tag[0] == 'X0' and int(tag[1]) > 1:
          numMultipleTopHits += 1

      readsMappedTo16S.add(read.qname)

    bam.close()

    print '  Reads with multiple top hits: ' + str(numMultipleTopHits)

    return readsMappedTo16S

  def processPairs(self, pairs, outputDir, prefix):
    for i in xrange(0, len(pairs), 2):
      pair1 = pairs[i]
      pair2 = pairs[i+1]

      print 'Identifying 16S sequences in paired-end reads: ' + pair1 + ', ' + pair2

      outputPrefix = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')]
      outputPrefix1 = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')] + '.1'
      outputPrefix2 = prefix + '.' + pair2[pair2.rfind('/')+1:pair2.rfind('.')] + '.2'

      readsMappedTo16S_1, readsMappedTo16S_2  = self.readPairedBAM(outputDir + ntpath.basename(pair1) + '.bam')

      print '  Hits in ' + pair1 + ': ' + str(len(readsMappedTo16S_1))
      print '  Hits in ' + pair2 + ': ' + str(len(readsMappedTo16S_2))
      print ''

      # extract all pairs where at least one read is on a 16S sequence
      print 'Extracting putative 16S reads: '
      readsMappedTo16S = readsMappedTo16S_1.union(readsMappedTo16S_2)

      seqs1 = extractSeqs(pair1, readsMappedTo16S)
      seqs2 = extractSeqs(pair2, readsMappedTo16S)

      # create file with all 16S sequences
      allSeqFile = outputPrefix + '.all.16S.fasta'
      fout = open(allSeqFile, 'w')
      for seqId in readsMappedTo16S_1:
        if seqs1[seqId][0].endswith('/1'):
          fout.write('>' + seqs1[seqId][0] + '\n')
        else:
          fout.write('>' + seqs1[seqId][0] + '/1\n')
        fout.write(seqs1[seqId][1] + '\n')

      for seqId in readsMappedTo16S_2:
        if seqs2[seqId][0].endswith('/1'):
          fout.write('>' + seqs2[seqId][0] + '\n')
        else:
          fout.write('>' + seqs2[seqId][0] + '/2\n')
        fout.write(seqs2[seqId][1] + '\n')

      fout.close()

      print '  All identified 16S reads written to: ' + allSeqFile

      # create paired-end files where at least one read maps to a 16S
      pair1File = outputPrefix1 + '.16S.fasta'
      fout = open(pair1File, 'w')
      for seqId in readsMappedTo16S:
        if seqs1[seqId][0].endswith('/1'):
          fout.write('>' + seqs1[seqId][0] + '\n')
        else:
          fout.write('>' + seqs1[seqId][0] + '/1\n')
        fout.write(seqs1[seqId][1] + '\n')
      fout.close()

      pair2File = outputPrefix2 + '.16S.fasta'
      fout = open(pair2File, 'w')
      for seqId in readsMappedTo16S:
        if seqs2[seqId][0].endswith('/1'):
          fout.write('>' + seqs2[seqId][0] + '\n')
        else:
          fout.write('>' + seqs2[seqId][0] + '/2\n')
        fout.write(seqs2[seqId][1] + '\n')
      fout.close()

      print '  Pairs with at least one read identified as 16S written to: ' + pair1File + ', ' + pair2File
      print '    Pairs with at least one read identified as 16S: ' + str(len(readsMappedTo16S))
      print ''

  def processSingles(self, singles, outputDir, prefix):
    for i in xrange(0, len(singles)):
      seqFile = singles[i]

      print 'Identifying 16S sequences in single-end reads: ' + seqFile

      outputPrefix = prefix + '.' + seqFile[seqFile.rfind('/')+1:seqFile.rfind('.')]

      readsMappedTo16S = self.readSingleBAM(outputDir + ntpath.basename(seqFile) + '.bam')

      print '  Hits in ' + seqFile + ': ' + str(len(readsMappedTo16S))

      # extract reads with hits
      seqs = extractSeqs(seqFile, readsMappedTo16S)

      # create file with all 16S sequences
      allSeqFile = outputPrefix + '.all.16S.fasta'
      fout = open(allSeqFile, 'w')
      for seqId in readsMappedTo16S:
        fout.write('>' + seqs[seqId][0] + '\n')
        fout.write(seqs[seqId][1] + '\n')

      fout.close()

      print '  Identified 16S reads written to: ' + allSeqFile
      print ''

  def run(self, configFile, mappingQual, minLength):
    self.mappingQualityThreshold = mappingQual
    self.minLength = minLength

    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    for sample in sampleParams:
      outputDir = projectParams['output_dir']
      prefix = outputDir + sample
      pairs = sampleParams[sample]['pairs']
      singles = sampleParams[sample]['singles']

      # identify 16S sequences in paired-end reads
      self.processPairs(pairs, outputDir, prefix)

      # identify 16S sequences in single-end reads
      self.processSingles(singles, outputDir, prefix)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Extract reads that mapped to the GreenGenes database using BWA.")

  parser.add_argument('config_file', help='project config file')
  parser.add_argument('-m','--mapping_qual', help='mapping quality required to consider read (default = 10).', type=int, default=10)
  parser.add_argument('-l','--min_length', help='minimum alignment length required to retain read (default = 90).', type=int, default=90)

  parser.add_argument('--version', help='show version number of program', action='version', version='Extract 16S using BWA v0.0.1')

  args = parser.parse_args()

  extractBWA = ExtractBWA()
  extractBWA.run(args.config_file, args.mapping_qual, args.min_length)
