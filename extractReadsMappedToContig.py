#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


"""
Extract reads mapped to contigs.
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
import pysam

class MapReads(object):
  def __init__(self):
    pass

  def extractReads(self, bamFile, pairedEnd1, pairedEnd2, outputDir):
    read1 = {}
    for line in open(pairedEnd1):
      if line[0] == '>':
        readId = line[1:].split(' ')[0]
        readId = readId.split('/')[0]
        read1[readId] = ''
      else:
        read1[readId] += line.strip()

    read2 = {}
    for line in open(pairedEnd2):
      if line[0] == '>':
        readId = line[1:].split(' ')[0]
        readId = readId.split('/')[0]
        read2[readId] = ''
      else:
        read2[readId] += line.strip()

    bam = pysam.Samfile(bamFile, 'rb')

    # get length of each reference sequences
    for ref in bam.references:
      fout1 = open(outputDir + '/' + ref + '.1.fna', 'w')
      fout2 = open(outputDir + '/' + ref + '.2.fna', 'w')
      for read in bam.fetch(reference = ref):
        if read.tid != read.rnext:
          continue

        print read.qname

        if read.is_read1:
          fout1.write('>' + read.qname + '\n')
          fout1.write(read1[read.qname])
          fout2.write('>' + read.qname + '\n')
          fout2.write(read2[read.qname])

      fout1.close()
      fout2.close()

    bam.close()


  def run(self, seqFile, pairedEnd1, pairedEnd2, threads, outputDir):
    outputFile = seqFile[0:seqFile.rfind('.')]

    print 'Indexing sequence file:'
    os.system('bwa index -a is ' + seqFile)

    print ''
    print 'Mapping reads to sequences:'
    os.system('bwa mem -M -t ' + str(threads) + ' ' + seqFile + ' ' + pairedEnd1 + ' ' + pairedEnd2 + '| samtools view -SubhF 4 - | samtools sort - ' + outputFile)

    print ''
    print 'Indexing BAM file.'
    os.system('samtools index ' + outputFile + '.bam')

    print ''
    print 'Read mapping statistics:'
    os.system('samtools flagstat ' + outputFile + '.bam')

    self.extractReads(outputFile + '.bam', pairedEnd1, pairedEnd2, outputDir)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Extract reads mapped to contigs.")
  parser.add_argument('seq_file', help='Sequence file in FASTA format.')
  parser.add_argument('pair1', help='First paired-end read file in FASTA format.')
  parser.add_argument('pair2', help='Second paired-end read file in FASTA format.')
  parser.add_argument('output_dir', help='Output directory.')
  parser.add_argument('-t', '--threads', help='Number of threads.', type=int)

  args = parser.parse_args()

  mapReads = MapReads()
  mapReads.run(args.seq_file, args.pair1, args.pair2, args.threads, args.output_dir)
