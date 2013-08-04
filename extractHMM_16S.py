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
Extract 16S sequences from metagenomic data using HMMs.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse

from readConfig import ReadConfig
from seqUtils import extractSeqs

class Extract16S(object):
  def __init__(self):
    self.bacteriaModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_bacteria.hmm'
    self.bacteriaRevCompModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_bacteria.revComp.hmm'

    self.archaeaModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_archaea.hmm'
    self.archaeaRevCompModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_archaea.revComp.hmm'

    self.eukModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_euk.hmm'
    self.eukRevCompModelFile = '/srv/whitlam/bio/apps/12.04/sw/communitym/dev/models/SSU_euk.revComp.hmm'

    pass

  def getHits(self, hitTable):
    seqIds = set()
    for line in open(hitTable):
      if line[0] == '#' or line.strip() == '':
        continue

      seqId = line.split()[0].split('/')[0]   # remove read identifier
      seqIds.add(seqId)

    return seqIds

  def hmmSearch(self, seqFile, threads, evalue, outputPrefix):
    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.bacteria.txt --tblout ' + outputPrefix + '.bacteria.table.txt -E ' + evalue + ' ' + self.bacteriaModelFile + ' ' + seqFile)
    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.bacteria.rev_comp.txt --tblout ' + outputPrefix + '.bacteria.table.rev_comp.txt -E ' + evalue + ' ' + self.bacteriaRevCompModelFile + ' ' + seqFile)

    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.archeae.txt --tblout ' + outputPrefix + '.archeae.table.txt -E ' + evalue + ' ' + self.archaeaModelFile + ' ' + seqFile)
    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.archeae.rev_comp.txt --tblout ' + outputPrefix + '.archeae.table.rev_comp.txt -E ' + evalue + ' ' + self.archaeaRevCompModelFile + ' ' + seqFile)

    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.euk.txt --tblout ' + outputPrefix + '.euk.table.txt -E ' + evalue + ' ' + self.eukModelFile + ' ' + seqFile)
    rtn = os.system('hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.euk.rev_comp.txt --tblout ' + outputPrefix + '.euk.table.rev_comp.txt -E ' + evalue + ' ' + self.eukRevCompModelFile + ' ' + seqFile)

  def processPairs(self, pairs, threads, evalue, prefix, bQuiet):
    for i in xrange(0, len(pairs), 2):
      pair1 = pairs[i]
      pair2 = pairs[i+1]

      if not bQuiet:
        print 'Identifying 16S sequences in paired-end reads: ' + pair1 + ', ' + pair2

      outputPrefix = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')]
      outputPrefix1 = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')] + '.1'
      outputPrefix2 = prefix + '.' + pair2[pair2.rfind('/')+1:pair2.rfind('.')] + '.2'

      self.hmmSearch(pair1, threads, evalue, outputPrefix1)
      self.hmmSearch(pair2, threads, evalue, outputPrefix2)

      # reads hits
      hitsBacteria1 = self.getHits(outputPrefix1 + '.bacteria.table.txt')
      hitsRevCompBacteria1 = self.getHits(outputPrefix1 + '.bacteria.table.rev_comp.txt')

      hitsArchaea1 = self.getHits(outputPrefix1 + '.archeae.table.txt')
      hitsRevCompArcheae1 = self.getHits(outputPrefix1 + '.archeae.table.rev_comp.txt')

      hitsEuk1 = self.getHits(outputPrefix1 + '.euk.table.txt')
      hitsRevCompEuk1 = self.getHits(outputPrefix1 + '.euk.table.rev_comp.txt')

      hitsBacteria2 = self.getHits(outputPrefix2 + '.bacteria.table.txt')
      hitsRevCompBacteria2 = self.getHits(outputPrefix2 + '.bacteria.table.rev_comp.txt')

      hitsArchaea2 = self.getHits(outputPrefix2 + '.archeae.table.txt')
      hitsRevCompArcheae2 = self.getHits(outputPrefix2 + '.archeae.table.rev_comp.txt')

      hitsEuk2 = self.getHits(outputPrefix2 + '.euk.table.txt')
      hitsRevCompEuk2 = self.getHits(outputPrefix2 + '.euk.table.rev_comp.txt')

      # combine hits
      hits1 = hitsBacteria1.union(hitsRevCompBacteria1).union(hitsArchaea1).union(hitsRevCompArcheae1).union(hitsEuk1).union(hitsRevCompEuk1)
      hits2 = hitsBacteria2.union(hitsRevCompBacteria2).union(hitsArchaea2).union(hitsRevCompArcheae2).union(hitsEuk2).union(hitsRevCompEuk2)

      if not bQuiet:
        print '  Hits in ' + pair1 + ': ' + str(len(hits1))
        print '    Fwd. bacterial hits: ' + str(len(hitsBacteria1))
        print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria1))
        print '    Fwd. archeael hits: ' + str(len(hitsArchaea1))
        print '    Rev. comp. archeael hits: ' + str(len(hitsRevCompArcheae1))
        print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk1))
        print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk1))
        print ''

        print '  Hits in ' + pair2 + ': ' + str(len(hits2))
        print '    Fwd. bacterial hits: ' + str(len(hitsBacteria2))
        print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria2))
        print '    Fwd. archeael hits: ' + str(len(hitsArchaea2))
        print '    Rev. comp. archeael hits: ' + str(len(hitsRevCompArcheae2))
        print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk2))
        print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk2))
        print ''

      # extract reads with hits
      if not bQuiet:
        print 'Extracting putative 16S reads:'
      hitUnion = hits1.union(hits2)

      seqs1 = extractSeqs(pair1, hitUnion)
      seqs2 = extractSeqs(pair2, hitUnion)

      # create file with all 16S sequences
      allSeqFile = outputPrefix + '.all.16S.fasta'
      fout = open(allSeqFile, 'w')
      for seqId in hits1:
        fout.write('>' + seqs1[seqId][0] + '\n')
        fout.write(seqs1[seqId][1] + '\n')

      for seqId in hits2:
        fout.write('>' + seqs2[seqId][0] + '\n')
        fout.write(seqs2[seqId][1] + '\n')

      fout.close()

      # create paired-end files where at least one read maps to a 16S
      pair1FileUnion = outputPrefix1 + '.union.16S.fasta'
      fout = open(pair1FileUnion, 'w')
      for seqId in hitUnion:
        fout.write('>' + seqs1[seqId][0] + '\n')
        fout.write(seqs1[seqId][1] + '\n')
      fout.close()

      pair2FileUnion = outputPrefix2 + '.union.16S.fasta'
      fout = open(pair2FileUnion, 'w')
      for seqId in hitUnion:
        fout.write('>' + seqs2[seqId][0] + '\n')
        fout.write(seqs2[seqId][1] + '\n')
      fout.close()

      # create paired-end files where at least one read maps to a 16S
      hitIntersection = hits1.intersection(hits2)
      pair1FileIntersect = outputPrefix1 + '.intersect.16S.fasta'
      fout = open(pair1FileIntersect, 'w')
      for seqId in hitIntersection:
        fout.write('>' + seqs1[seqId][0] + '\n')
        fout.write(seqs1[seqId][1] + '\n')
      fout.close()

      pair2FileIntersect = outputPrefix2 + '.intersect.16S.fasta'
      fout = open(pair2FileIntersect, 'w')
      for seqId in hitIntersection:
        fout.write('>' + seqs2[seqId][0] + '\n')
        fout.write(seqs2[seqId][1] + '\n')
      fout.close()

      if not bQuiet:
        print '  Hits to left reads: ' + str(len(hits1))
        print '  Hits to right reads: ' + str(len(hits2)) + ' reads'
        print '  Pairs with at least one read identified as 16S: ' + str(len(hitUnion)) + ' reads'
        print '  Pairs with both read identified as 16S: ' + str(len(hitIntersection)) + ' reads'
        print ''
        print '  All identified 16S reads: ' + allSeqFile
        print '  Pairs with at least one read identified as 16S: '
        print '      ' + pair1FileUnion
        print '      ' + pair2FileUnion
        print '  Pairs with both read identified as 16S: '
        print '      ' + pair1FileIntersect
        print '      ' + pair2FileIntersect

  def processSingles(self, singles, threads, evalue, prefix, bQuiet):
    for i in xrange(0, len(singles)):
      seqFile = singles[i]

      if not bQuiet:
        print 'Identifying 16S sequences in single-end reads: ' + seqFile

      outputPrefix = prefix + '.' + seqFile[seqFile.rfind('/')+1:seqFile.rfind('.')]

      self.hmmSearch(seqFile, threads, evalue, outputPrefix)

      # reads hits
      hitsBacteria = self.getHits(outputPrefix + '.bacteria.table.txt')
      hitsRevCompBacteria = self.getHits(outputPrefix + '.bacteria.table.rev_comp.txt')

      hitsArchaea = self.getHits(outputPrefix + '.archeae.table.txt')
      hitsRevCompArcheae = self.getHits(outputPrefix + '.archeae.table.rev_comp.txt')

      hitsEuk = self.getHits(outputPrefix + '.euk.table.txt')
      hitsRevCompEuk = self.getHits(outputPrefix + '.euk.table.rev_comp.txt')

      hits = hitsBacteria.union(hitsRevCompBacteria).union(hitsArchaea).union(hitsRevCompArcheae).union(hitsEuk).union(hitsRevCompEuk)

      if not bQuiet:
        print '  Hits in ' + seqFile + ': ' + str(len(hits))
        print '    Fwd. bacterial hits: ' + str(len(hitsBacteria))
        print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria))
        print '    Fwd. archeael hits: ' + str(len(hitsArchaea))
        print '    Rev. comp. archeael hits: ' + str(len(hitsRevCompArcheae))
        print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk))
        print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk))

      # extract reads with hits
      seqs = extractSeqs(seqFile, hits)

      # create file with all 16S sequences
      allSeqFile = outputPrefix + '.all.16S.fasta'
      fout = open(allSeqFile, 'w')
      for seqId in hits:
        fout.write('>' + seqs[seqId][0] + '\n')
        fout.write(seqs[seqId][1] + '\n')

      fout.close()

      if not bQuiet:
        print '  Identified 16S reads written to: ' + allSeqFile
        print ''

  def run(self, configFile, threads, evalue, bQuiet):
    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = False)

    for sample in sampleParams:
      prefix = projectParams['output_dir'] + sample
      pairs = sampleParams[sample]['pairs']
      singles = sampleParams[sample]['singles']

      # identify 16S sequences in paired-end reads
      self.processPairs(pairs, threads, evalue, prefix, bQuiet)

      # identify 16S sequences in single-end reads
      self.processSingles(singles, threads, evalue, prefix, bQuiet)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Extract 16S sequences from metagenomic data using HMMs.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('config_file', help='project config file.')

  parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)
  parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')
  parser.add_argument('-q', '--quiet', help='Surpress all output', action='store_true')

  parser.add_argument('--version', help='Show version number of program', action='version', version='Extract 16S using HMMs v0.0.1')

  args = parser.parse_args()

  extract16S = Extract16S()
  extract16S.run(args.config_file, args.threads, args.evalue, args.quiet)
