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
Identify 16S sequences within bins.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import argparse

from extractHMM_16S import Extract16S

class IdentifyBinned16S(object):
  def __init__(self):
    pass

  def readFasta(self, fastaFile):
    try:
      fh = file(fastaFile)
    except IOError:
      print "File '" + fastaFile + "' does not exist."
      sys.exit()

    seqs = {}
    for line in fh:
      if line.startswith('>'):
        seqId = line[1:].split()[0]
        seqs[seqId] = ''
      else:
        seqs[seqId] += line.rstrip('\n').rstrip('*')

    return seqs

  def readHits(self, resultsFile, domain, evalueThreshold, revComp = False):
    seqInfo = {}

    bReadHit = False
    for line in open(resultsFile):
      if line[0:2] == '>>':
        lineSplit = line.split()
        seqId = lineSplit[1]
        bReadHit = True
        hitCounter = 0
        continue
      elif line.strip() == '':
        bReadHit = False
      elif bReadHit:
        hitCounter += 1
        if hitCounter >= 3:
          lineSplit = line.split()

          iEvalue = lineSplit[5]
          aliFrom = lineSplit[9]
          aliTo = lineSplit[10]
          alignLen = int(aliTo) - int(aliFrom) + 1

          if float(iEvalue) <= evalueThreshold:
            seqInfo[seqId] = seqInfo.get(seqId, []) + [[domain, iEvalue, aliFrom, aliTo, str(alignLen), str(revComp)]]

    return seqInfo

  def addHit(self, hits, seqId, info, concatenateThreshold):
    if seqId not in hits:
      hits[seqId] = info
      return

    baseSeqId = seqId
    index = 1
    bConcatenate = False
    concateSeqId = seqId
    while(True):
      # if hits overlap then retain only the longest
      startNew = int(info[2])
      endNew = int(info[3])
      lengthNew = int(info[4])

      start = int(hits[seqId][2])
      end = int(hits[seqId][3])
      length = int(hits[seqId][4])

      # check if hits should be concatenated
      if abs(start - endNew) < concatenateThreshold:
        # new hit closely preceded old hit
        del hits[seqId]
        info[2] = str(startNew)
        info[3] = str(end)
        info[4] = str(end - startNew + 1)
        hits[concateSeqId] = info
        bConcatenate = True
      elif abs(startNew - end) < concatenateThreshold:
        # new hit closely follows old hit
        del hits[seqId]
        info[2] = str(start)
        info[3] = str(endNew)
        info[4] = str(endNew - start + 1)
        hits[concateSeqId] = info
        bConcatenate = True

      index += 1
      newSeqId = baseSeqId + '-#' + str(index)
      if bConcatenate:
        if newSeqId in hits:
          seqId = newSeqId # see if other sequences concatenate
        else:
          break
      else:
        # hits are not close enough to concatenate
        if newSeqId in hits:
          seqId = newSeqId # see if the new hit overlaps with this
          concateSeqId = newSeqId
        else:
          hits[newSeqId] = info
          break

  def addDomainHit(self, hits, seqId, info):
    if seqId not in hits:
      hits[seqId] = info
      return

    baseSeqId = seqId
    overlapSeqId = seqId

    index = 1
    bOverlap = False
    while(True):
      # if hits overlap then retain only the longest
      startNew = int(info[2])
      endNew = int(info[3])
      lengthNew = int(info[4])

      start = int(hits[seqId][2])
      end = int(hits[seqId][3])
      length = int(hits[seqId][4])

      if (startNew <= start and endNew >= start) or (start <= startNew and end >= startNew):
        bOverlap = True

        if lengthNew > length:
          hits[overlapSeqId] = info
        else:
          hits[overlapSeqId] = hits[seqId]

        if overlapSeqId != seqId:
          del hits[seqId]

      index += 1
      newSeqId = baseSeqId + '-#' + str(index)
      if newSeqId in hits:
        seqId = newSeqId # see if the new hit overlaps with this
        if not bOverlap:
          overlapSeqId = seqId
      else:
        break

    if not bOverlap:
      hits[newSeqId] = info

  def run(self, contigFile, binDir, outputDir, extension, threads, evalueThreshold, concatenateThreshold):
    # make sure output directory exists
    if not os.path.exists(outputDir):
      os.makedirs(outputDir)

    # get bin id of binned contigs
    seqIdToBinId = {}
    files = os.listdir(binDir)
    for f in files:
      if f.endswith(extension):
        binId = f[0:f.rfind('.')]
        for line in open(binDir + '/' + f):
          if line[0] == '>':
            seqId = line[1:].split()[0].strip()
            seqIdToBinId[seqId] = binId

    if len(seqIdToBinId) == 0:
      print '[Error] No contigs/scaffolds identified in bins. Check the extension of used to identify genome bins.'
      sys.exit()

    # identify 16S reads from contigs/scaffolds
    print 'Identifying 16S genes on assembled contigs/scaffolds.'
    extract16S = Extract16S()
    extract16S.hmmSearch(contigFile, threads, evalueThreshold, outputDir + '/identified16S')

    # read HMM hits
    print 'Parsing HMM results.'
    hitsPerDomain = {}
    for domain in ['archaea', 'bacteria', 'euk']:
      hits = {}

      # forward hits
      seqInfo = self.readHits(outputDir + '/identified16S' + '.' + domain + '.txt', domain, evalueThreshold)
      if len(seqInfo) > 0:
        for seqId, seqHits in seqInfo.iteritems():
          for hit in seqHits:
            self.addHit(hits, seqId, hit, concatenateThreshold)

      # reverse complement hits
      seqInfo = self.readHits(outputDir + '/identified16S' + '.' + domain + '.rev_comp.txt', domain, evalueThreshold, True)
      if len(seqInfo) > 0:
        for seqId, seqHits in seqInfo.iteritems():
          for hit in seqHits:
            self.addHit(hits, seqId, hit, concatenateThreshold)

      hitsPerDomain[domain] = hits

    # find best domain hit for each sequence
    bestHits = {}
    for _, hits in hitsPerDomain.iteritems():
      for seqId, info in hits.iteritems():
        if '-#' in seqId:
          seqId = seqId[0:seqId.rfind('-#')]

        self.addDomainHit(bestHits, seqId, info)

    # write summary file and putative 16S genes to file
    print 'Writing results to file.'
    summaryFile = outputDir + '/identified16S.tsv'
    summaryOut = open(summaryFile, 'w')
    summaryOut.write('Bin Id\tSeq. Id\tHMM model\ti-Evalue\tStart hit\tEnd hit\t16S/18S gene length\tRev. Complement\tContig/Scaffold length\n')

    seqFile = outputDir + '/identified16S.fna'
    seqOut = open(seqFile, 'w')

    seqs = self.readFasta(contigFile)

    hitsToBins = {}
    for seqId in bestHits:
      origSeqId = seqId
      if '-#' in seqId:
        seqId = seqId[0:seqId.rfind('-#')]

      if seqId in seqIdToBinId:
        binId = seqIdToBinId[seqId]
      else:
        binId = 'Unbinned'

      seqInfo = [origSeqId] + bestHits[origSeqId]
      hitsToBins[binId] = hitsToBins.get(binId, []) + [seqInfo]

    for binId in sorted(hitsToBins.keys()):
      for seqInfo in hitsToBins[binId]:
        seqId = seqInfo[0]
        if '-#' in seqId:
          seqId = seqId[0:seqId.rfind('-#')]

        seq = seqs[seqId]
        summaryOut.write(binId + '\t' + '\t'.join(seqInfo) + '\t' + str(len(seq)) + '\n')
        seqOut.write('>' + binId + '_' + seqInfo[0] + '\n')
        seqOut.write(seq[int(seqInfo[3]):int(seqInfo[4])] + '\n')

    summaryOut.close()
    seqOut.close()

    print ''
    print 'Identified ' + str(len(bestHits)) + ' putative 16S genes:'
    print '  Summary of identified hits written to: ' + summaryFile
    print '  Putative 16S sequences written to: ' + seqFile

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Identify 16S sequences within bins.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('contig_file', help='FASTA file of assembled contigs/scaffolds')
  parser.add_argument('bin_dir', help='directory containing bin')
  parser.add_argument('output_dir', help='output directory')

  parser.add_argument('-x', '--extension', help='extension of bins', default = 'fna')
  parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)
  parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', type=float, default = 1e-5)
  parser.add_argument('-c', '--concatenate', help='concatenate hits that are within the specified number of base pairs', type=int, default = 100)

  parser.add_argument('--version', help='show version number of program', action='version', version='Extract 16S using HMMs v0.0.1')

  args = parser.parse_args()

  identifyBinned16S = IdentifyBinned16S()
  identifyBinned16S.run(args.contig_file, args.bin_dir, args.output_dir, args.extension, args.threads, args.evalue, args.concatenate)
