#!/usr/bin/env python

"""
Identify 16S sequences within bins.
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

  def readHitTable(self, resultsFile):
    seqIds = set()
    for line in open(resultsFile):
      if line.strip() == '' or line[0] == '#':
        continue

      seqIds.add(line.split()[0])

    return seqIds

  def readHits(self, resultsFile, domain, revComp = False):
    seqInfo = {}

    bReadHit = False
    for line in open(resultsFile):
      if line[0:2] == '>>':
        lineSplit = line.split()
        seqId = lineSplit[1]
        bReadHit = True
        hitCounter = 0
        continue
      elif bReadHit:
        hitCounter += 1
        if hitCounter == 3:
          lineSplit = line.split()

          iEvalue = lineSplit[5]
          alignFrom = lineSplit[9]
          alignTo = lineSplit[10]
          alignLen = int(alignTo) - int(alignFrom) + 1

          seqInfo[seqId]  = [domain, iEvalue, alignFrom, alignTo, str(alignLen), str(revComp)]
          bReadHit = False

    return seqInfo

  def addHit(self, hits, seqId, info):
    if seqId not in hits:
      hits[seqId] = info
      return

    index = 1
    while(True):
      # if hits overlap then retain only the longest
      startNew = int(info[2])
      endNew = int(info[3])
      lengthNew = int(info[4])

      start = int(hits[seqId][2])
      end = int(hits[seqId][3])
      length = int(hits[seqId][4])

      if (start >= startNew and start <= endNew) or (startNew >= start and startNew <= endNew):
        if lengthNew > length:
          hits[seqId] = info
        break
      else:
        # multiple hits on the same contig/scaffold
        index += 1
        newSeqId = seqId + '-#' + str(index)
        if newSeqId in hits:
          seqId = newSeqId # see if the new hit overlaps with this
        else:
          hits[newSeqId] = info
          break

  def run(self, contigFile, binDir, outputDir, extension, threads, evalue):
    # make sure output directory exists
    if not os.path.exists(outputDir):
      os.makedirs(outputDir)

    # get bin id of binned contigs
    seqIdToBinId = {}
    files = os.listdir(binDir)
    for file in files:
      if file.endswith(extension):
        binId = file[0:file.rfind('.')]
        for line in open(binDir + '/' + file):
          if line[0] == '>':
            seqId = line[1:].split(' ')[0].strip()
            seqIdToBinId[seqId] = binId

    # identify 16S reads from contigs/scaffolds
    print 'Identifying 16S genes on assembled contigs/scaffolds.'
    extract16S = Extract16S()
    extract16S.hmmSearch(contigFile, threads, evalue, outputDir + '/identified16S')

    # read HMM hits
    print 'Parsing HMM results.'
    hits = {}
    for domain in ['archeae', 'bacteria', 'euk']:
      # forward hits
      seqIds = self.readHitTable(outputDir + '/identified16S' + '.' + domain + '.table.txt')
      if len(seqIds) > 0:
        seqInfo = self.readHits(outputDir + '/identified16S' + '.' + domain + '.txt', domain)
        for seqId, info in seqInfo.iteritems():
          self.addHit(hits, seqId, info)

      # reverse complement hits
      seqIds = self.readHitTable(outputDir + '/identified16S' + '.' + domain + '.table.rev_comp.txt')
      if len(seqIds) > 0:
        seqInfo = self.readHits(outputDir + '/identified16S' + '.' + domain + '.rev_comp.txt', domain, True)
        for seqId, info in seqInfo.iteritems():
          self.addHit(hits, seqId, info)

    # write summary file and putative 16S genes to file
    print 'Writing results to file.'
    summaryFile = outputDir + '/identified16S.tsv'
    summaryOut = open(summaryFile, 'w')
    summaryOut.write('Bin Id\tSeq. Id\tHMM model\ti-Evalue\tStart hit\tEnd hit\t16S/18S gene length\tRev. Complement\tContig/Scaffold length\n\n')

    seqFile = outputDir + '/identified16S.fna'
    seqOut = open(seqFile, 'w')

    seqs = self.readFasta(contigFile)

    hitsToBins = {}
    for seqId in hits:
      origSeqId = seqId
      if '-#' in seqId:
        seqId = seqId[0:seqId.rfind('-#')]

      if seqId in seqIdToBinId:
        binId = seqIdToBinId[seqId]
      else:
        binId = 'Unbinned'

      seqInfo = [origSeqId] + hits[origSeqId]
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
    print 'Identified ' + str(len(hits)) + ' putative 16S genes:'
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
  parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')

  parser.add_argument('--version', help='show version number of program', action='version', version='Extract 16S using HMMs v0.0.1')

  args = parser.parse_args()

  identifyBinned16S = IdentifyBinned16S()
  identifyBinned16S.run(args.contig_file, args.bin_dir, args.output_dir, args.extension, args.threads, args.evalue)
