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
Create files containing 16S reads inferred to be from the same gene.
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
import operator
import ntpath

from readConfig import ReadConfig
from seqUtils import extractReads

class ReferenceSeqHit(object):
  def __init__(self, refId):
    self.refId = refId
    self.numHits = 0
    self.pairs = {}
    self.singles = {}

  def addPairHit(self, filename1, filename2, queryId1, queryId2):
    self.numHits += 2

    if filename1 not in self.pairs:
      self.pairs[filename1] = set()

    if filename2 not in self.pairs:
      self.pairs[filename2] = set()

    self.pairs[filename1].add(queryId1)
    self.pairs[filename2].add(queryId2)

  def addSingleHit(self, filename, queryId):
    self.numHits += 1

    if filename not in self.singles:
      self.singles[filename] = set()
    self.singles[filename].add(queryId)

class IdentifyRecoverable16S(object):
  def __init__(self):
    self.ggRefDist = '/srv/db/gg/2013_05/gg_13_5_otus/dist/dist_##_otus.tsv'

    self.bQuiet = False
    pass

  def ggIdFromTaxonomy(self, taxonomy):
    if '(' in taxonomy[7]:
      return taxonomy[7][4:taxonomy[7].rfind('(')]
    else:
      return taxonomy[7][4:].rstrip()

  def readClassifications(self, classificationFile):
    classifications = {}

    for line in open(classificationFile):
      lineSplit = line.split()
      seqId = lineSplit[0]
      taxa = [x.strip() for x in lineSplit[1].split(';') if x.strip() != '']

      classifications[seqId] = taxa

    return classifications

  def identifyConsistentPairs(self, referenceSeqHits, pairFile1, pairFile2, classificationFile, neighbours, bPairsAsSingles, bSingleEnded):
    if not self.bQuiet:
      print '  Reading classification file.'

    classifications = self.readClassifications(classificationFile)

    # get read ids
    if not self.bQuiet:
      print '  Identifying consistent pairs.'

    readIds = set()
    for seqId in classifications:
      readId = seqId
      if '/' in seqId:
        readId = seqId[0:seqId.find('/')]
      readIds.add(readId)

    # find pairs that agree on classification
    numSingletons = 0
    numPairsInAgreement = 0
    numPairsInDisagreement = 0
    numUnclassified = 0
    for readId in readIds:
      ggId1 = ggId2 = None

      seqId1 = readId + '/1'
      seqId2 = readId + '/2'

      if seqId1 in classifications and self.ggIdFromTaxonomy(classifications[seqId1]) != 'unclassified':
        ggId1 = self.ggIdFromTaxonomy(classifications[seqId1])
      else:
        if seqId1 in classifications and self.ggIdFromTaxonomy(classifications[seqId1]) == 'unclassified':
          numUnclassified += 1

        if seqId2 in classifications:
          ggId = self.ggIdFromTaxonomy(classifications[seqId2])
          if ggId != 'unclassified':
            numSingletons += 1
            if bSingleEnded:
              referenceSeqHits[ggId] = referenceSeqHits.get(ggId, ReferenceSeqHit(ggId))
              referenceSeqHits[ggId].addSingleHit(pairFile2, seqId2)
          else:
            numUnclassified += 1

        continue

      if seqId2 in classifications and self.ggIdFromTaxonomy(classifications[seqId2]) != 'unclassified':
        ggId2 = self.ggIdFromTaxonomy(classifications[seqId2])
      else:
        if seqId2 in classifications and self.ggIdFromTaxonomy(classifications[seqId2]) == 'unclassified':
          numUnclassified += 1

        if seqId1 in classifications:
          ggId = self.ggIdFromTaxonomy(classifications[seqId1])
          if ggId != 'unclassified':
            numSingletons += 1
            if bSingleEnded:
              referenceSeqHits[ggId] = referenceSeqHits.get(ggId, ReferenceSeqHit(ggId))
              referenceSeqHits[ggId].addSingleHit(pairFile1, seqId1)
          else:
            numUnclassified += 1

        continue

      if ggId1 == ggId2 or ggId1 in neighbours[ggId2]:
        referenceSeqHits[ggId1] = referenceSeqHits.get(ggId1, ReferenceSeqHit(ggId1))
        referenceSeqHits[ggId1].addPairHit(pairFile1, pairFile2, seqId1, seqId2)
        numPairsInAgreement += 1
      else:
        numPairsInDisagreement += 1
        if bPairsAsSingles:
          referenceSeqHits[ggId1] = referenceSeqHits.get(ggId1, ReferenceSeqHit(ggId1))
          referenceSeqHits[ggId1].addSingleHit(pairFile1, seqId1)
          referenceSeqHits[ggId2] = referenceSeqHits.get(ggId2, ReferenceSeqHit(ggId2))
          referenceSeqHits[ggId2].addSingleHit(pairFile2, seqId2)

    if not self.bQuiet:
      print '    Classified reads: ' + str(len(classifications))
      print '      Singletons: ' + str(numSingletons)
      print '      Reads in pairs with similar classifications: ' + str(2*numPairsInAgreement)
      print '      Reads in pairs with different classifications: ' + str(2*numPairsInDisagreement)
      print '      Unclassified reads: ' + str(numUnclassified)
      print ''

  def addSingletons(self, referenceSeqHits, single, classificationFile):
    if not self.bQuiet:
      print '  Reading classifications.'

    classifications = self.readClassifications(classificationFile)

    if not self.bQuiet:
      print '  Processing single-ended reads.'

    numUnclassified = 0
    for seqId in classifications:
      ggId = self.ggIdFromTaxonomy(classifications[seqId])

      if ggId != 'unclassified':
        referenceSeqHits[ggId] = referenceSeqHits.get(ggId, ReferenceSeqHit(ggId))
        referenceSeqHits[ggId].addSingleHit(single, seqId)
      else:
        numUnclassified += 1

    if not self.bQuiet:
      print '    Classified single-ended reads: ' + str(len(classifications) - numUnclassified)
      print '    Unclassified single-ended reads: ' + str(numUnclassified)
      print ''

  def sortHits(self, referenceSeqHits):
    ggHits = {}
    for ggId in referenceSeqHits:
      ggHits[ggId] = referenceSeqHits[ggId].numHits

    sortedHits = sorted(ggHits.iteritems(), key=operator.itemgetter(1), reverse=True)

    return sortedHits

  def getNeighbours(self, ggRefDistFile, seqIdentityThreshld):
    # determine GG reference genes within sequence identity threshold
    neighbours = {}
    for line in open(ggRefDistFile):
      lineSplit = line.split('\t')
      refId = lineSplit[0].rstrip()
      similarRefSeqs = [x.rstrip() for x in lineSplit[1:]]

      clusteredSeqs = set()
      for item in similarRefSeqs:
        itemSplit = item.split(':')
        if float(itemSplit[1]) < seqIdentityThreshld:
          clusteredSeqs.add(itemSplit[0])

      neighbours[refId] = clusteredSeqs

    return neighbours

  def clusterHits(self, sortedHits, neighbours, minSeqCutoff):
    # clusters reference sequences within sequence identity threshold
    ggClusters = []
    processedIds = set()
    for i in xrange(0, len(sortedHits)):
      ggIdI = sortedHits[i][0]
      if ggIdI in processedIds:
        continue

      processedIds.add(ggIdI)

      clusteredIds = [ggIdI]
      totalHits = sortedHits[i][1]

      if not self.bQuiet:
        reportStr = '  Combining ' + str(ggIdI) + ' with ' + str(totalHits) + ' hits to: '

      for j in xrange(i+1, len(sortedHits)):
        ggIdJ = sortedHits[j][0]
        if ggIdJ in processedIds:
          continue

        if ggIdJ in neighbours[ggIdI]:
          clusteredIds.append(ggIdJ)
          totalHits += sortedHits[j][1]
          processedIds.add(ggIdJ)

          if not self.bQuiet:
            reportStr += str(ggIdJ) + ' (' + str(sortedHits[j][1]) + ' hits) '

      if not self.bQuiet and totalHits > minSeqCutoff:
        print reportStr

      ggClusters.append([totalHits, clusteredIds])

    return ggClusters

  def extractClusteredReads(self, clusters, referenceSeqHits, minSeqCutoff):
    readsInFiles = {}
    for cluster in clusters:
      hits = cluster[0]
      ggIds = cluster[1]

      if hits > minSeqCutoff:
        for ggId in ggIds:
          refSeqHit = referenceSeqHits[ggId]

          for pairFile in refSeqHit.pairs:
            readsInFiles[pairFile] = readsInFiles.get(pairFile, set()).union(refSeqHit.pairs[pairFile])

          for singleFile in refSeqHit.singles:
            readsInFiles[singleFile] = readsInFiles.get(singleFile, set()).union(refSeqHit.singles[singleFile])

    seqsInFiles = {}
    for filename in readsInFiles:
      seqsInFiles[filename] = extractReads(filename, readsInFiles[filename])

    return seqsInFiles

  def extractRecoverable16S(self, referenceSeqHits, neighbours, minSeqCutoff, outputDir):
    # sort hits to GreenGene reference sequences
    if not self.bQuiet:
      print 'Sorting reference sequences by number of hits.'

    sortedHits = self.sortHits(referenceSeqHits)

    if not self.bQuiet:
      print '  Initial clusters: ' + str(len(sortedHits)) + '\n'

    # greedily combine hits to similar GreenGene reference sequences
    if not self.bQuiet:
      print 'Clustering similar reference sequences.'

    ggClusters = self.clusterHits(sortedHits, neighbours, minSeqCutoff)

    if not self.bQuiet:
      print '\n  Final clusters: ' + str(len(ggClusters)) + '\n'

    # get sequences within clusters
    if not self.bQuiet:
      print 'Extracting clustered reads from file.\n'

    seqsInFiles = self.extractClusteredReads(ggClusters, referenceSeqHits, minSeqCutoff)

    # write out clusters containing a sufficient number of hits
    if not self.bQuiet:
      print 'Writing out reads belonging to putative 16S genes: '

    numRecoverable16S = 0
    for cluster in ggClusters:
      hits = cluster[0]
      ggIds = cluster[1]

      if hits > minSeqCutoff:
        numRecoverable16S += 1

        if not self.bQuiet:
          print '  Cluster ' + str(numRecoverable16S) + ': ' + ggIds[0] + ' (' + str(hits) + ' reads)'

        # write out singletons
        singletonsOut = open(outputDir + ggIds[0] + '.singletons.fasta', 'w')
        pairOut1 = open(outputDir + ggIds[0] + '.1.fasta', 'w')
        pairOut2 = open(outputDir + ggIds[0] + '.2.fasta', 'w')
        allOut = open(outputDir + ggIds[0] + '.all.fasta', 'w')

        for ggId in ggIds:
          refSeqHit = referenceSeqHits[ggId]

          for singleFile in refSeqHit.singles:
            for readId in refSeqHit.singles[singleFile]:
              singletonsOut.write('>' + readId + '\n')
              singletonsOut.write(seqsInFiles[singleFile][readId] + '\n')

              allOut.write('>' + readId + '\n')
              allOut.write(seqsInFiles[singleFile][readId] + '\n')

          for pairFile in refSeqHit.pairs:
            for readId in sorted(list(refSeqHit.pairs[pairFile])): # ensure reads are in the same order for both files
              if '/1' in readId:
                pairOut1.write('>' + readId + '\n')
                pairOut1.write(seqsInFiles[pairFile][readId] + '\n')
              elif '/2' in readId:
                pairOut2.write('>' + readId + '\n')
                pairOut2.write(seqsInFiles[pairFile][readId] + '\n')
              else:
                print "[Error] Unrecognized file format. Pairs must be denoted by '/1' and '/2'"
                sys.exit()

              allOut.write('>' + readId + '\n')
              allOut.write(seqsInFiles[pairFile][readId] + '\n')

        singletonsOut.close()
        pairOut1.close()
        pairOut2.close()
        allOut.close()

    if not self.bQuiet:
      print ''
      print '  Number of recoverable 16S genes identified: ' + str(numRecoverable16S)

  def run(self, configFile, otu, seqIdentityThreshold, minSeqCutoff, bPairsAsSingles, bSingleEnded, bQuiet):
    self.bQuiet = bQuiet

    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    ggRefDistFile = self.ggRefDist.replace('##', str(otu))
    neighbours = self.getNeighbours(ggRefDistFile, seqIdentityThreshold)

    # create directory to store putative 16S genes
    dirPutative16S = projectParams['output_dir'] + 'putativeSSU/'
    if not os.path.exists(dirPutative16S):
      os.makedirs(dirPutative16S)
    else:
      rtn = raw_input('Remove previously recovered 16S reads (Y or N)? ')
      if rtn.lower() == 'y' or rtn.lower() == 'yes':
        files = os.listdir(dirPutative16S)
        for f in files:
          if f.endswith('fasta'):
            os.remove(dirPutative16S + '/' + f)
      else:
        sys.exit()

    referenceSeqHits = {}
    for sample in sampleParams:
      if not self.bQuiet:
        print ''
        print sample + ':'

      extractedPrefix = projectParams['output_dir'] + 'extracted/' + sample
      classifiedPrefix = projectParams['output_dir'] + 'classified/' + sample
      pairs = sampleParams[sample]['pairs']
      singles = sampleParams[sample]['singles']

      for i in xrange(0, len(pairs), 2):
        pair1Base = ntpath.basename(pairs[i])
        pair2Base = ntpath.basename(pairs[i+1])
        
        classificationFile = classifiedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.16S.tsv'

        if not self.bQuiet:
          print '  Processing file: ' + classificationFile

        pairFile1 = extractedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.SSU.fasta'
        pairFile2 = extractedPrefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.SSU.fasta'

        self.identifyConsistentPairs(referenceSeqHits, pairFile1, pairFile2, classificationFile, neighbours, bPairsAsSingles, bSingleEnded)
        
        if bSingleEnded:
          classificationFile = classifiedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.16S.tsv'
          
          if not self.bQuiet:
            print '  Processing file: ' + classificationFile
          
          singleFile = extractedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.SSU.fasta'
          self.addSingletons(referenceSeqHits, singleFile, classificationFile)
          
      if bSingleEnded:
        for single in singles:
          singleBase = ntpath.basename(single)
          classificationFile = classifiedPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S.tsv'

          if not self.bQuiet:
            print '  Processing file: ' + classificationFile

          singleFile = extractedPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.SSU.fasta'
          self.addSingletons(referenceSeqHits, singleFile, classificationFile)

    self.extractRecoverable16S(referenceSeqHits, neighbours, minSeqCutoff, dirPutative16S)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Create files containing 16S reads inferred to be from the same gene.")
  parser.add_argument('config_file', help='project config file')
  parser.add_argument('otu', help='GreenGenes reference database used to classify 16S reads (choices: 94, 97, 99)', type=int, choices=[94, 97, 99])
  parser.add_argument('-i', '--seq_identity', help='sequence identity threshold for combining GreenGene hits (default = 0.03)', type=float, default=0.03)
  parser.add_argument('-m', '--min_seqs', help='reads required to create a putative 16S gene file (default = 100).', type=int, default = 100)
  parser.add_argument('-p', '--pairs_as_singles', help='treat paired reads with different classifications as singletons', action="store_true")
  parser.add_argument('-s', '--singletons', help='use single-ended reads and pairs with only one mapped read', action="store_true")
  parser.add_argument('-q', '--quiet', help='suppress output', action='store_true')

  parser.add_argument('--version', help='show version number of program', action='version', version='Identify recoverable 16S v0.0.1')

  args = parser.parse_args()

  identifyRecoverable16S = IdentifyRecoverable16S()
  identifyRecoverable16S.run(args.config_file, args.otu, args.seq_identity, args.min_seqs, args.pairs_as_singles, args.singletons, args.quiet)
