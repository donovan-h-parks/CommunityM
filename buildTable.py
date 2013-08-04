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
Build table in BIOM format summarizing classified 16S sequences.
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

from biom.table import table_factory, SparseOTUTable
import numpy

from readConfig import ReadConfig

class BuildTable(object):
  def __init__(self):
    self.ranksByLabel = {'Domain':0, 'Phylum':1, 'Class':2, 'Order':3, 'Family':4, 'Genus':5, 'Species':6, 'GG_ID':7}
    self.ranksByLevel = {0:'Domain', 1:'Phylum', 2:'Class', 3:'Order', 4:'Family', 5:'Genus', 6:'Species', 7:'GG_ID'}
    pass

  def parseClassificationFile(self, classificationFile, bootstrapThreshold, rankIndex, counts, taxonomy):
    with open(classificationFile) as fin:
      for line in fin:
        lineSplit = line.split('\t')
        taxa = lineSplit[1].split(';')

        if '(' in taxa[rankIndex]:
          taxonSplit = taxa[rankIndex].split('(')
          taxonId = taxonSplit[0]
          taxonId = taxonId[taxonId.find('__')+2:].strip()
          bootstrapSupport = int(taxonSplit[1][0:taxonSplit[1].find(')')])
        else:
          taxonId = taxa[rankIndex].strip()
          bootstrapSupport = 0

        if bootstrapSupport >= bootstrapThreshold:
          counts[taxonId] = counts.get(taxonId, 0) + 1

          taxaDict = {}
          taxaList = []
          for t in xrange(0, rankIndex+1):
            if '(' in taxa[t]:
              taxaName = taxa[t][0:taxa[t].find('(')]
            else:
              taxaName = taxa[t].strip()
            taxaDict[self.ranksByLevel[t]] = taxaName
            taxaList.append(taxaName)

          taxaDict['taxonomy'] = taxaList
          taxonomy[taxonId] = taxaDict
        else:
          counts['Unclassified'] = counts.get('Unclassified', 0) + 1

  def write(self, fout, sampleCounts, taxonomy):
    sampleIds = sorted(sampleCounts.keys())
    otuIds = sorted(taxonomy.keys())

    otuMetadata = []
    sparseData = []
    for rowIndex, otuId in enumerate(otuIds):
      otuMetadata.append(taxonomy[otuId])

      for colIndex, sampleId in enumerate(sampleIds):
        if otuId in sampleCounts[sampleId]:
          sparseData.append([rowIndex, colIndex, sampleCounts[sampleId][otuId]])

    t = table_factory(sparseData, sampleIds, otuIds, None, otuMetadata, constructor=SparseOTUTable)

    t.getBiomFormatJsonString("CommunityM", direct_io=fout)

  def run(self, configFile, bootstrap, rank, bAbsoluteValues, fout):
    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    # read classification results for all sequence files in each sample
    sampleCounts = {}
    taxonomy = {}
    rankIndex = self.ranksByLabel[rank]
    for sample in sampleParams:
      prefix = projectParams['output_dir'] + sample

      counts = {}

      pairs = sampleParams[sample]['pairs']
      for i in xrange(0, len(pairs), 2):
        pair = pairs[i]
        classificationFile = prefix + '.' + pair[pair.rfind('/')+1:pair.rfind('.')] + '.classified.16S.tsv'
        self.parseClassificationFile(classificationFile, bootstrap, rankIndex, counts, taxonomy)

      singles = sampleParams[sample]['singles']
      for single in singles:
        classificationFile = prefix + '.' + single[single.rfind('/')+1:single.rfind('.')] + '.classified.16S.tsv'
        self.parseClassificationFile(classificationFile, bootstrap, rankIndex, counts, taxonomy)

      if not bAbsoluteValues:
        sumCounts = 0
        for taxa, count in counts.iteritems():
          sumCounts += count

        for taxa in counts:
          counts[taxa] /= float(sumCounts)

      sampleCounts[sampleParams[sample]['name']] = counts

    # write out results in BIOM format
    self.write(fout, sampleCounts, taxonomy)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Build table in BIOM format summarizing classified 16S sequences.")

  parser.add_argument('config_file', help='project config file')
  parser.add_argument('-o', '--output', help='output file (default = stdout).', type=argparse.FileType('w'), default=sys.stdout)
  parser.add_argument('-b', '--bootstrap', help='bootstrap threshold required to accept classification (default = 0)', type=int, default=0)
  parser.add_argument('-a', '--absolute', help='write absolute values instead of relative values', action='store_true')
  parser.add_argument('-r', '--rank', help='taxonomic rank of table (choices: Domain, Phylum, Class, Order, Family, Genus, Species, GG_ID), (default = GG_ID)',
                            choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'GG_ID'], default='GG_ID')

  parser.add_argument('--version', help='show version number of program', action='version', version='Build classification table v0.0.1')

  args = parser.parse_args()

  buildTable = BuildTable()
  buildTable.run(args.config_file, args.bootstrap, args.rank, args.absolute, args.output)
