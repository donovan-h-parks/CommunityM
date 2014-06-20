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

import argparse
import ntpath
import os

from readConfig import ReadConfig
from taxonomyUtils import ranksByLabel, ranksByLevel, rankPrefixes, LCA, parseTaxon

class BuildTable(object):
    def __init__(self):
        pass

    def parseClassification(self, seqClassification, bIgnoreUnmapped, bootstrapThreshold, rankIndex, counts, taxonomy):
        countUnmapped = 0
        for seqId in seqClassification:
            taxa = seqClassification[seqId]

            taxaName, _ = parseTaxon(taxa[rankIndex])

            if bIgnoreUnmapped:
                # check if taxa is unmapped, or if we are building a table at the sequence id level check if the species is unmapped
                if 'unmapped' in taxaName or (rankIndex == ranksByLabel['SEQ_ID'] and 'unmapped' in parseTaxon(taxa[ranksByLabel['Species']])[0]):
                    countUnmapped += 1
                    continue

            taxaDict = {}
            taxaList = []
            for r in xrange(0, rankIndex+1):
                taxaName, bootstrapSupport = parseTaxon(taxa[r])

                if bootstrapSupport < bootstrapThreshold:
                    taxaName = rankPrefixes[r] + 'unclassified'

                taxaDict[ranksByLevel[r]] = taxaName
                taxaList.append(taxaName)

            taxaDict['taxonomy'] = taxaList

            taxaStr = ';'.join(taxaList)
            taxonomy[taxaStr] = taxaDict
            counts[taxaStr] = counts.get(taxaStr, 0) + 1

        if bIgnoreUnmapped:
            print '  [Note] Ignoring ' + str(countUnmapped) + ' unmapped reads.'
        print ''

    def parsePairedClassificationFiles(self, classificationFile1, classificationFile2, bIgnoreUnmapped, bTreatPairsAsSingles, bootstrapThreshold, rankIndex, counts, taxonomy):
        print 'Processing: '
        print '  ' + classificationFile1
        print '  ' + classificationFile2

        seqClassification = {}
        with open(classificationFile1) as fin:
            for line in fin:
                lineSplit = line.split('\t')
                seqId = lineSplit[0]
                taxa = lineSplit[1].split(';')

                if not bTreatPairsAsSingles:
                    seqId = seqId[0:seqId.rfind('/')]

                seqClassification[seqId] = taxa

        identicalSeq = 0
        nonIdenticalSeq = 0
        with open(classificationFile2) as fin:
            for line in fin:
                lineSplit = line.split('\t')
                seqId = lineSplit[0]
                taxa = lineSplit[1].split(';')

                if not bTreatPairsAsSingles:
                    seqId = seqId[0:seqId.rfind('/')]
                    
                    
                    
                    if 'unmapped' not in ';'.join(seqClassification[seqId]) and 'unmapped' not in ';'.join(taxa):
                        if seqClassification[seqId][7] == taxa[7]:
                            identicalSeq += 1
                        else:
                            nonIdenticalSeq += 1
                    
                    
                    seqClassification[seqId] = LCA(taxa, seqClassification[seqId])                           
                else:
                    seqClassification[seqId] = taxa

        self.parseClassification(seqClassification, bIgnoreUnmapped, bootstrapThreshold, rankIndex, counts, taxonomy)

    def parseSingleClassificationFile(self, classificationFile, bIgnoreUnmapped, bootstrapThreshold, rankIndex, counts, taxonomy):
        print 'Processing: '
        print '  ' + classificationFile

        seqClassification = {}
        with open(classificationFile) as fin:
            for line in fin:
                lineSplit = line.split('\t')
                seqId = lineSplit[0]
                taxa = lineSplit[1].split(';')
                seqClassification[seqId] = taxa

        self.parseClassification(seqClassification, bIgnoreUnmapped, bootstrapThreshold, rankIndex, counts, taxonomy)

    def writeTable(self, output, sampleCounts, taxonomy):
        fout = open(output, 'w')
        
        fout.write("#OTU_ID")
        fout.write("\t" + "\t".join(sampleCounts.keys()))
        fout.write("\tConsensus Lineage\n")

        for outId, taxonomy in enumerate(taxonomy.keys()):
            fout.write(str(outId))
            
            for sampleId in sampleCounts.keys():
                if taxonomy in sampleCounts[sampleId]:
                    fout.write('\t%f' % sampleCounts[sampleId][taxonomy])
                else:
                    fout.write('\t0')
                    
            fout.write('\t' + taxonomy + '\n')
                      
        fout.close()

    def paramsToCountsAndTaxonomy(self, projectParams, sampleParams, bIgnoreUnmapped, bTreatPairsAsSingles, bUseSingletons, bootstrap, rank, bMode):
        # read classification results for all sequence files in each sample
        sampleCounts = {}
        taxonomy = {}
        rankIndex = ranksByLabel[rank]
        for sample in sampleParams:
            prefix = os.path.join(projectParams['output_dir'],'classified/',sample)

            counts = {}

            pairs = sampleParams[sample]['pairs']
            for i in xrange(0, len(pairs), 2):
                pair1Base = ntpath.basename(pairs[i])
                pair2Base = ntpath.basename(pairs[i+1])

                classificationFile1 = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.16S.tsv'
                classificationFile2 = prefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.16S.tsv'
                self.parsePairedClassificationFiles(classificationFile1, classificationFile2, bIgnoreUnmapped, bTreatPairsAsSingles, bootstrap, rankIndex, counts, taxonomy)

                if bUseSingletons:
                    classificationFile = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.16S.tsv'
                    self.parseSingleClassificationFile(classificationFile, bIgnoreUnmapped, bootstrap, rankIndex, counts, taxonomy)

                    singles = sampleParams[sample]['singles']
                    for single in singles:
                        singleBase = ntpath.basename(single)
                        classificationFile = prefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S.tsv'
                        self.parseSingleClassificationFile(classificationFile, bIgnoreUnmapped, bootstrap, rankIndex, counts, taxonomy)

            if bMode == "rel":
                # relative values
                sumCounts = 0
                for taxa, count in counts.iteritems():
                    sumCounts += count

                for taxa in counts:
                    counts[taxa] /= float(sumCounts)
            elif bMode == "pre":
                # presence absence
                for taxa in counts:
                    if counts[taxa] > 0:
                        counts[taxa] = 1.
                    else:
                        counts[taxa] = 0.

            sampleCounts[sampleParams[sample]['name']] = counts

        return sampleCounts, taxonomy


    def run(self, projectParams, sampleParams, bIgnoreUnmapped, bTreatPairsAsSingles, bUseSingletons, bootstrap, rank, bMode, output):
        sampleCounts, taxonomy = self.paramsToCountsAndTaxonomy(projectParams, sampleParams, bIgnoreUnmapped, bTreatPairsAsSingles, bUseSingletons, bootstrap, rank, bMode)
        self.writeTable(output, sampleCounts, taxonomy)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Build table summarizing classified 16S sequences.")

    parser.add_argument('config_file', help='project config file')
    parser.add_argument('output_table', help='output file')
    parser.add_argument('-u', '--ignore_unmapped', help='do not consider unmapped reads', action="store_true")
    parser.add_argument('-p', '--pairs_as_singles', help='treat paired reads as singletons', action="store_true")
    parser.add_argument('-s', '--singletons', help='use singleton 16S/18S reads', action="store_true")
    parser.add_argument('-b', '--bootstrap', help='bootstrap threshold required to accept classification (default = 0)', type=int, default=0)
    parser.add_argument('-m', '--mode', help='write values as "rel"ative, "abs"olute or "pre"sence/absense (default = abs)', choices=['rel', 'abs', 'pre'], default="abs")
    parser.add_argument('-r', '--rank', help='taxonomic rank of table (choices: Domain, Phylum, Class, Order, Family, Genus, Species, SEQ_ID), (default = SEQ_ID)',
                              choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SEQ_ID'], default='SEQ_ID')

    args = parser.parse_args()

    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(args.config_file, outputDirExists = True)

    buildTable = BuildTable()
    buildTable.run(projectParams, sampleParams, args.ignore_unmapped, args.pairs_as_singles, args.singletons, args.bootstrap, args.rank, args.mode, args.output_table)
