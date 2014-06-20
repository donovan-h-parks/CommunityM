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
Summarize OTU table at a specific taxonomic rank.
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

from taxonomyUtils import ranksByLabel

class SummarizeTableAtRank(object):
    def __init__(self):
        pass

    def run(self, inputTable, summaryRank, otherThreshold, outputFile):
        rankIndex = ranksByLabel[summaryRank]

        # summarize results at specified rank
        fout = open(outputFile, 'w')
        sampleIds = []
        table = {}
        for line in open(inputTable):
            lineSplit = line.split('\t')

            if line[0] == '#':
                fout.write(line)

                if lineSplit[-1].rstrip() == 'Consensus Lineage' or lineSplit[-1].rstrip() == 'taxonomy':
                    taxonomyIndex = len(lineSplit)-1
                    sampleIds = lineSplit[1:taxonomyIndex]
            else:
                taxonomy = [x.strip() for x in lineSplit[taxonomyIndex].split(';')]
                summaryLineage = ';'.join(taxonomy[0:rankIndex+1])

                if summaryLineage not in table:
                    table[summaryLineage] = [0 for _ in xrange(0, len(sampleIds))]

                for i in xrange(0, len(sampleIds)):
                    table[summaryLineage][i] += float(lineSplit[i+1])

        # assign low abundant taxa to 'other'
        other = [0 for _ in xrange(0, len(sampleIds))]
        for summaryLineage, row in table.iteritems():
            bOther = True
            for r in row:
                if r >= otherThreshold:
                    bOther = False
                    break

            if bOther and 'unmapped' not in summaryLineage:
                for i in xrange(0, len(row)):
                    other[i] += row[i]
            else:
                fout.write(summaryLineage)
                for r in row:
                    fout.write('\t' + str(r))
                fout.write('\n')

        fout.write('Other')
        for r in other:
            fout.write('\t' + str(r))
        fout.write('\n')

        fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize OTU table at a specific taxonomic rank.")

    parser.add_argument('table', help='QIIME-format OTU table')
    parser.add_argument('output', help='output file')
    parser.add_argument('-r', '--rank', help='taxonomic rank to summarize table at (choices: Domain, Phylum, Class, Order, Family, Genus, Species), (default = Genus)',
                                        choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'GG_ID'], default='GG_ID')
    parser.add_argument('-o', '--other', help="assign OTUs below this percent to a category called 'other' (default = 0)", type=float, default=0.0)

    args = parser.parse_args()

    summarizeTableAtRank = SummarizeTableAtRank()
    summarizeTableAtRank.run(args.table, args.rank, args.other, args.output)
