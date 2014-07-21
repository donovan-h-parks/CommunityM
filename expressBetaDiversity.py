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
Calculate phylogenetic beta diversity between samples.
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
import argparse
import operator

class ExpressBetaDiversity(object):
    def __init__(self):
        self.treeFiles = {'GG94':'/srv/whitlam/bio/db/communitym/201308_gg/94_otus.tree',
                          'GG97':'/srv/whitlam/bio/db/communitym/201308_gg/97_otus.tree',
                          'GG99':'/srv/whitlam/bio/db/communitym/201308_gg/99_otus.tree' }
    
    def __covertTableToEBD(self, inputTable, outputTable):
        """Convert table to format compatible with EBD."""
        
        sampleOTUs = {}
        otuIds = set([])
        
        for line in open(inputTable):
            if '#OTU_ID' in line:
                lineSplit = line.split('\t')
                sampleIds = [x.strip() for x in lineSplit[1:-1]]
                
                for sampleId in sampleIds:
                    sampleOTUs[sampleId] = {}
                
                continue
            
            elif line[0] == '#' or line.strip() == '':
                continue
            else:
                lineSplit = line.split('\t')
                otuId = lineSplit[-1].split(';')[7].strip()
                otuId = otuId.replace('id__', '')
                if otuId == 'unmapped':
                    continue
                
                counts = [float(x) for x in lineSplit[1:len(sampleIds) + 1]]
                
                for i in xrange(0, len(sampleIds)):
                    sampleOTUs[sampleIds[i]][otuId] = sampleOTUs[sampleIds[i]].get(otuId, 0) + counts[i]
                otuIds.add(otuId)

        # write EBD OTU table
        fout = open(outputTable, 'w')
        sortedOTUs = sorted(list(otuIds))
        for otuId in sortedOTUs:
            fout.write('\t' + otuId)
        fout.write('\n')
        
        for sampleId in sampleOTUs:
            fout.write(sampleId)
            
            for otuId in sortedOTUs:
                if otuId in sampleOTUs[sampleId]:
                    fout.write('\t' + str(sampleOTUs[sampleId][otuId]))
                else:
                    fout.write('\t0')
            fout.write('\n')
      
        fout.close()
    
    def run(self, table, refDB, outputDir, perc):
        """Run Express Beta Diversity (EBD)."""
        
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        
        # convert table to format compatible with EBD
        ebdTable = os.path.join(outputDir, 'otu_table.ebd.tsv')
        self.__covertTableToEBD(table, ebdTable)
        
        # get number of sequences in each sample
        sampleCounts = os.path.join(outputDir, 'sample_counts.tsv')
        os.system('ebd -s %s -z > %s' % (ebdTable, sampleCounts))
        
        counts = {}
        with open(sampleCounts) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')
                if len(lineSplit) == 2:
                    counts[lineSplit[0]] = float(lineSplit[1])
                
        sortedCounts = sorted(counts.iteritems(), key=operator.itemgetter(1))
        jackknifeSeqs = int(sortedCounts[0][1] * perc)
        
        print ''
        print '  Sample %s contains the fewest sequences with %d' % (sortedCounts[0][0], sortedCounts[0][1])
        print '  Performing jackknife analysis with %d sequences.' % jackknifeSeqs
        
        # run EBD
        treeFile = self.treeFiles[refDB]
        
        print ''
        print '  Calculating weighted Soergel dissimilarity values.'
        prefix = os.path.join(outputDir, 'soergel')
        os.system('ebd -w -t %s -s %s -c %s -p %s -j %d -d %d' % (treeFile, ebdTable, 'Soergel', prefix, 100, jackknifeSeqs))
        
        print '  Calculating unweighted Soergel dissimilarity values.'
        prefix = os.path.join(outputDir, 'usoergel')
        os.system('ebd -t %s -s %s -c %s -p %s -j %d -d %d' % (treeFile, ebdTable, 'Soergel', prefix, 100, jackknifeSeqs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate phylogenetic beta diversity between samples",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('table', help='CommunityM table built with rank=SEQ_ID')
    parser.add_argument('ref_db', help='reference DB used to classify sequences (choices: GG94, GG97, GG99)', choices=['GG94', 'GG97', 'GG99'])
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('-p', '--perc', help='percentage of sequences to retain during jackknife analysis', type=float, default=0.9)

    args = parser.parse_args()

    expressBetaDiversity = ExpressBetaDiversity()
    expressBetaDiversity.run(args.table, args.ref_db, args.output_dir, args.perc)
