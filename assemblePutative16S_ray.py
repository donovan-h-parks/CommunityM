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
Assemble putative 16S genes.
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
import shutil

from readConfig import ReadConfig


class AssemblePutative16S(object):
    def __init__(self):
        pass

    def parseContigInfo(self, outputDir):
        contigLens = []
        for line in open(outputDir + '/ContigLengths.txt'):
            lineSplit = line.split('\t')
            contigLens.append(lineSplit[1].rstrip())

        return contigLens

    def run(self, configFile, threads, kmerLen, minContigLen):
        rc = ReadConfig()
        projectParams, _ = rc.readConfig(configFile, outputDirExists=True)

        # create directory to store putative 16S genes
        dirPutative16S = os.path.join(projectParams['output_dir'], 'putativeSSU/')
        if not os.path.exists(dirPutative16S):
            print '[Error] Putative 16S gene reads expected in: ' + dirPutative16S
            sys.exit()

        # extract GreenGene Ids of putative 16S genes
        ggIds = set()
        files = os.listdir(dirPutative16S)
        for f in files:
            if f.endswith('fasta'):
                ggIds.add(int(f.split('.')[0]))

        print 'Putative 16S genes to assemble: ' + str(len(ggIds))

        contigInfo = {}
        for ggId in ggIds:
            print 'Assembling ' + str(ggId) + ': '
            print ''

            pair1 = dirPutative16S + str(ggId) + '.1.fasta'
            pair2 = dirPutative16S + str(ggId) + '.2.fasta'
            single = dirPutative16S + str(ggId) + '.singletons.fasta'

            outputDir = dirPutative16S + str(ggId) + '_assembly'
            if os.path.exists(outputDir):
                shutil.rmtree(outputDir)

            cmd = 'mpiexec -n ' + str(threads) + ' Ray -k ' + str(kmerLen) + ' -minimum-contig-length ' + str(minContigLen) + ' -o ' + outputDir
            if os.stat(single).st_size > 0:  # check if file contains any sequences
                cmd += ' -s ' + single
            if os.stat(pair1).st_size > 0:
                cmd += ' -p ' + pair1 + ' ' + pair2

            os.system(cmd)

            contigInfo[ggId] = self.parseContigInfo(outputDir)

        print '\n*********************************'
        allContigsFile = projectParams['output_dir'] + 'assembled_contigs.16S.fasta'
        fout = open(allContigsFile, 'w')
        print 'Assembly results: '
        for ggId in contigInfo:
            print '  Assembly of ' + str(ggId) + ' produce ' + str(len(contigInfo[ggId])) + ' contig(s): ' + ' '.join(contigInfo[ggId])

            index = 0
            for line in open(dirPutative16S + str(ggId) + '_assembly/Contigs.fasta'):
                if line[0] == '>':
                    lineSplit = line.split()
                    seqLen = lineSplit[1]

                    fout.write('>16S_' + str(ggId) + '-' + str(index) + ' ' + seqLen + '\n')

                    index += 1
                else:
                    fout.write(line)
        fout.close()

        print ''
        print '  All assembled 16S contigs written to: ' + allContigsFile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Assemble putative 16S genes.",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config_file', help='project config file')

    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=1)
    parser.add_argument('-k', '--kmer_len', help='kmer length used for assembly', type=int, default=21)
    parser.add_argument('-m', '--min_contig_len', help='minimum contig length to retain', type=int, default=800)
    parser.add_argument('--version', help='show version number of program', action='version', version='Assemble 16S v0.0.1')

    args = parser.parse_args()

    assemblePutative16S = AssemblePutative16S()
    assemblePutative16S.run(args.config_file, args.threads, args.kmer_len, args.min_contig_len)
