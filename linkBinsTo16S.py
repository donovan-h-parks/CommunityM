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
Link assembled 16S sequences to bins. This file makes extensive use of the pairm project.
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
import ntpath

from readConfig import ReadConfig

class LinkBinsTo16S(object):
    def __init__(self):
        pass

    def link16S(self, combinedFile, ssuReads1, ssuReads2, binDir, threads, outputDir):
        # determine links between 16S reads to reference sequences and assembled 16S sequences
        print 'Mapping 16S reads to reference and assembled 16S sequences.'
        os.system('mapReads.py -t ' + str(threads) + ' ' + combinedFile + ' ' + ssuReads1 + ' ' + ssuReads2)

        print '\nCreating graph showing linked 16S reads.'
        bamFile = combinedFile[0:combinedFile.rfind('.')] + '.bam'
        graphFile = outputDir + 'linksToBins.gdf'
        os.system('createPairedEndGraph.py ' + bamFile + ' ' + binDir + ' ' + graphFile)
        print '  Graph written to: ' + graphFile

        print 'Tabulating links.'
        linkFile = outputDir + 'linkedContigs.tsv'
        os.system('associateUnbinnedSeqs.py ' + graphFile + ' ' + linkFile)
        print '  Link file written reads to: ' + linkFile

    def run(self, configFile, contigFile, assemblies16S, binDir, threads):
        rc = ReadConfig()
        projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

        # check if links directory already exists
        if not os.path.exists(projectParams['output_dir'] + 'linksToBin'):
            os.makedirs(projectParams['output_dir'] + 'linksToBin')
        else:
            rtn = raw_input('Remove previously identified links (Y or N)? ')
            if rtn.lower() == 'y' or rtn.lower() == 'yes':
                files = os.listdir(projectParams['output_dir'] + 'linksToBin')
                for f in files:
                    os.remove(projectParams['output_dir'] + 'linksToBin/' + f)
            else:
                sys.exit()

        outputDir = projectParams['output_dir'] + 'linksToBin/'

        # create combined file with reference sequences and assembled 16S sequences
        print 'Combining unbinned reference sequences with de novo assembled 16S sequences.'
        combinedFile = outputDir + 'scaffolds.combined.fasta'
        os.system('cat ' + contigFile + ' ' + assemblies16S + ' > ' + combinedFile)

        # create combined 16S read files
        print 'Combining 16S/18S reads from all samples.'
        reads1 = ''
        reads2 = ''
        for sample in sampleParams:
            extractedPrefix = projectParams['output_dir'] + 'extracted/' + sample
            pairs = sampleParams[sample]['pairs']
            for i in xrange(0, len(pairs), 2):
                pair1Base = ntpath.basename(pairs[i])
                pair2Base = ntpath.basename(pairs[i+1])

                classificationFile1 = extractedPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.union.SSU.fasta'
                classificationFile2 = extractedPrefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.union.SSU.fasta'

                reads1 += classificationFile1 + ' '
                reads2 += classificationFile2 + ' '

        os.system('cat ' + reads1 + ' > ' + outputDir + 'ssu.1.fasta')
        os.system('cat ' + reads2 + ' > ' + outputDir + 'ssu.2.fasta')

        # identify 16S sequences in paired-end reads
        self.link16S(combinedFile, outputDir + 'ssu.1.fasta', outputDir + 'ssu.2.fasta', binDir, threads, outputDir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Link assembled 16S sequences to bins.",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config_file', help='project config file')
    parser.add_argument('contig_file', help='FASTA file of assembled contigs/scaffolds')
    parser.add_argument('assemblies16S', help='FASTA file of de novo 16S assemblies')
    parser.add_argument('bin_dir', help='directory containing bin')

    parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)

    args = parser.parse_args()

    linkBinsTo16S = LinkBinsTo16S()
    linkBinsTo16S.run(args.config_file, args.contig_file, args.assemblies16S, args.bin_dir, args.threads)
