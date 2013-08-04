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

from readConfig import ReadConfig

class LinkBinsTo16S(object):
  def __init__(self):
    pass

  def link16S(self, pairs, prefix, outputDir, contigFile, assemblies16S, binDir, threads):
    for i in xrange(0, len(pairs), 2):
      pair1 = pairs[i]
      pair2 = pairs[i+1]

      reads16S_1 = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')] + '.1.16S.fasta'
      reads16S_2 = prefix + '.' + pair2[pair2.rfind('/')+1:pair2.rfind('.')] + '.2.16S.fasta'

      # create combined file with reference sequences and assembled 16S sequences
      combinedFile = prefix + '.combined.fasta'
      os.system('cat ' + contigFile + ' ' + assemblies16S + ' > ' + combinedFile)

      # determine links between 16S reads to reference sequences and assembled 16S sequences
      print '\nMapping 16S reads to reference and assembled 16S sequences.'
      rtn = os.system('mapReads.py -t ' + str(threads) + ' ' + combinedFile + ' ' + reads16S_1 + ' ' + reads16S_2)

      print '\nCreating graph showing linked 16S reads.'
      bamFile = combinedFile[0:combinedFile.rfind('.')] + '.bam'
      graphFile = prefix + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')] + '.gdf'
      rtn = os.system('createPairedEndGraph.py ' + bamFile + ' ' + binDir + ' ' + graphFile)
      print '  Graph written to: ' + graphFile

      print '\nTabulating links.'
      linkFile = prefix + '.linkedUnbinnedContigs.tsv'
      rtn = os.system('associateUnbinnedSeqs.py ' + graphFile + ' ' + linkFile)
      print '  Writting linked 16S reads to: ' + linkFile

  def run(self, configFile, contigFile, assemblies16S, binDir, threads):
    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    for sample in sampleParams:
      outputDir = projectParams['output_dir']
      prefix = outputDir + sample
      pairs = sampleParams[sample]['pairs']

      # identify 16S sequences in paired-end reads
      self.link16S(pairs, prefix, outputDir, contigFile, assemblies16S, binDir, threads)

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
