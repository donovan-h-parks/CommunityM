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
Assemble putative 16S genes with the SPADES assembler.
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

  def parseScaffoldInfo(self, outputDir):
    scaffoldLens = []

    if not os.path.isfile(outputDir + '/scaffolds.fasta'):
      return []

    for line in open(outputDir + '/scaffolds.fasta'):
      if line[0] == '>':
        lineSplit = line.split('_')
        seqLen = lineSplit[3]

        scaffoldLens.append(seqLen)

    return scaffoldLens

  def run(self, configFile, threads):
    rc = ReadConfig()
    projectParams, _ = rc.readConfig(configFile, outputDirExists = True)

    # create directory to store putative 16S genes
    dirPutative16S = projectParams['output_dir'] + 'putativeSSU/'
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

    scaffoldInfo = {}
    for ggId in ggIds:
      print 'Assembling ' + str(ggId) + ': '
      print ''

      pair1 = dirPutative16S + str(ggId) + '.1.fasta'
      pair2 = dirPutative16S + str(ggId) + '.2.fasta'
      single = dirPutative16S + str(ggId) + '.singletons.fasta'

      outputDir = dirPutative16S + str(ggId) + '_assembly_spades'
      if os.path.exists(outputDir):
        shutil.rmtree(outputDir)

      cmd = 'spades.py --only-assembler -o ' + outputDir + ' -t ' + str(threads)
      if os.stat(single).st_size > 0: # check if file contains any sequences
        cmd += ' -s ' + single
      if os.stat(pair1).st_size > 0:
        cmd += ' -1 ' + pair1 + ' -2 ' + pair2

      os.system(cmd)

      scaffoldInfo[ggId] = self.parseScaffoldInfo(outputDir)

    print '\n*********************************'
    allScaffoldsFile = projectParams['output_dir'] + 'assembled_scaffolds.16S.fasta'
    fout = open(allScaffoldsFile, 'w')
    print 'Assembly results: '
    for ggId in scaffoldInfo:
      print '  Assembly of ' + str(ggId) + ' produce ' + str(len(scaffoldInfo[ggId])) + ' scaffold(s): ' + ' '.join(scaffoldInfo[ggId])

      if not os.path.isfile(dirPutative16S + str(ggId) + '_assembly_spades/scaffolds.fasta'):
        print '    Failed to build scaffolds for ' + str(ggId)
        continue

      index = 0
      for line in open(dirPutative16S + str(ggId) + '_assembly_spades/scaffolds.fasta'):
        if line[0] == '>':
          lineSplit = line.split('_')
          seqLen = lineSplit[3]

          fout.write('>16S_' + str(ggId) + '-' + str(index) + ' ' + seqLen + '\n')

          index += 1
        else:
          fout.write(line)
    fout.close()

    print ''
    print '  All assembled 16S contigs written to: ' + allScaffoldsFile

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Assemble putative 16S genes with the SPADES assembler.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('config_file', help='project config file')

  parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)

  args = parser.parse_args()

  assemblePutative16S = AssemblePutative16S()
  assemblePutative16S.run(args.config_file, args.threads)
