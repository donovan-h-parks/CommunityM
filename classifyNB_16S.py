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
Classify 16S fragments using mothur's implementation of Wang's naive Bayes classifier.
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
import tempfile
import argparse
import ntpath

from readConfig import ReadConfig

class Classify16S(object):
    def __init__(self):
        self.dbFiles = {'GG94':'/srv/db/gg/2013_05/gg_13_5_otus/rep_set_aligned/94_otus.fasta',
                          'GG97':'/srv/db/gg/2013_05/gg_13_5_otus/rep_set_aligned/97_otus.fasta',
                          'GG99':'/srv/db/gg/2013_05/gg_13_5_otus/rep_set_aligned//99_otus.fasta',
                          'SILVA98':'/srv/whitlam/bio/db/mothur/silva/SSURef_111_NR_trunc.fna' }

        self.taxonomyFiles = {'GG94':'/srv/db/gg/2013_05/gg_13_5_otus/taxonomy/94_otu_taxonomy.full.txt',
                              'GG97':'/srv/db/gg/2013_05/gg_13_5_otus/taxonomy/97_otu_taxonomy.full.txt',
                              'GG99':'/srv/db/gg/2013_05/gg_13_5_otus/taxonomy/99_otu_taxonomy.full.txt',
                              'SILVA98':'/srv/whitlam/bio/db/mothur/silva/SSURef_111_NR_taxonomy.txt' }

    def classify(self, seqFile, dbFile, taxonomyFile, threads, bQuiet):
        tempFD, tempFilePath = tempfile.mkstemp(dir='.')
        fout = os.fdopen(tempFD,'w')
        fout.write("classify.seqs(fasta=" + seqFile + ", template=" + dbFile + ",taxonomy=" + taxonomyFile + ",processors=" + str(threads) + ")")
        fout.close()

        # classify seqs
        if bQuiet:
            os.system('mothur ' + tempFilePath + ' > /dev/null')
        else:
            os.system('mothur ' + tempFilePath)

        os.remove(tempFilePath)

    def run(self, configFile, db, threads, bQuiet):
        rc = ReadConfig()
        projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

        # check if classification directory already exists
        if not os.path.exists(projectParams['output_dir'] + 'classified'):
            os.makedirs(projectParams['output_dir'] + 'classified')
        else:
            rtn = raw_input('Remove previously classified reads (Y or N)? ')
            if rtn.lower() == 'y' or rtn.lower() == 'yes':
                files = os.listdir(projectParams['output_dir'] + 'classified')
                for f in files:
                    os.remove(projectParams['output_dir'] + 'classified/' + f)
            else:
                sys.exit()

        dbFile = self.dbFiles[db]
        taxonomyFile = self.taxonomyFiles[db]

        if not bQuiet:
            print 'Classifying reads with: ' + dbFile
            print 'Assigning taxonomy with: ' + taxonomyFile
            print 'Threads: ' + str(threads)
            print ''

        # create list of all sequence to classify
        mothurSeqFileList = ''
        for sample in sampleParams:
            prefix = projectParams['output_dir'] + 'extracted/' + sample
            pairs = sampleParams[sample]['pairs']
            singles = sampleParams[sample]['singles']

            for i in xrange(0, len(pairs), 2):
                pair1Base = ntpath.basename(pairs[i])
                pair1File = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.SSU.fasta'

                pair2Base = ntpath.basename(pairs[i+1])
                pair2File = prefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.SSU.fasta'

                diffFile = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.SSU.fasta'

                mothurSeqFileList += pair1File + '-' + pair2File + '-' + diffFile + '-'

            for single in singles:
                singleBase = ntpath.basename(single)
                singleFile = prefix + '.' + singleBase[0:singleBase.rfind('.')] + '.SSU.fasta'

                mothurSeqFileList += singleFile + '-'

        # classify with mothur
        mothurSeqFileList = mothurSeqFileList[0:-1] # remove trailing dash
        self.classify(mothurSeqFileList, dbFile, taxonomyFile, threads, bQuiet)

        # rename classification file for consistency with down-stream processing
        print 'Final classifications written to: '
        for filename in mothurSeqFileList.split('-'):
            if 'GG' in db:
                inputName = filename[0:filename.rfind('.')] + '.full.wang.taxonomy'
            else:
                inputName = filename[0:filename.rfind('.')] + '.SSURef_111_NR_taxonomy.wang.taxonomy'
            outputName = inputName.replace('/extracted/','/classified/')
            outputName = outputName.replace('SSU.full.wang.taxonomy','16S.tsv')
            os.system('mv ' + inputName + ' ' + outputName)
            print '  ' + outputName

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Classify 16S fragments using mothur's naive Bayes classifier.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config_file', help='project config file')
    parser.add_argument('--db', help='database to use for classification (choices: GG94, GG97, GG99, SILVA98)', choices=['GG94', 'GG97', 'GG99', 'SILVA98'], default='GG97')
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=1)
    parser.add_argument('-q', '--quiet', help='supress all output', action='store_true')

    parser.add_argument('--version', help='Show version number of program', action='version', version='Classify 16S using NB v0.0.1')

    args = parser.parse_args()

    classify16S = Classify16S()
    classify16S.run(args.config_file, args.db, args.threads, args.quiet)
