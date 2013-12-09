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
Classify 16S fragments by mapping them to the GreenGenes DB with BWA.
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
from bwaUtils import mapPair, mapSingle
from taxonomyUtils import LCA, readTaxonomy

import pysam

class ClassifyBWA(object):
    def __init__(self):
        self.unmappedStr = 'k__unmapped;p__unmapped;c__unmapped;o__unmapped;f__unmapped;g__unmapped;s__unmapped;id__unmapped;'
        self.ggDB = '/srv/db/gg/2013_05/gg_13_5_otus/rep_set/##_otus.fasta'
        self.taxonomyFile = '/srv/db/gg/2013_05/gg_13_5_otus/taxonomy/##_otu_taxonomy.full.txt'

    def processRead(self, bam, read, ggIdToTaxonomy, maxEditDistance, minLength, counts = None):
        bMapped = False
        if read.is_unmapped:
            if counts != None:
                counts['unmapped'] += 1
            return self.unmappedStr, bMapped
        elif (read.alen < minLength*read.rlen):
            if counts != None:
                counts['align len'] += 1
            return self.unmappedStr, bMapped
        elif (read.opt('NM') > maxEditDistance*read.rlen):
            if counts != None:
                counts['edit dist'] += 1
            return self.unmappedStr, bMapped
        else:
            taxonomy = ggIdToTaxonomy[bam.getrname(read.tid)]
            bMapped = True

        return taxonomy, bMapped

    def readPairedBAM(self, bamFile, ggIdToTaxonomy, maxEditDistance, minLength):
        # read compressed BAM file and report basic statistics
        bam = pysam.Samfile(bamFile, 'rb')

        # find primary mappings for each query read
        readsMappedTo16S_1 = {}
        readsMappedTo16S_2 = {}
        editDists = {}
        counts = {'unmapped':0, 'edit dist':0, 'align len':0}

        for read in bam.fetch(until_eof=True):
            if not read.is_secondary:
                taxonomy, bMapped = self.processRead(bam, read, ggIdToTaxonomy, maxEditDistance, minLength, counts)

                if bMapped:
                    editDist = read.opt('NM')
                else:
                    editDist = -1 # flag as unmapped

                if read.is_read1:
                    qname = read.qname + '/1'
                    readsMappedTo16S = readsMappedTo16S_1
                elif read.is_read2:
                    qname = read.qname + '/2'
                    readsMappedTo16S = readsMappedTo16S_2

                if qname in readsMappedTo16S and editDists[qname] != -1:
                    # read has multiple primary alignments for different parts of the query sequence
                    # which may indicate it is chimeric. For classification purposes, the LCA of
                    # all primary alignments is taken.
                    lca = LCA(readsMappedTo16S[qname].split(';'), taxonomy.split(';'))
                    readsMappedTo16S[qname] = ';'.join(lca)
                    editDists[qname] = max(editDist, editDists[qname])
                else:
                    readsMappedTo16S[qname] = taxonomy
                    editDists[qname] = editDist

        # process secondary mappings for each query read
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                # process primary read
                taxonomy, bMapped = self.processRead(bam, read, ggIdToTaxonomy, maxEditDistance, minLength)
                editDist = read.opt('NM')

                if read.is_read1:
                    qname = read.qname + '/1'
                    readsMappedTo16S = readsMappedTo16S_1
                elif read.is_read2:
                    qname = read.qname + '/2'
                    readsMappedTo16S = readsMappedTo16S_2

                if bMapped and editDist <= editDists[qname]:
                    lca = LCA(readsMappedTo16S[qname].split(';'), taxonomy.split(';'))
                    readsMappedTo16S[qname] = ';'.join(lca)

        bam.close()

        if len(readsMappedTo16S_1) != len(readsMappedTo16S_2):
            print '[Error] Paired files do not have the same number of reads.'
            sys.exit()

        print '  Number of reads: ' + str(len(readsMappedTo16S))
        print '    Reads unmapped: ' + str(counts['unmapped'])
        print '    Reads failing edit distance threshold: ' + str(counts['edit dist'])
        print '    Reads failing alignment length threshold: ' + str(counts['align len'])

        return readsMappedTo16S_1, readsMappedTo16S_2

    def readSingleBAM(self, bamFile, ggIdToTaxonomy, maxEditDistance, minLength):
        # read compressed BAM file
        bam = pysam.Samfile(bamFile, 'rb')

        # find primary mappings for each query read
        readsMappedTo16S = {}
        editDists = {}
        counts = {'unmapped':0, 'edit dist':0, 'align len':0}

        for read in bam.fetch(until_eof=True):
            if not read.is_secondary:
                taxonomy, bMapped = self.processRead(bam, read, ggIdToTaxonomy, maxEditDistance, minLength, counts)

                if bMapped:
                    editDist = read.opt('NM')
                else:
                    editDist = -1 # flag as unmapped

                if read.qname in readsMappedTo16S and editDists[read.qname] != -1:
                    # read has multiple primary alignments for different parts of the query sequence
                    # which may indicate it is chimeric. For classification purposes, the LCA of
                    # all primary alignments is taken.
                    lca = LCA(readsMappedTo16S[read.qname].split(';'), taxonomy.split(';'))
                    readsMappedTo16S[read.qname] = ';'.join(lca)
                    editDists[read.qname] = max(editDist, editDists[read.qname])
                else:
                    readsMappedTo16S[read.qname] = taxonomy
                    editDists[read.qname] = editDist

        # process secondary mappings for each query read
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                # process primary read
                taxonomy, bMapped = self.processRead(bam, read, ggIdToTaxonomy, maxEditDistance, minLength)
                editDist = read.opt('NM')

                if bMapped and editDist <= editDists[read.qname]:
                    lca = LCA(readsMappedTo16S[read.qname].split(';'), taxonomy.split(';'))
                    readsMappedTo16S[read.qname] = ';'.join(lca)

        bam.close()

        print '  Number of reads: ' + str(len(readsMappedTo16S))
        print '    Reads unmapped: ' + str(counts['unmapped'])
        print '    Reads failing edit distance threshold: ' + str(counts['edit dist'])
        print '    Reads failing alignment length threshold: ' + str(counts['align len'])

        return readsMappedTo16S

    def writeClassification(self, filename, mappedReads):
        fout = open(filename, 'w')
        for refName, taxonomy in mappedReads.iteritems():
            fout.write(refName + '\t' + taxonomy + '\n')
        fout.close()

    def processPairs(self, pairs, ggIdToTaxonomy, maxEditDistance, minLength, outputDir, prefix):
        for i in xrange(0, len(pairs), 2):
            pair1 = pairs[i]
            pair2 = pairs[i+1]

            pair1Base = ntpath.basename(pair1)
            pair2Base = ntpath.basename(pair2)

            print 'Identifying 16S sequences in paired-end reads: ' + pair1 + ', ' + pair2

            # write out classifications for paired-end reads with both ends identified as 16S
            bamFile = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.16S.bam'
            readsMappedTo16S_1, readsMappedTo16S_2  = self.readPairedBAM(bamFile, ggIdToTaxonomy, maxEditDistance, minLength)

            output1 = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.16S.tsv'
            output2 = prefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.16S.tsv'
            print '  Paired results written to: '
            print '    ' + output1
            print '    ' + output2 + '\n'
            self.writeClassification(output1, readsMappedTo16S_1)
            self.writeClassification(output2, readsMappedTo16S_2)

            # write out classifications for paired-ends reads with only one end identified as 16S
            bamFile = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.16S.bam'
            readsMappedTo16S = self.readSingleBAM(bamFile, ggIdToTaxonomy, maxEditDistance, minLength)

            output = prefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.16S.tsv'
            print '  Singleton results written to: ' + output + '\n'
            self.writeClassification(output, readsMappedTo16S)

    def processSingles(self, singles, ggIdToTaxonomy, maxEditDistance, minLength, outputDir, prefix):
        for i in xrange(0, len(singles)):
            seqFile = singles[i]

            print 'Identifying 16S sequences in single-end reads: ' + seqFile

            singleBase = ntpath.basename(seqFile)
            bamFile = prefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S.bam'
            readsMappedTo16S = self.readSingleBAM(bamFile, ggIdToTaxonomy, maxEditDistance, minLength)

            output = prefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S.tsv'
            print '  Classification results written to: ' + output + '\n'
            self.writeClassification(output, readsMappedTo16S)

    def run(self, projectParams, sampleParams, otu, threads):
        # check if classification directory already exists
        dir_path = os.path.join(projectParams['output_dir'],'classified')
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        else:
            rtn = raw_input('Remove previously classified reads (Y or N)? ')
            if rtn.lower() == 'y' or rtn.lower() == 'yes':
                files = os.listdir(dir_path)
                for f in files:
                    os.remove(os.path.join(dir_path, f))
            else:
                sys.exit()

        taxonomyFile = self.taxonomyFile.replace('##', str(otu), 1)
        ggIdToTaxonomy = readTaxonomy(taxonomyFile)

        ggDB = self.ggDB.replace('##', str(otu), 1)

        print 'Classifying reads with: ' + ggDB
        print 'Assigning taxonomy with: ' + taxonomyFile
        print 'Threads: ' + str(threads)
        print ''

        if not os.path.exists(ggDB + '.amb'):
            print 'Indexing GreenGenes DB:'
            os.system('bwa index -a is ' + ggDB)
            print ''

        # map reads
        for sample in sampleParams:
            print 'Mapping sample: ' + sample
            outputDir = projectParams['output_dir']
            inputPrefix = os.path.join(outputDir,'extracted',sample)
            outputPrefix = os.path.join(outputDir,'classified',sample)
            pairs = sampleParams[sample]['pairs']
            singles = sampleParams[sample]['singles']

            for i in xrange(0, len(pairs), 2):
                pair1Base = ntpath.basename(pairs[i])
                pair1File = inputPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.intersect.SSU.fasta'

                pair2Base = ntpath.basename(pairs[i+1])
                pair2File = inputPrefix + '.' + pair2Base[0:pair2Base.rfind('.')] + '.intersect.SSU.fasta'

                bamPrefix = ntpath.basename(pairs[i])
                bamPrefixFile = outputPrefix + '.' + bamPrefix[0:bamPrefix.rfind('.')] + '.intersect.16S'
                mapPair(ggDB, pair1File, pair2File, bamPrefixFile, threads)

                diffFile = inputPrefix + '.' + pair1Base[0:pair1Base.rfind('.')] + '.difference.SSU.fasta'
                bamPrefixFile = outputPrefix + '.' + bamPrefix[0:bamPrefix.rfind('.')] + '.difference.16S'
                mapSingle(ggDB, diffFile, bamPrefixFile, threads)

            for i in xrange(0, len(singles)):
                singleBase = ntpath.basename(singles[i])
                singleFile = inputPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.SSU.fasta'

                bamPrefixFile = outputPrefix + '.' + singleBase[0:singleBase.rfind('.')] + '.16S'
                mapSingle(ggDB, singleFile, bamPrefixFile, threads)

            print '************************************************************'

        # classify reads
        for sample in sampleParams:
            print 'Classifying sample: ' + sample
            outputDir = os.path.join(projectParams['output_dir'],'classified')
            prefix = os.path.join(outputDir,sample)
            pairs = sampleParams[sample]['pairs']
            singles = sampleParams[sample]['singles']
            maxEditDistance = sampleParams[sample]['edit_dist']
            minLength = sampleParams[sample]['min_align_len']

            # identify 16S sequences in paired-end reads
            self.processPairs(pairs, ggIdToTaxonomy, maxEditDistance, minLength, outputDir, prefix)

            # identify 16S sequences in single-end reads
            self.processSingles(singles, ggIdToTaxonomy, maxEditDistance, minLength, outputDir, prefix)

            print ''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Classify 16S fragments by mapping them to the GreenGenes DB with BWA.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config_file', help='project config file.')
    parser.add_argument('otu', help='GreenGenes DB to use for classification (choices: 94, 97, 99)', type=int, choices=[94, 97, 99], default=97)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)

    args = parser.parse_args()

    classifyBWA = ClassifyBWA()

    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)

    classifyBWA.run(projectParams, sampleParams, args.otu, args.threads)
