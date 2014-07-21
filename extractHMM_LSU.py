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
Extract LSU sequences from metagenomic data using HMMs.
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

from readConfig import ReadConfig
from seqUtils import extractSeqs

class ExtractLSU(object):
    def __init__(self):
        self.bacteriaModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_bacteria.hmm'
        self.bacteriaRevCompModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_bacteria.revComp.hmm'

        self.archaeaModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_archaea.hmm'
        self.archaeaRevCompModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_archaea.revComp.hmm'

        self.eukModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_euk.hmm'
        self.eukRevCompModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/LSU_euk.revComp.hmm'

    def getHits(self, hitTable):
        seqIds = set()
        for line in open(hitTable):
            if line[0] == '#' or line.strip() == '':
                continue

            seqId = line.split()[0].split('/')[0]   # remove read identifier
            seqIds.add(seqId)

        return seqIds

    def hmmSearch(self, seqFile, threads, evalue, outputPrefix):
        if not self.bQuiet:
            print '    Identifying bacterial LSU sequences'

        if seqFile.endswith('gz'):
            pipe = 'zcat ' + seqFile + ' | '
        else:
            pipe = 'cat ' + seqFile + ' | '

        if '.fq' in seqFile or '.fastq' in seqFile:
            pipe += 'fastq2fasta - | '

        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.bacteria.txt --tblout ' + outputPrefix + '.lsu.bacteria.table.txt -E ' + str(evalue) + ' ' + self.bacteriaModelFile + ' -')
        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.bacteria.rev_comp.txt --tblout ' + outputPrefix + '.lsu.bacteria.table.rev_comp.txt -E ' + str(evalue) + ' ' + self.bacteriaRevCompModelFile + ' -')

        if not self.bQuiet:
            print '    Identifying archaeal LSU sequences'
        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.archaea.txt --tblout ' + outputPrefix + '.lsu.archaea.table.txt -E ' + str(evalue) + ' ' + self.archaeaModelFile + ' -')
        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.archaea.rev_comp.txt --tblout ' + outputPrefix + '.lsu.archaea.table.rev_comp.txt -E ' + str(evalue) + ' ' + self.archaeaRevCompModelFile + ' -')

        if not self.bQuiet:
            print '    Identifying eukaryotic LSU sequences'
        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.euk.txt --tblout ' + outputPrefix + '.lsu.euk.table.txt -E ' + str(evalue) + ' ' + self.eukModelFile + ' -')
        os.system(pipe + 'hmmsearch --noali --cpu ' + str(threads) + ' -o ' + outputPrefix + '.lsu.euk.rev_comp.txt --tblout ' + outputPrefix + '.lsu.euk.table.rev_comp.txt -E ' + str(evalue) + ' ' + self.eukRevCompModelFile + ' -')

        if not self.bQuiet:
            print ''

    def processPairs(self, pairs, threads, evalue, outputDir, sample):
        for i in xrange(0, len(pairs), 2):
            pair1 = pairs[i]
            pair2 = pairs[i+1]

            if not self.bQuiet:
                print 'Identifying LSU sequences in paired-end reads: ' + pair1 + ', ' + pair2

            outputPrefix1 = os.path.join(outputDir, 'extracted_lsu', sample + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')])
            outputPrefix2 = os.path.join(outputDir, 'extracted_lsu', sample + '.' + pair2[pair2.rfind('/')+1:pair2.rfind('.')])

            if not self.bQuiet:
                print '  Processing file: ' + pair1
            self.hmmSearch(pair1, threads, evalue, outputPrefix1)

            if not self.bQuiet:
                print '  Processing file: ' + pair2
            self.hmmSearch(pair2, threads, evalue, outputPrefix2)

            # reads hits
            hitsBacteria1 = self.getHits(outputPrefix1 + '.lsu.bacteria.table.txt')
            hitsRevCompBacteria1 = self.getHits(outputPrefix1 + '.lsu.bacteria.table.rev_comp.txt')

            hitsArchaea1 = self.getHits(outputPrefix1 + '.lsu.archaea.table.txt')
            hitsRevCompArcheae1 = self.getHits(outputPrefix1 + '.lsu.archaea.table.rev_comp.txt')

            hitsEuk1 = self.getHits(outputPrefix1 + '.lsu.euk.table.txt')
            hitsRevCompEuk1 = self.getHits(outputPrefix1 + '.lsu.euk.table.rev_comp.txt')

            hitsBacteria2 = self.getHits(outputPrefix2 + '.lsu.bacteria.table.txt')
            hitsRevCompBacteria2 = self.getHits(outputPrefix2 + '.lsu.bacteria.table.rev_comp.txt')

            hitsArchaea2 = self.getHits(outputPrefix2 + '.lsu.archaea.table.txt')
            hitsRevCompArcheae2 = self.getHits(outputPrefix2 + '.lsu.archaea.table.rev_comp.txt')

            hitsEuk2 = self.getHits(outputPrefix2 + '.lsu.euk.table.txt')
            hitsRevCompEuk2 = self.getHits(outputPrefix2 + '.lsu.euk.table.rev_comp.txt')

            # combine hits
            hits1 = hitsBacteria1.union(hitsRevCompBacteria1).union(hitsArchaea1).union(hitsRevCompArcheae1).union(hitsEuk1).union(hitsRevCompEuk1)
            hits2 = hitsBacteria2.union(hitsRevCompBacteria2).union(hitsArchaea2).union(hitsRevCompArcheae2).union(hitsEuk2).union(hitsRevCompEuk2)

            if not self.bQuiet:
                print '  Hits in ' + pair1 + ': ' + str(len(hits1))
                print '    Fwd. bacterial hits: ' + str(len(hitsBacteria1))
                print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria1))
                print '    Fwd. archaeal hits: ' + str(len(hitsArchaea1))
                print '    Rev. comp. archaeal hits: ' + str(len(hitsRevCompArcheae1))
                print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk1))
                print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk1))
                print ''

                print '  Hits in ' + pair2 + ': ' + str(len(hits2))
                print '    Fwd. bacterial hits: ' + str(len(hitsBacteria2))
                print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria2))
                print '    Fwd. archaeal hits: ' + str(len(hitsArchaea2))
                print '    Rev. comp. archaeal hits: ' + str(len(hitsRevCompArcheae2))
                print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk2))
                print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk2))
                print ''

            # extract reads with hits
            if not self.bQuiet:
                print '  Extracting putative 16S/18S reads:'
            hitUnion = hits1.union(hits2)

            seqs1 = extractSeqs(pair1, hitUnion)
            seqs2 = extractSeqs(pair2, hitUnion)

            # create file with all 16S/18S sequences
            allSeqFile = outputPrefix1 + '.lsu.all.SSU.fasta'
            fout = open(allSeqFile, 'w')
            for seqId in hits1:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')

            for seqId in hits2:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')

            fout.close()

            # create paired-end files where at least one read maps to a 16S/18S
            pair1FileUnion = outputPrefix1 + '.lsu.union.SSU.fasta'
            fout = open(pair1FileUnion, 'w')
            for seqId in hitUnion:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')
            fout.close()

            pair2FileUnion = outputPrefix2 + '.lsu.union.SSU.fasta'
            fout = open(pair2FileUnion, 'w')
            for seqId in hitUnion:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')
            fout.close()

            # create paired-end files where at least one read maps to a 16S/18S
            hitIntersection = hits1.intersection(hits2)
            pair1FileIntersect = outputPrefix1 + '.lsu.intersect.SSU.fasta'
            fout = open(pair1FileIntersect, 'w')
            for seqId in hitIntersection:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')
            fout.close()

            pair2FileIntersect = outputPrefix2 + '.lsu.intersect.SSU.fasta'
            fout = open(pair2FileIntersect, 'w')
            for seqId in hitIntersection:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')
            fout.close()

            # create file where only one read maps to a 16S/18S
            hitDiff1 = hits1.difference(hits2)
            diffFile = outputPrefix1 + '.lsu.difference.SSU.fasta'
            fout = open(diffFile, 'w')
            for seqId in hitDiff1:
                r = seqs1[seqId][0]     # strip read identifier as these will be mapped as singletons
                if '/' in r:
                    r = r[0:r.rfind('/')]
                fout.write('>' + r + '\n')
                fout.write(seqs1[seqId][1] + '\n')

            hitDiff2 = hits2.difference(hits1)
            for seqId in hitDiff2:
                r = seqs2[seqId][0]
                if '/' in r:
                    r = r[0:r.rfind('/')]
                fout.write('>' + r + '\n')
                fout.write(seqs2[seqId][1] + '\n')
            fout.close()

            if not self.bQuiet:
                print '    Hits to left reads: ' + str(len(hits1))
                print '    Hits to right reads: ' + str(len(hits2))
                print '    Pairs with at least one read identified as 16S/18S: ' + str(len(hitUnion)) + ' pairs'
                print '    Pairs with both read identified as 16S/18S: ' + str(len(hitIntersection)) + ' pairs'
                print '    Pairs with only one read identified as 16S/18S: ' + str(len(hitDiff1) + len(hitDiff2)) + ' reads'
                print ''
                print '    All identified 16S reads: ' + allSeqFile
                print '    Pairs with at least one read identified as 16S/18S written to: '
                print '      ' + pair1FileUnion
                print '      ' + pair2FileUnion
                print '    Pairs with both read identified as 16S/18S written to: '
                print '      ' + pair1FileIntersect
                print '      ' + pair2FileIntersect
                print '    Pairs with only one read identified as 16S/18S written to: '
                print '      ' + diffFile
                print ''

    def processSingles(self, singles, threads, evalue, outputDir, sample):
        for i in xrange(0, len(singles)):
            seqFile = singles[i]

            if not self.bQuiet:
                print 'Identifying 16S/18S sequences in single-end reads: ' + seqFile

            outputPrefix = os.path.join(outputDir, 'extracted_lsu', + sample + '.' + seqFile[seqFile.rfind('/')+1:seqFile.rfind('.')])

            self.hmmSearch(seqFile, threads, evalue, outputPrefix)

            # reads hits
            hitsBacteria = self.getHits(outputPrefix + '.lsu.bacteria.table.txt')
            hitsRevCompBacteria = self.getHits(outputPrefix + '.lsu.bacteria.table.rev_comp.txt')

            hitsArchaea = self.getHits(outputPrefix + '.lsu.archaea.table.txt')
            hitsRevCompArcheae = self.getHits(outputPrefix + '.lsu.archaea.table.rev_comp.txt')

            hitsEuk = self.getHits(outputPrefix + '.lsu.euk.table.txt')
            hitsRevCompEuk = self.getHits(outputPrefix + '.lsu.euk.table.rev_comp.txt')

            hits = hitsBacteria.union(hitsRevCompBacteria).union(hitsArchaea).union(hitsRevCompArcheae).union(hitsEuk).union(hitsRevCompEuk)

            if not self.bQuiet:
                print '  Hits in ' + seqFile
                print '    Fwd. bacterial hits: ' + str(len(hitsBacteria))
                print '    Rev. comp. bacterial hits: ' + str(len(hitsRevCompBacteria))
                print '    Fwd. archaeal hits: ' + str(len(hitsArchaea))
                print '    Rev. comp. archaeal hits: ' + str(len(hitsRevCompArcheae))
                print '    Fwd. eukaryotic hits: ' + str(len(hitsEuk))
                print '    Rev. comp. eukaryotic hits: ' + str(len(hitsRevCompEuk))
                print ''
                print '  Identified 16S/18S reads: ' + str(len(hits)) + ' reads'
                print ''

            # extract reads with hits
            seqs = extractSeqs(seqFile, hits)

            # create file with all 16S sequences
            allSeqFile = outputPrefix + '.lsu.SSU.fasta'
            fout = open(allSeqFile, 'w')
            for seqId in hits:
                r = seqs[seqId][0]     # strip read identifier as these will be mapped as singletons
                if '/' in r:
                    r = r[0:r.rfind('/')]
                fout.write('>' + r + '\n')
                fout.write(seqs[seqId][1] + '\n')

            fout.close()

            if not self.bQuiet:
                print '  Identified 16S/18 reads written to: ' + allSeqFile
                print ''

    def run(self, configFile, threads, evalue, bQuiet):
        rc = ReadConfig()
        projectParams, sampleParams = rc.readConfig(configFile, outputDirExists = True)
        os.makedirs(projectParams['output_dir'] + 'extracted_lsu')

        self.bQuiet = bQuiet

        for sample in sampleParams:
            pairs = sampleParams[sample]['pairs']
            singles = sampleParams[sample]['singles']

            # identify 16S sequences in paired-end reads
            self.processPairs(pairs, threads, evalue, projectParams['output_dir'], sample)

            # identify 16S sequences in single-end reads
            self.processSingles(singles, threads, evalue, projectParams['output_dir'], sample)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract LSU sequences from metagenomic data using HMMs.",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config_file', help='project config file.')

    parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)
    parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')
    parser.add_argument('-q', '--quiet', help='Surpress all output', action='store_true')

    args = parser.parse_args()

    extractLSU = ExtractLSU()
    extractLSU.run(args.config_file, args.threads, args.evalue, args.quiet)
