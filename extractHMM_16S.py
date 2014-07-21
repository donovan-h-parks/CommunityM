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
Extract 16S/18S sequences from metagenomic data using HMMs.
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
import multiprocessing as mp

from readConfig import ReadConfig
from seqUtils import extractSeqs

class HitRecord(object):
    def __init__(self):
        self.pair1 = None
        self.pair2 = None
        
        self.hits1 = None
        self.hitsBacteria1 = None
        self.hitsArchaea1 = None
        self.hitsEuk1 = None
        
        self.hits2 = None
        self.hitsBacteria2 = None
        self.hitsArchaea2 = None
        self.hitsEuk2 = None
        
        self.uniqueBacterialHits = None
        self.uniqueArchaealHits = None
        self.uniqueEukHits = None
        
        self.hitUnion = None
        self.hitIntersect = None
        self.hitDiff = None
        
        self.allSeqFile = None
        self.pair1FileUnion = None
        self.pair2FileUnion = None
        self.pair1FileIntersect = None
        self.pair2FileIntersect = None
        self.diffFile = None

class Extract16S(object):
    def __init__(self):
        self.bQuiet = False
        
        self.bacteriaModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/SSU_bacteria.hmm'
        self.archaeaModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/SSU_archaea.hmm'
        self.eukModelFile = '/srv/whitlam/bio/db/communitym/ssu_hmm/SSU_euk.hmm'
        
    def hmmSearch(self, seqFile, evalue, threadsPerSample, outputPrefix):
        if seqFile.endswith('gz'):
            pipe = 'zcat ' + seqFile + ' | '
        else:
            pipe = 'cat ' + seqFile + ' | '

        if '.fq' in seqFile or '.fastq' in seqFile:
            pipe += 'fastq2fasta - | '

        os.system(pipe + 'nhmmer --cpu ' + str(threadsPerSample) + ' --noali -o ' + outputPrefix + '.bacteria.txt --tblout ' + outputPrefix + '.bacteria.table.txt -E ' + str(evalue) + ' ' + self.bacteriaModelFile + ' - > /dev/null')
        os.system(pipe + 'nhmmer --cpu ' + str(threadsPerSample) + ' --noali -o ' + outputPrefix + '.archaea.txt --tblout ' + outputPrefix + '.archaea.table.txt -E ' + str(evalue) + ' ' + self.archaeaModelFile + ' - > /dev/null')
        os.system(pipe + 'nhmmer --cpu ' + str(threadsPerSample) + ' --noali -o ' + outputPrefix + '.euk.txt --tblout ' + outputPrefix + '.euk.table.txt -E ' + str(evalue) + ' ' + self.eukModelFile + ' - > /dev/null')

    def __getHits(self, hitTable, alignLenThreshold):
        seqIds = set()
        for line in open(hitTable):
            if line[0] == '#' or line.strip() == '':
                continue

            lineSplit =  line.split()
            
            seqId = lineSplit[0].split('/')[0]   # remove read identifier
            alignLen = abs(int(lineSplit[6]) - int(lineSplit[7]))  
            seqLen = int(lineSplit[10])
            
            if alignLen >= alignLenThreshold * seqLen:
                seqIds.add(seqId)

        return seqIds

    def __processPairs(self, pairs, evalue, alignLenThreshold, outputDir, sample, threadsPerSample, queueOut):
        for i in xrange(0, len(pairs), 2):
            pair1 = pairs[i]
            pair2 = pairs[i+1]

            outputPrefix1 = os.path.join(outputDir,'extracted',sample + '.' + pair1[pair1.rfind('/')+1:pair1.rfind('.')])
            outputPrefix2 = os.path.join(outputDir,'extracted',sample + '.' + pair2[pair2.rfind('/')+1:pair2.rfind('.')])

            self.hmmSearch(pair1, evalue, threadsPerSample, outputPrefix1)
            self.hmmSearch(pair2, evalue, threadsPerSample, outputPrefix2)

            # reads hits
            hitsBacteria1 = self.__getHits(outputPrefix1 + '.bacteria.table.txt', alignLenThreshold)
            hitsArchaea1 = self.__getHits(outputPrefix1 + '.archaea.table.txt', alignLenThreshold)
            hitsEuk1 = self.__getHits(outputPrefix1 + '.euk.table.txt', alignLenThreshold)
            
            hitsBacteria2 = self.__getHits(outputPrefix2 + '.bacteria.table.txt', alignLenThreshold)
            hitsArchaea2 = self.__getHits(outputPrefix2 + '.archaea.table.txt', alignLenThreshold)
            hitsEuk2 = self.__getHits(outputPrefix2 + '.euk.table.txt', alignLenThreshold)

            # combine hits
            hits1 = hitsBacteria1.union(hitsArchaea1).union(hitsEuk1)
            hits2 = hitsBacteria2.union(hitsArchaea2).union(hitsEuk2)

            # extract reads with hits
            hitUnion = hits1.union(hits2)

            seqs1 = extractSeqs(pair1, hitUnion)
            seqs2 = extractSeqs(pair2, hitUnion)

            # create file with all 16S/18S sequences
            allSeqFile = outputPrefix1 + '.all.SSU.fasta'
            fout = open(allSeqFile, 'w')
            for seqId in hits1:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')

            for seqId in hits2:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')

            fout.close()

            # create paired-end files where at least one read maps to a 16S/18S
            pair1FileUnion = outputPrefix1 + '.union.SSU.fasta'
            fout = open(pair1FileUnion, 'w')
            for seqId in hitUnion:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')
            fout.close()

            pair2FileUnion = outputPrefix2 + '.union.SSU.fasta'
            fout = open(pair2FileUnion, 'w')
            for seqId in hitUnion:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')
            fout.close()

            # create paired-end files where both read maps to a 16S/18S
            hitIntersection = hits1.intersection(hits2)
            pair1FileIntersect = outputPrefix1 + '.intersect.SSU.fasta'
            fout = open(pair1FileIntersect, 'w')
            for seqId in hitIntersection:
                fout.write('>' + seqs1[seqId][0] + '\n')
                fout.write(seqs1[seqId][1] + '\n')
            fout.close()

            pair2FileIntersect = outputPrefix2 + '.intersect.SSU.fasta'
            fout = open(pair2FileIntersect, 'w')
            for seqId in hitIntersection:
                fout.write('>' + seqs2[seqId][0] + '\n')
                fout.write(seqs2[seqId][1] + '\n')
            fout.close()

            # create file where only one read maps to a 16S/18S
            hitDiff1 = hits1.difference(hits2)
            diffFile = outputPrefix1 + '.difference.SSU.fasta'
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

            # gather output data
            bacHits = hitsBacteria1.union(hitsBacteria2)
            arHits = hitsArchaea1.union(hitsArchaea2)
            eukHits = hitsEuk1.union(hitsEuk2)
            
            hitRecord = HitRecord()
            hitRecord.pair1 = pair1
            hitRecord.pair2 = pair2
            
            hitRecord.hits1 = len(hits1)
            hitRecord.hitsBacteria1 = len(hitsBacteria1)
            hitRecord.hitsArchaea1 = len(hitsArchaea1)
            hitRecord.hitsEuk1 = len(hitsEuk1)
            
            hitRecord.hits2 = len(hits2)
            hitRecord.hitsBacteria2 = len(hitsBacteria2)
            hitRecord.hitsArchaea2 = len(hitsArchaea2)
            hitRecord.hitsEuk2 = len(hitsEuk2)
            
            hitRecord.uniqueBacterialHits = len(bacHits.difference(arHits.union(eukHits)))
            hitRecord.uniqueArchaealHits = len(arHits.difference(bacHits.union(eukHits)))
            hitRecord.uniqueEukHits = len(eukHits.difference(bacHits.union(arHits)))
            
            hitRecord.hitUnion = len(hitUnion)
            hitRecord.hitIntersect = len(hitIntersection)
            hitRecord.hitDiff = len(hitDiff1) + len(hitDiff2)
            
            hitRecord.allSeqFile = allSeqFile
            hitRecord.pair1FileUnion = pair1FileUnion
            hitRecord.pair2FileUnion = pair2FileUnion
            hitRecord.pair1FileIntersect = pair1FileIntersect
            hitRecord.pair2FileIntersect = pair2FileIntersect
            hitRecord.diffFile = diffFile
               
            queueOut.put(hitRecord)

    def __processSingles(self, singles, evalue, alignLenThreshold, outputDir, sample, threadsPerSample, queueOut):
        for i in xrange(0, len(singles)):
            seqFile = singles[i]

            outputPrefix = os.path.join(outputDir, 'extracted', sample + '.' + seqFile[seqFile.rfind('/')+1:seqFile.rfind('.')])

            self.hmmSearch(seqFile, evalue, threadsPerSample, outputPrefix)

            # reads hits
            hitsBacteria = self.__getHits(outputPrefix + '.bacteria.table.txt', alignLenThreshold)
            hitsArchaea = self.__getHits(outputPrefix + '.archaea.table.txt', alignLenThreshold)
            hitsEuk = self.__getHits(outputPrefix + '.euk.table.txt', alignLenThreshold)

            hits = hitsBacteria.union(hitsArchaea).union(hitsEuk)

            # extract reads with hits
            seqs = extractSeqs(seqFile, hits)

            # create file with all 16S sequences
            allSeqFile = outputPrefix + '.SSU.fasta'
            fout = open(allSeqFile, 'w')
            for seqId in hits:
                r = seqs[seqId][0]     # strip read identifier as these will be mapped as singletons
                if '/' in r:
                    r = r[0:r.rfind('/')]
                fout.write('>' + r + '\n')
                fout.write(seqs[seqId][1] + '\n')

            fout.close()

            # gather output data
            hitRecord = HitRecord()
            hitRecord.pair1 = seqFile
            hitRecord.pair2 = None
            
            hitRecord.hits1 = len(hits)
            hitRecord.hitsBacteria1 = len(hitsBacteria)
            hitRecord.hitsArchaea1 = len(hitsArchaea)
            hitRecord.hitsEuk1 = len(hitsEuk)
            
            hitRecord.hits2 = None
 
            hitRecord.uniqueBacterialHits = len(hitsBacteria.difference(hitsArchaea.union(hitsEuk)))
            hitRecord.uniqueArchaealHits = len(hitsArchaea.difference(hitsBacteria.union(hitsEuk)))
            hitRecord.uniqueEukHits = len(hitsEuk.difference(hitsBacteria.union(hitsArchaea)))
            
            hitRecord.allSeqFile = allSeqFile
     
            queueOut.put(hitRecord)
                
    def __storeResults(self, numFiles, writerQueue):
        """Write results from threads."""
        
        numFilesProcessed = 0
        while True:
            hitRecord = writerQueue.get(block=True, timeout=None)
            if hitRecord == None:
                break
     
            if not self.bQuiet:
                numFilesProcessed += 1
                print ''
                print '----------------------------------------------------------------------'
                print '  Finished processing %d of %d (%.2f%%) files.' % (numFilesProcessed, numFiles, float(numFilesProcessed)*100/numFiles)
                print ''
                print '  Hits in ' + hitRecord.pair1 + ': ' + str(hitRecord.hits1)
                print '    Bacterial hits: ' + str(hitRecord.hitsBacteria1)
                print '    Archaeal hits: ' + str(hitRecord.hitsArchaea1)
                print '    Eukaryotic hits: ' + str(hitRecord.hitsEuk1)
                print ''

                if hitRecord.pair2:
                    print '  Hits in ' + hitRecord.pair2 + ': ' + str(hitRecord.hits2)
                    print '    Bacterial hits: ' + str(hitRecord.hitsBacteria2)
                    print '    Archaeal hits: ' + str(hitRecord.hitsArchaea2)
                    print '    Eukaryotic hits: ' + str(hitRecord.hitsEuk2)
                    print ''

                print '  Unique bacterial hits: ' + str(hitRecord.uniqueBacterialHits)
                print '  Unique archaeal hits: ' + str(hitRecord.uniqueArchaealHits)
                print '  Unique eukaryotic hits: ' + str(hitRecord.uniqueEukHits)
                
                if hitRecord.pair2:
                    print ''
                    print '  Extracting putative 16S/18S reads:'
                    print '    Hits to left reads: ' + str(hitRecord.hits1)
                    print '    Hits to right reads: ' + str(hitRecord.hits2)
                    print '    Pairs with at least one read identified as 16S/18S: ' + str(hitRecord.hitUnion) + ' pairs'
                    print '    Pairs with both read identified as 16S/18S: ' + str(hitRecord.hitIntersect) + ' pairs'
                    print '    Pairs with only one read identified as 16S/18S: ' + str(hitRecord.hitDiff) + ' reads'
                    print ''
                    print '    All identified 16S reads: ' + hitRecord.allSeqFile
                    print '    Pairs with at least one read identified as 16S/18S written to: '
                    print '      ' + hitRecord.pair1FileUnion
                    print '      ' + hitRecord.pair2FileUnion
                    print '    Pairs with both read identified as 16S/18S written to: '
                    print '      ' + hitRecord.pair1FileIntersect
                    print '      ' + hitRecord.pair2FileIntersect
                    print '    Pairs with only one read identified as 16S/18S written to: '
                    print '      ' + hitRecord.diffFile
                    print ''
                else:
                    print '  Extracting putative 16S/18S reads:'
                    print '  Identified 16S/18 reads written to: ' + hitRecord.allSeqFile
  
    def __runSample(self, evalue, alignLenThreshold, outDir, threadsPerSample, queueIn, queueOut):
        """Run each sample in a separate thread."""
        
        while True:
            sample, pairs, singles = queueIn.get(block=True, timeout=None) 
            if sample == None:
                break 
            
            # identify 16S sequences in single-end reads
            self.__processSingles(singles, evalue, alignLenThreshold, outDir, sample, threadsPerSample, queueOut)
            
            # identify 16S sequences in paired-end reads
            self.__processPairs(pairs, evalue, alignLenThreshold, outDir, sample, threadsPerSample, queueOut)

    def run(self, projectParams, sampleParams, threads, evalue, alignLenThreshold, bQuiet):
        print '\nProcessing %s sample(s).\n' % len(sampleParams)

        os.makedirs(os.path.join(projectParams['output_dir'], 'extracted'))

        self.bQuiet = bQuiet
        
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        numFiles = 0
        for sample in sampleParams:
            pairs = sampleParams[sample]['pairs']
            singles = sampleParams[sample]['singles']
            
            numFiles += len(pairs)/2 + len(singles)
            
            workerQueue.put((sample, pairs, singles))
            
        for _ in range(threads):
            workerQueue.put((None, None, None))
            
        threadsPerSample = max(1, int(threads / len(sampleParams)))
        calcProc = [mp.Process(target = self.__runSample, args = (evalue, alignLenThreshold, projectParams['output_dir'], threadsPerSample, workerQueue, writerQueue)) for _ in range(threads)]
        writeProc = mp.Process(target = self.__storeResults, args = (numFiles, writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put(None)
        writeProc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract 16S/18S sequences from metagenomic data using HMMs.",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config_file', help='project config file.')

    parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)
    parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')
    parser.add_argument('-a', '--align_len', type=float, help='fraction of read that must align for identifying hits', default = '0.5')
    parser.add_argument('-q', '--quiet', help='suppress all output', action='store_true')

    args = parser.parse_args()

    # Read config file
    rc = ReadConfig()
    projectParams, sampleParams = rc.readConfig(args.config_file, outputDirExists = False)

    extract16S = Extract16S()
    extract16S.run(projectParams, sampleParams, args.threads, args.evalue, args.align_len, args.quiet)
