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
Useful methods for reading, writing, and processing sequence files.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import gzip

def canonicalHeader(header):
  lineSplit = header[1:].rstrip().split()
  if lineSplit[0][-2] == '/':
    # parse headers with read specified after backslash (e.g., seq001/1)
    readId = lineSplit[0][0:-2]
    readNum = lineSplit[0][-1]
  elif len(lineSplit) >= 2:
    # parse headers with read specified at start of second record (e.g., seq001 1:X:0)
    readId = lineSplit[0]
    readNum = int(lineSplit[1][0])
  else:
    # assume a singleton
    readId = lineSplit[0]
    readNum = None

  readId = readId.replace(':', '_') # follow mothur convention
  universalHeader = '>' + readId
  if readNum != None:
    universalHeader += '/' + str(readNum)

  return universalHeader

def readFasta(fastaFile):
  try:
    fh = file(fastaFile)
  except IOError:
    print "File '" + fastaFile + "' does not exist."
    sys.exit()

  seqs = {}
  for line in fh:
    if line.startswith('>'):
      universalHeader = canonicalHeader(line)
      seqId = universalHeader[1:]
      seqs[seqId] = ''
    else:
      seqs[seqId] += line.rstrip('\n').rstrip('*')

  return seqs

class FastqRecord:
  def __init__(self, seqId = ''):
    self.seqId = seqId
    self.seq = ''
    self.quality = ''

  def asString(self):
    return self.seqId + '\n' + self.seq + '\n' + '+\n' + self.quality + '\n'

def readFastq(fastqFile, seqs):
  try:
    fh = file(fastqFile)
  except IOError:
    print "File '" + fastqFile + "' does not exist."
    sys.exit()

  seqs = {}
  lineNum = 0
  for line in fh:
    if lineNum == 0:
      seqId = line[1:].rstrip('\n')
      seqs[seqId] = FastqRecord(seqId)
    elif lineNum == 1:
      seqs[seqId].seq = line.rstrip('\n')
    elif lineNum == 3:
      seqs[seqId].quality = line.rstrip('\n')

    lineNum += 1
    lineNum %= 4

  return seqs

def extractSeqs(f, seqIdsToExtract):
  ''' Sequence ids indicate the sequence from which reads were taken (i.e., do NOT contain a /1 or /2). '''
  if f.endswith('gz'):
    if '.fastq' in f or '.fq' in f:
      return extractSeqsFastqGz(f, seqIdsToExtract)
  else:
    return extractSeqsFasta(f, seqIdsToExtract)

  print '[Error] Unknown file type'
  sys.exit()

def extractSeqsFastqGz(fastqFile, seqIdsToExtract):
  seqs = {}

  lineNum = 0
  for line in gzip.open(fastqFile):
    if lineNum == 0:
      readId = line[1:].split()[0].rstrip()
      seqId = readId

      index = readId.rfind('/')
      if index != -1:
        seqId = readId[0:readId.rfind('/')]

      bKeep = (seqId in seqIdsToExtract)
    elif lineNum == 1:
      if bKeep:
        seqs[seqId] = [readId, line.rstrip()]

    lineNum += 1
    lineNum %= 4

  return seqs

def extractSeqsFasta(fastaFile, seqIdsToExtract):
  seqs = {}

  bKeep = False
  for line in open(fastaFile):
    if line[0] == '>':
      if bKeep:
        seqs[seqId] = [readId, seq]

      readId = line[1:].split()[0].rstrip()
      seqId = readId

      index = readId.rfind('/')
      if index != -1:
        seqId = readId[0:readId.rfind('/')]

      bKeep = (seqId in seqIdsToExtract)
      seq = ''
    elif bKeep:
      seq += line.rstrip()

  if bKeep:
    seqs[seqId] = [readId, seq]

  return seqs


def extractReads(f, readIdsToExtract):
  ''' Read ids identify individual reads taken from a sequence (i.e., contain a /1 or /2). '''
  if f.endswith('gz'):
    if '.fastq' in f or '.fq' in f:
      return extractReadsFastqGz(f, readIdsToExtract)
  else:
    return extractReadsFasta(f, readIdsToExtract)

  print '[Error] Unknown file type'
  sys.exit()

def extractReadsFastqGz(fastqFile, readIdsToExtract):
  reads = {}
  lineNum = 0
  for line in gzip.open(fastqFile):
    if lineNum == 0:
      readId = line[1:].split()[0].rstrip()
      bKeep = (readId in readIdsToExtract)
    elif lineNum == 1:
      if bKeep:
        reads[readId] = line.rstrip()

    lineNum += 1
    lineNum %= 4

  return reads

def extractReadsFasta(fastaFile, readIdsToExtract):
  reads = {}

  bKeep = False
  for line in open(fastaFile):
    if line[0] == '>':
      if bKeep:
        reads[readId] = read

      readId = line[1:].split()[0].rstrip()
      bKeep = (readId in readIdsToExtract)
      read = ''
    elif bKeep:
      read += line.rstrip()

  if bKeep:
    reads[readId] = read

  return reads
