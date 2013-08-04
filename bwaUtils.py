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
Useful methods for interacting with BWA
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
import os
import tempfile

def mapPair(db, file1, file2, bamPrefix, threads):
  tmpFD, tmpAlnFile1 = tempfile.mkstemp()
  tmpFD, tmpAlnFile2 = tempfile.mkstemp()

  print 'Processing pair: ' + file1 + ', ' + file2

  print ''
  print 'Aligning pairs: '
  os.system('bwa aln -t ' + str(threads) + ' ' + db + ' ' + file1 + ' > ' + tmpAlnFile1)
  os.system('bwa aln -t ' + str(threads) + ' ' + db + ' ' + file2 + ' > ' + tmpAlnFile2)

  print ''
  print 'Mapping reads:'
  os.system('bwa sampe ' + db + ' ' + tmpAlnFile1 + ' ' + tmpAlnFile2 + ' ' + file1 + ' ' + file2 + '| samtools view -SubhF 4 - | samtools sort - ' + bamPrefix)

  print ''
  print 'Indexing BAM file.'
  os.system('samtools index ' + bamPrefix + '.bam')
  print ''

  os.remove(tmpAlnFile1)
  os.remove(tmpAlnFile2)

def mapSingle(db, filename, bamPrefix, threads):
  tmpFD, tmpAlnFile = tempfile.mkstemp()

  print 'Processing single-ended read file: ' + filename

  print ''
  print 'Aligning reads: '
  os.system('bwa aln -t ' + str(threads) + ' ' + db + ' ' + filename + ' > ' + tmpAlnFile)

  print ''
  print 'Mapping reads:'
  os.system('bwa samse ' + db + ' ' + tmpAlnFile + ' ' + filename + '| samtools view -SubhF 4 - | samtools sort - ' + bamPrefix)

  print ''
  print 'Indexing BAM file.'
  os.system('samtools index ' + bamPrefix + '.bam')
  print ''

  os.remove(tmpAlnFile)
