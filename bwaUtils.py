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

import os

def mapPair(db, file1, file2, bamPrefix, threads):
    print 'Processing pair: ' + file1 + ', ' + file2

    if not os.path.exists(db + '.bwt'):
        print ''
        print 'Indexing database file: '
        os.system('bwa index ' + db)

    print ''
    print 'Aligning pairs with bwa-mem: '
    os.system('bwa mem -a -t ' + str(threads) + ' ' + db + ' ' + file1 + ' ' + file2 + ' | samtools view -Subh - | samtools sort - ' + bamPrefix)

    print ''
    print 'Indexing BAM file.'
    os.system('samtools index ' + bamPrefix + '.bam')
    print ''

def mapSingle(db, filename, bamPrefix, threads):
    print 'Processing single-ended read file: ' + filename

    if not os.path.exists(db + '.bwt'):
        print ''
        print 'Indexing database file: '
        os.system('bwa index ' + db)

    print ''
    print 'Aligning with bwa-mem: '
    os.system('bwa mem -a -t ' + str(threads) + ' ' + db + ' ' + filename + ' | samtools view -Subh - | samtools sort - ' + bamPrefix)

    print ''
    print 'Indexing BAM file.'
    os.system('samtools index ' + bamPrefix + '.bam')
    print ''
