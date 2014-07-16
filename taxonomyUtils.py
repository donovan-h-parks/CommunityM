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
Useful methods for processing taxonomy strings.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

ranksByLabel = {'Domain':0, 'Phylum':1, 'Class':2, 'Order':3, 'Family':4, 'Genus':5, 'Species':6, 'SEQ_ID':7}
ranksByLevel = {0:'Domain', 1:'Phylum', 2:'Class', 3:'Order', 4:'Family', 5:'Genus', 6:'Species', 7:'SEQ_ID'}
rankPrefixes = {0:'k__', 1:'p__', 2:'c__', 3:'o__', 4:'f__', 5:'g__', 6:'s__', 7:'id__'}

def readTaxonomy(taxonomyFile):
    ggIdToTaxonomy = {}
    for line in open(taxonomyFile):
        lineSplit = line.split('\t')
        
        seqId = lineSplit[0]
        taxonomy = [x.strip() for x in lineSplit[1].split(';')]
        
        if 'id__' not in taxonomy[-1]:
            # missing sequence ID field
            taxonomy += ['id__' + str(seqId)]

        ggIdToTaxonomy[seqId] = taxonomy
        
    return ggIdToTaxonomy

def parseTaxon(taxon):
    if '(' in taxon:
        taxonSplit = taxon.split('(')
        taxonId = taxonSplit[0]
        taxonId = taxonId.strip()
        bootstrapSupport = int(taxonSplit[1][0:taxonSplit[1].find(')')])
    else:
        taxonId = taxon.strip()
        bootstrapSupport = 0

    return taxonId, bootstrapSupport

def LCA(taxonomy1, taxonomy2):
    taxonomy = []
    bUnmapped = False
    for i in xrange(0, len(ranksByLevel)-1):
        t1, b1 = parseTaxon(taxonomy1[i])
        t2, b2 = parseTaxon(taxonomy2[i])

        if t1 != t2:
            if 'unmapped' in t1 or 'unmapped' in t2:
                bUnmapped = True
                taxonomy.append(rankPrefixes[i] + 'unmapped')
            else:
                taxonomy.append(rankPrefixes[i] + 'unresolved_by_lca')
        else:
            if b1 == 0 and b2 == 0:
                taxonomy.append(t1)
            else:
                taxonomy.append(t1 + '(' + str(min(b1, b2)) + ')')

    # return reference sequence id
    if not bUnmapped:
        t1, b1 = parseTaxon(taxonomy1[len(ranksByLevel)-1])
        if b1 == 0:
            taxonomy.append(t1)
        else:
            taxonomy.append(t1 + '(' + str(min(b1, b2)) + ')')
    else:
        taxonomy.append('id__unmapped')

    return taxonomy
