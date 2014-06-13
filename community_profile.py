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

import sys
import os
import argparse

# Add the current directory to the python path so local libraries can be used
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import extractHMM_16S
import classifyBWA_16S
import buildTable

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a community profile based on 16S reads')
    parser.add_argument('-p','--paired-reads', help='metagenomic reads to use, comma separated',required=True)
    parser.add_argument('-o','--output_biom', help='output BIOM file', required=True)
    parser.add_argument('--output_dir', help='output directory', default='community_profile')
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 16)
    parser.add_argument('-q', '--quiet', help='Surpress all output', action='store_true')
    parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')
    parser.add_argument('--edit_distance', help='edit distance for LCA calculations', type = float, default = 0.15)
    parser.add_argument('--min_align_len', help='minimum alignment length for LCA calculations', type = float, default = 0.85)
    parser.add_argument('-u', '--ignore_unmapped', help='do not consider unmapped reads', action="store_true")
    parser.add_argument('--pairs_as_singles', help='treat paired reads as singletons', action="store_true")
    parser.add_argument('-s', '--singletons', help='use singleton 16S/18S reads identified within paired reads', action="store_true")
    parser.add_argument('-b', '--bootstrap', help='bootstrap threshold required to accept classification (default = 0)', type=int, default=0)
    parser.add_argument('-a', '--absolute', help='write absolute values instead of relative values', action='store_true')
    parser.add_argument('-r', '--rank', help='taxonomic rank of table (choices: Domain, Phylum, Class, Order, Family, Genus, Species, GG_ID), (default = GG_ID)',
                              choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'GG_ID'], default='GG_ID')


    args = parser.parse_args()

    # Parameters (simulates the output from ConfigFile class)
    projectParams = {}
    projectParams['output_dir'] = args.output_dir

    sample = {}
    splits = args.paired_reads.split(',')
    pair1 = splits[0]
    pair2 = splits[1]
    print "Using read files",pair1,"and",pair2
    sample['pairs'] = splits
    sample['singles'] = []
    if len(sample['pairs']) != 2:
        raise ValueError("Need 2 read files comma separated, found %i",len(splits))
    sample_name = 'sample0'
    sample['name'] = sample_name
    sample['edit_dist'] = args.edit_distance
    sample['min_align_len'] = args.min_align_len

    sampleParams = {}
    sampleParams[sample_name] = sample

    # Extract the 16S reads
    print "Extracting 16S reads.."
    extractor = extractHMM_16S.Extract16S()
    extractor.run(projectParams, sampleParams, args.threads, args.evalue, args.quiet)

    # Classify the reads taxonomically
    print "Classifying reads..."
    classifier = classifyBWA_16S.ClassifyBWA()
    classifier.run(projectParams, sampleParams, '97', args.threads)#TODO: remove the hardcoded 97 here, and add ability to use Silva etc.

    print "Extracting last common ancestors of pairs..."
    lcaer = buildTable.BuildTable()
    lcaer.run(projectParams, sampleParams, args.ignore_unmapped, args.pairs_as_singles, args.singletons, args.bootstrap, args.rank, args.absolute, args.output_biom)
