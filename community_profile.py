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
import tempfile
import shutil

# Add the current directory to the python path so local libraries can be used
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import extractHMM_16S
import classifyBWA_16S
import buildTable

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a community profile based on 16S reads')
    parser.add_argument('-p','--paired-reads', help='metagenomic reads to use, comma separated',required=True)
    parser.add_argument('-o','--output_tsv', help='output OTU table file name', required=True)
    parser.add_argument('-r', '--ref_db', help='Reference DB to use for classification (choices: GG94, GG97, GG99, SILVA98)', choices=['GG94', 'GG97', 'GG99', 'SILVA98'],required=True)

    parser.add_argument('--output_dir', help='output directory [default: use temporary direcory, delete afterwards]')
    parser.add_argument('-t', '--threads', help='number of threads', type = int, default = 1)
    parser.add_argument('-q', '--quiet', help='Surpress all output', action='store_true')
    parser.add_argument('-e', '--evalue', help='e-value threshold for identifying hits', default = '1e-5')
    parser.add_argument('--edit_distance', help='edit distance for LCA calculations', type = float, default = 0.15)
    parser.add_argument('--min_align_len', help='minimum alignment length for LCA calculations', type = float, default = 0.85)
    parser.add_argument('-u', '--ignore_unmapped', help='do not consider unmapped reads', action="store_true")
    parser.add_argument('--pairs_as_singles', help='treat paired reads as singletons', action="store_true")
    parser.add_argument('-s', '--singletons', help='use singleton 16S/18S reads identified within paired reads', action="store_true")
    parser.add_argument('-b', '--bootstrap', help='bootstrap threshold required to accept classification (default = 0)', type=int, default=0)
    parser.add_argument('-m', '--mode', help='write values as "rel"ative, "abs"olute or "pre"sence/absense (default = abs)', default="abs")
    parser.add_argument('--rank', help='taxonomic rank of table (choices: Domain, Phylum, Class, Order, Family, Genus, Species, SEQ_ID), (default = SEQ_ID)',
                              choices=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SEQ_ID'], default='SEQ_ID')
    parser.add_argument('-a', '--align_len', type=float, help='fraction of read that must align for identifying hits', default = '0.5')

    args = parser.parse_args()

    # Parameters (simulates the output from ConfigFile class)
    projectParams = {}

    # output directory
    outdir = args.output_dir
    if outdir is None:
        outdir = tempfile.mkdtemp('','CommunityM-')

    projectParams['output_dir'] = outdir

    sample = {}
    splits = args.paired_reads.split(',')
    pair1 = splits[0]
    pair2 = splits[1]
    print "Using read files",pair1,"and",pair2
    sample['pairs'] = splits
    sample['singles'] = []
    if len(sample['pairs']) != 2:
        raise ValueError("Need 2 read files comma separated, found %i",len(splits))
    sample_name = os.path.splitext(os.path.splitext(os.path.basename(pair1))[0])[0]
    sample['name'] = sample_name
    sample['edit_dist'] = args.edit_distance
    sample['min_align_len'] = args.min_align_len

    sampleParams = {}
    sampleParams[sample_name] = sample

    # Extract the 16S reads
    print "Extracting 16S reads.."
    extractor = extractHMM_16S.Extract16S()
    extractor.run(projectParams, sampleParams, args.threads, args.evalue, args.align_len, args.quiet)

    # Classify the reads taxonomically
    print "Classifying reads..."
    classifier = classifyBWA_16S.ClassifyBWA()
    classifier.run(projectParams, sampleParams, args.ref_db, args.threads)

    print "Extracting last common ancestors of pairs..."
    lcaer = buildTable.BuildTable()
    sampleCounts, taxonomy = lcaer.paramsToCountsAndTaxonomy(
        projectParams, sampleParams, args.ignore_unmapped, args.pairs_as_singles, args.singletons, args.bootstrap, args.rank, args.mode)

    # Delete temporary directory if one was created
    if args.output_dir is None:
        shutil.rmtree(outdir)

    # Write out results to file
    fout = open(args.output_tsv, 'w')
    keys = sampleCounts.keys()
    if len(keys) != 1: raise
    key = keys[0]
    fout.write("\t".join(["#OTU_ID", key, 'Consensus Lineage']))
    fout.write("\n")
    i=0
    for taxonomy in sampleCounts[key].keys():
        fout.write("\t".join([str(i), str(sampleCounts[key][taxonomy]), taxonomy]))
        fout.write("\n")
        i += 1

    fout.close()


