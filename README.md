CommunityM
============================

Tools for identifying, classifying, and assembling 16S reads in metagenomic data.


##Overview

CommunityM is a set of scripts aimed at identifying, classifying, and assembling 16S/18S reads within metagenomes and on assembled contigs. The overall goal of CommunityM is to provide a complete picture of the 16S/18S data within your metagenomes and how this relates to your putative genomes (i.e., obtained from GroopM). 


##SOP

Many of the CommunityM scripts assume you have a set of genome bins. These can be obtained with [GroopM](http://minillinim.github.io/GroopM/). Please note that these scripts do not QC your sequences.

**Note:** If you wish to compare community profiles obtained from a metagenome against those obtain with pyrotag data see the section ''Comparing Community Profiles''.

##Identifying and Classifying Binned 16S/18S rRNA Genes

The following scripts identify and classify bacterial, archaeal, and eukaryotic 16S/18S genes on contigs or scaffolds within a set of putative genome. This can often provides a refined estimate of taxonomy. To identify 16S/18S genes run:

<code>
identifyBinned16S.py -t 16 -x fa ./contigs.fna ./bins ./output
</code>

The script ''identifyBinned16S.py'' uses hidden Markov models (HMMs) to identify 16S/18S genes. The number of threads to use is specified with the -t argument. Genes are identified in all contigs/scaffolds within the specified contig file (./contigs.fna in this example). Binning information is taken from the FASTA files within the specified directory (./bins in this example) that have the extension indicated by the -x flag (fa in this example). Results are written to a user specified directory (./output in this example). Use -h for additional information on command line arguments. 

Two files of interest are produced by ''identifyBinned16S.py'':
  * identified16S.tsv: a summary of the identified 16S/18S genes
  * identified16S.fna: a FASTA file containing the identified 16S/18S genes

A number of additional intermediate files are written to the output directory and may be of interest to user familiar with HMMs.

The 16S/18S genes in ''identified16S.fna'' can be taxonomically classified using the ''classifyFastaWithNB.py'' script:

<code>
classifyFastaWithNB.py -t 16 --db GG97 ./output ./output/identified16S.fna
</code>

Genes are classified with [mothur's](http://www.mothur.org/wiki/Classify.seqs) naive Bayes classifier. The reference database to classify against is specified with the ''db'' argument. I recommend using the GreenGenes DB clustered at 97% identitiy (GG97) unless you have Eukaryotic genes to classify in which case the SILVA98 DB is more appropriate. An output directory must be specified (./output in this example) along with one or more fasta files (./output/identified16S.fna in this example). The classification of ''identified16S.fna'' will be written to:
  * identified16S.fna.classified.16S.tsv: taxonomic classification of each sequence in the FASTA file

See the [mothur](http://www.mothur.org/wiki/Classify.seqs) page for information on interpreting the support values provided in parentheses after each classified rank. My rule-of-thumb is to ignore anything below 60 and treat anything between 60-80 as highly putative.

##Identifying and Classifying 16S/18S Reads (config file approach)

**NB:** The notes below are for when using the config file approach. To run each sample on an independent (& scriptable) command line, see below.

The above scripts find 16S/18S genes that existing on binned contigs/scaffolds. Here we are interested in trying to do a more focused search for all 16S/18S reads in a set of metagenome(s) with the goal of producing community profiles or assembling de novo 16S/18S sequences. I have had success with this on an enrichment culture characterized with 2x100bp reads (500bp insert size) (i.e., Steve's coal community). These scripts make use of a configuration file for specifying the sequence files associate with each sample in your project. An example configuration file can be found at /srv/whitlam/bio/apps/12.04/sw/communitym/0.0.10/test/test.cfg:

```
[project]
output_dir = ./test_output

[sample1] <- Note: this must be labelled as sampleX, where X is an integer!
name = Sample 1
edit_dist = 0.1
min_align_len = 0.9
pairs = ./sample1/pair.1.fna, ./sample1/pair.2.fna
singles = 

[sample2]
name = Sample 2
edit_dist = 0.1
min_align_len = 0.9
pairs = ./sample2/pair.1.fna, ./sample2/pair.2.fna
singles = ./sample2/single1.fna, ./sample2/single2.fna, ./sample2/single3.fna

[sample3]
name = pyrotags
edit_dist = 0.1
min_align_len = 0.9
pairs = 
singles = pyrotags.fna
``` 

The config file specified properties for the entire project along with the sequence files for each sample. The sample section headers in square brackets must be labelled as sampleX, where X is an integer (no spaces!). Samples are given an identifying name and any paired-end or single-ended sequence files specified. Input files can be compressed with gzip. In general, I recommend ignoring singletons (i.e., merged pairs) and just processing the paired reads. The parameters ''edit_dist'' and ''min_align_len'' must also be given for each sample. These parameters are used by the ''classifyBWA_16S.py'' script which classifies reads by mapping them with BWA-mem to a database of 16S sequences. These parameters indicate the maximum edit distance and minimum alignment length (as a percentage of sequence length) required to accept a mapping, respectively. 

To identify 16S/18S reads run:
<code>
extractHMM_16S.py -t 16 project.cfg
</code> 

This script uses HMMs to identify 16S/18S reads within all the samples specified in the provided configuration file. Results will be written to the output directory specified in the configuration file. Identified 16S/18S reads can then be classified by mapping them to a reference DB with BWA:
<code>
classifyBWA_16S.py -t 16 project.cfg GG97
</code>

This maps each identified 16S/18S read to the SILVA DB or the GreenGenes DB clustered at either 94, 97, or 99%. Any reads that fail to map are marked as 'unmapped' in the classification file. A read fails to map if BWA-mem reports it as unmapped, or if it fails the maximum edit distance and minimum alignment length settings specified in the configuration file. 

Alternatively, the script ''classifyNB_16S.py'' can be used to classify reads with [mothur's](http://www.mothur.org/wiki/Classify.seqs) naive Bayes classifier. This gives an alternative classification and is provided for convenience, but is not supported by the ''buildTable.py'' script detailed below. I believe the BWA approach is a bit more sensible for short 16S/18S reads. 

Note that the GreenGenes database only has 16S sequences so any 18S sequences in your file will likely get marked as unmapped or unclassified. If you suspect you have an appreciable number of 18S reads use the SILVA database.

At this point, you can examine the classification of these reads using the scripts ''buildTable.py''. This script has a number of arguments that influence the resulting community profile. The ''-r'' arguments specify the desired taxonomic rank of the OTU table. By default, results are given as absolute counts. If relative abundances are desired, use the ''-m rel'' argument. Perhaps the most critical decision is how paired-end reads should be treated. By default, each pair contributes only a single classification and the lowest-common ancestor in common between the two classifications is used. Specifying the ''-p'' argument causes each read in a pair to be treated separately. Currently, I am advocating treating paired reads as a true pairs. You may also include singletons (e.g., merged pairs or pairs where only one end was identified as being from the SSU gene) using the '-s' flag. At present, I do not recommend using singletons as there mapping is less accurate. Finally, the ''-u'' argument can be used to ignore unmapped reads in the resulting OTU table. I suggest looking at your OTU tables both with and without unmapped reads considered. Typical usage:

<code>
buildTable.py -r Class project.cfg myTable.tsv
</code>

Use ''-h'' to get a full list of arguments supported by the ''buildTable.py'' scripts. If all you are after is community profiles your done! If you are interested in trying to assemble these into full length 16S sequences check out the 'Assembling 16S/18S Reads' section below.

##Identifying and Classifying 16S/18S Reads (command-line approach)

To run paired-end reads through CommunityM using only command line parameters:
<code>
community_profile.py -r GG97 -p reads.1.fq.gz,reads.2.fq.gz -o otu_table.tsv
</code>
This command line can then be parallelised with e.g. GNU parallel. Below assumes the bam files are named with the reads e.g. 'reads.bam', '/path/to/reads/reads.1.fq.gz' and '/path/to/reads/reads.2.fq.gz'
<code>
ls *.bam |parallel community_profile.py -r GG97 -p /path/to/reads/{}.1.fq.gz,/path/to/reads/{}.2.fq.gz -o {}.tsv
</code>

##Generating Krona plots 
To make a krona out of communityM tsv files:
<code>
krona_from_community_profiles.py -i readset1.tsv,readset2.tsv -o krona.html
</code>
Then open krona.html in your web browser.

##Phylogenetic beta-diversity

The script expressBetaDiversity.py can be used to calculate unweighted and weighted Soergel dissimilarity values between your communities. This script requires an OTU table built with buildTable.py (see above) with the rank (-r) set to SEQ_ID. Phylogenetic beta-diversity is calculated with respect to a reference tree. You can specify one of the GG reference trees. I recommend GG97. Please note that reads which can't be mapped are effectively ignored when calculating phylogenetic beta-diversity. If you have large skews in the number of unmapped sequences between your samples this needs to be considered. To calculate phylogenetic beta-diversity run:
<code>
expressBetaDiversity.py my_otu_table.tsv GG97 ./output
</code>

The output directory will contain 6 files:
  * otu_table.ebd.tsv: a reformatted version of your OTU table (you can ignore this)
  * sample_counts.tsv: number of reads mapped to each of your samples
  * soergel.diss/usoergel.diss: Phylip-style dissimilarity matrix indicating the Soergel and unweighted Soergel values between your communities
  * soergel.tre/usoergel.tre: a UPGMA tree depicting the clustering of your samples based on the Soergel and unweighted Soergel values

The *.tre files can be viewed with any standard tree viewer (e.g., [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). The dissimilarity matrices can be used to generate PCoA plots. If you have a workflow to do this, please add in the details here. 

For more information of phylogenetic-beta diversity, feel free to visit the desk of Donovan. You can also read my ramblings on the topic:

  * [Measures of phylogenetic differentiation provide robust and complementary insights into microbial communities](http://www.nature.com/ismej/journal/v7/n1/full/ismej201288a.html)

##Assembling 16S/18S Reads

Once you have identified and classified your 16S/18S reads, you can try to assemble them. Not all identified 16S/18S reads are suitable for assembly. The strategy taken here is similar to AmpliShot. Only 16S/18S reads that map to similar reference genes are grouped together for assembly:

<code>
identifyRecoverable16S.py project.cfg 97
</code>

This script also has a number of parameters that can be set to influence which reads will be assembled together. In particular, the sequence identity between reference genes required for reads to be clustered together for assembly is set with -i (default = 0.03). By default, only pairs where both ends map to reference genes within the specified identity threshold will be retained for assembly. This can be relaxed by using the -p flag which indicates each end of a pair should be treated separately and the -s flag which specify that single ended reads should also be used. 

These sets of 16S/18S reads can then be assembled with either [SPAdes](http://bioinf.spbau.ru/spades), [Ray](http://denovoassembler.sourceforge.net), or CLC:

<code>
assemblePutative16S_spades.py -t 16 project.cfg
</code>

In my experience, Ray often fails to properly assemble 16S/18S reads and SPAdes will fail whenever the insert size cannot be estimated from the data. My current strategy is to run the SPAdes assembler and see how many putative 16S/18S genes fail to assemble (this is indicates by the script). If a lot of putative 16S/18S fail to assembly, CLC can often be used to obtain improved results. Unfortunately, assembly in CLC must be done manually be loading each of the <GG Id>.1.fasta and <GG Id>.2.fasta files in your <output dir>/putative16S folder into CLC and doing a de novo assembly for each pair. I expect the next version of SPAdes will allow one to input an estimate of the insert size which should drastically reduce the number of cases where it failed to produce an assembly.

We are primarily interested in 16S/18S assemblies that are not already on an assembled contig/scaffold. To determine which of these de novo 16S/18S assemblies match those on our previously assembled contigs/scaffolds a simple blastn can be performed:

<code>
module load blast
formatdb -i contigs.fna -p F
blastall -p blastn -i clc_assembly.fna -d contigs.fna -m 9 -e 1e-10 > blast.tsv
</code>

This will report all hits between the 16S/18S genes in clc_assembly.fasta and the contigs/scaffolds in contigs.fna. The results in blast.tsv can be examined to determine which of the de novo 16S/18S genes are located on existing assemblies.

As a final step, paired end data where at least one end was identified as being on a 16S gene can be used to try and link de novo 16S assemblies to binned contigs:

<code>
linkBinsTo16S.py -t 16 project.cfg contigs.fna clc_assembly.fna ./bins
</code>

The ''linkBinsTo16S.py'' script requires the configuration file, a FASTA file containing your assembles, a FASTA file containing the de novo assembled 16S/18S genes, and the directory containing your putative genomes. This script produces an output table indicating any unbinned contigs/scaffolds (including your de novo assembled 16S/18S sequences) that can be linked to binned contigs/scaffolds. The format of this table is a bit confusing so just come talk to me. :)

To do: write script to link de novo 16S assemblies to binned contigs via mate pair data. If you need this, let me know so I have some motivation to write this script! :)

##Comparing Community Profiles

CommunityM is flexible in terms on the input files it can process. In particular, you can specify a pyrotag sequence file as a separate sample if you wish to process it in a manner identical to the metagenomic data in order to compare community profiles. An example configuration file might be as follows:

```
[project]
output_dir = ./test_output

[sample1]
name = metagenomic
edit_dist = 0.15
min_align_len = 0.85
pairs = ./sample1/pair.1.fna, ./sample1/pair.2.fna
singles = ./sample1/single1.fna, ./sample1/single2.fna, ./sample1/single3.fna

[sample2]
name = pyrotag
edit_dist = 0.15
min_align_len = 0.85
pairs = 
singles = my_pyrotags.fna
``` 

Recall that you need to perform QC on your pyrotags before processing them with CommunityM. The [http://qiime.org/scripts/split_libraries.html](split_libraries.py) script in QIIME can handle most QC jobs. I also recommend that you make sure the HMMs used by CommunityM correctly identify all your pyrotags as being 16S. To do this, compare the number of sequences in your pyrotag file (my_pyrotags.fna in this example) to the number of sequences in the output of ''extractHMM_16S.py'' (./test_output/extracted/sample2.my_pyrotags.SSU.fasta in this example). If these files don't contain the same number of sequences (or at least nearly the same), then something is up. Either the QC of your pyrotags is insufficient or the HMMs are inadequate. If this happens, please let me know so I have a better sense of how sensitive the HMMs are and/or we can figure out if it is a QC issue.
