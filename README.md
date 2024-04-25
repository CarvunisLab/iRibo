------------------------------------------------------------------------------

# iRibo

A comprehensive tool for integrating ribosome profiling data to detect genome wide translation.

------------------------------------------------------------------------------

## Table of Contents:
1. Installation and Prerequisites
2. Data Collection and Preprocessing
3. GetCandidateORFs
4. GenerateTranslationProfile
5. GenerateTranslatome

------------------------------------------------------------------------------

## Installation and Prerequisites:

- iRibo Download: 
  Obtain iRibo from https://github.com/CarvunisLab/iRibo. 
  1. Download the zip file.
  2. cd into the directory.
  3. Run make.
  Note: Compiled using g++ compiler and c++17.

- System Requirements:
  - Memory: At least 64GB
  - Storage: At least 100GB

- Dependencies: 
  R version 4.2.2+. Download and documentation are available at https://www.r-project.org/.
  Samtools. Download and documentation are available at https://www.htslib.org/. 

------------------------------------------------------------------------------

## Data Collection and Preprocessing:

iRibo inputs:
1. Genome in FASTA format.
2. Genome annotations in GFF3 or GTF format.
3. Transcriptome (optional) in GFF3 or GTF format.
4. Aligned ribo-seq reads in SAM or BAM format produced by a program like STAR (https://github.com/alexdobin/STAR)

Recommendation: 
- For best alignment, trim low-quality reads and adapters using a program like cutadapt or trim-galore. 
- To find noncanonical translated ORFs, it is ideal to use a comprehensive transcriptome; i.e. a transcriptome that includes transcripts that may not be annnotated genes. 

Note:
- Ensure genome, annotations, and transcriptome have consistent chromosome identifiers, e.g., >chr1 in genome should match chr1 in annotations.

------------------------------------------------------------------------------

## GetCandidateORFs:

Generates a list of candidate ORFs to assess for translation by identifying open reading frames across a genome or transcriptome. Completion of this step will generate three output files:
- all_orfs: A list of ORFs that exists in the genome or transcriptome.
- candidate_orfs: A list of all candidate ORFs that will be assessed for translation. A subset of ORFs from all_orfs are removed so that the list of ORFs in candidate_orfs do not overlap in the same frame.
- candidate_orfs.gff3: Annotations of all candidate ORFs to enable visualization of the ORF structure in a genome broswer such as IGV.

./iRibo --RunMode=GetCandidateORFs --Genome=path/to/genome.fa --Annotations=path/to/annotations.gtf

Options:
- --Transcriptome=path/to/transcriptome.gtf: If a transcriptome is given, candidate ORFs will be generated from transcriptome sequences rather than directly from genomic sequence.
- --Output=path/to/output_folder: Specify output directory.
- --Threads=1: Set number of threads.

The columns of candidate_orfs and all_orfs are:
- CandidateORF_ID: index of the ORF.
- Transcript_ID: transcript_id of the transcript the ORF is on, as labeled in the transcriptome input file.
- Gene_ID: If the coordinates of the ORF exactly line up with a canonical gene CDS, then the name of that gene is indicated. Otherwise, X.
- contig: Index that corresponds to the contigs or chromosomes in the genome annotation file in order.
- strand: Strand the ORF is on.
- ORF_coord1: First genomic coordinate of the ORF.
- ORF_coord2: Last genomic coordinate of the ORF.
- genomic_coordinates: Comma-separated list of the genomic coordinates of the exons of each ORF. 
- ORF_length: Exonic length of the ORF.
- antisense_gene: If the ORF overlaps an exon of an annotated protein-coding gene on the opposite strand, its name is given here. X otherwise.
- gene_intersect: If the ORF overlaps the CDS of an annotated protein-coding gene on either strand, its name is given here. X otherwise.
- contig_str: Name of the contig/chromosome the ORF is located on.


------------------------------------------------------------------------------

## GenerateTranslationProfile:

Generate a genome-wide translation profile using aligned ribo-seq reads, which will be used in the next step to derive a list of translated ORFs. This step generates several output files:
- translation_calls: Statistics on ribo-seq reads mapping to each ORF.
- null_distribution: Statistics on scrambled reads for each ORF, used to generate a null distribution and calculate the false discovery rate. 
- all_passed_reads_f: Information on plus strand reads that pass quality control. 
- all_passed_reads_r: Information on minus strand reads that pass quality control. 
- riboseq_reads_plus.wig: Tracks of plus strand ribo-seq reads. Can be visualized in a genome browser.
- riboseq_reads_minus.wig: Tracks of minus strand ribo-seq reads. Can be visualized in a genome browser.

To run:
./iRibo --RunMode=GenerateTranslationProfile --Genome=path/to/genome.fa --Riboseq=path/to/sams.txt --CandidateORFs=path/to/candidate_orfs

SAMs File:
sams.txt should list paths to all SAM/BAM files, separated by lines:

e.g., 

sam_dir/SRR1042853_aligned.out.bam

sam_dir/SRR1042855_aligned.out.bam

Options:
- --Output=path/to/output_folder: Define output directory.
- --Threads=1: Set thread count.
- --Min_Length=25: Minimum ribo-seq read length for inclusion in analysis. Most riboseq reads are between length 25-35nt, but can vary by experiment.
- --Max_Length=35: Maximum ribo-seq read length for inclusion in analysis. Most riboseq reads are between length 25-35nt, but can vary by experiment.
- --P_Site_Distance=20: Maximum distance to check for a P-site from the start of the read. 
- --QC_Count=10000: Minimum number of reads mapping to the first codon position of canonical coding sequences for inclusion of any read length from an experiment in the analysis. 
- --QC_Periodicity=2.0: To pass quality control if QC_positions if false, the first codon position among canonical coding sequences must contain more reads than both the second and third positions by at least this factor. To pass quality control if QC_positions is true, there must be more first codon positions with at least one read among canonical coding sequences than second and third positions by at least this factor. 
- --QC_Positions=false: Use positions or read counts in quality control.
- --Scrambles=1: Set number of scrambles to calculate FDR. 

Description of translation_calls output file:
- index: The index of the ORF, as in the candidate_orfs file.
- frame0: The number of codons in the ORF in which the first position of the codon has more mapped ribo-seq reads than the second or third position. This is the number of "successes" in the binomial test for three nucleotide periodicity.
- frame_sum: The number of codons in the ORF in which either the first, second, or third positions in the codon have more reads than the other two positions. This is the number of trials in the binomial test for three nucleotide periodicity.
- reads0, reads1, reads2: The total count of reads mapping to the first, second, and third positions in each codon, respectively.

Description of null_distribution output file:
- index: The index of the ORF, as in the candidate_orfs file
- scrambledN: For each scrambled replicate N, the number of codons in the scrambled ORF in which the first position of the codon has more mapped ribo-seq reads than the second or third position. This is the number of "successes" in the binomial test for three nucleotide periodicity.
- scrambled_sumN: For each scrambled replicate N, the number of codons in the scrambled ORF in which either the first, second, or third positions in the codon have more reads than the other two positions. This is the number of trials in the binomial test for three nucleotide periodicity.

Description of all_passed_reads_f and all_passed_reads_r files:
- chr: Index of the contig or chromosome the read maps to. 
- strand: The strand the read aligns to, 0 for plis, 1 for minus.
- pos: The genomic position of the start of the read.
- count: How many reads map to this contig, strand, and position.

------------------------------------------------------------------------------

## GenerateTranslatome:

This step produces a list of translated ORFs. Running this step creates three output files:
- translated_orfs.csv: Data on all inferred translated ORFs. The same output format as candidate_orfs, but with 3 new columns described below. 
- nORF_discovery.png: A figure showing the number of actual vs. scrambled noncanonical (unannotated) ORFs detected at a range of p-value thresholds, indicating the specified FDR cutoff. 
- cORF_discovery.png: A figure showing the number of actual vs. scrambled canonical ORFs detected at a range of p-value thresholds.
- translated_orfs.gff3: An annotation file giving the bounds of exons of all translated ORFs, which can be visualized in a genome browser like IGV.

To run:
Rscript GenerateTranslatome.R --TranslationCalls=path/to/translation_calls --NullDistribution=path/to/null_distribution --CandidateORFs=path/to/candidate_orfs

Options:
- --Output=path/to/output_folder: Designate output directory.
- --Threads=1: Specify thread count.
- --ExcludeCHR=none: Comma-delimited list of chromosomes/contigs to exclude from analysis. Default is none. Example: --ExcludeCHR=chr1,chr8,chrM
- --ExcludeOverlapGene=True: Exclude noncanonical ORFs overlapping canonical CDSs on the same strand from analysis.
- --FDR=0.05: Define desired false discovery rate.
- --Scrambles=100: Set number of scrambles to calculate FDR. Must be lower than scrambles set in GenerateTranslationProfile step.

Description of translated_orfs.csv:
This file contains the same information as the candidate_orfs file for ORFs inferred to be translated at the specified FDR, but with three additional columns.

- in_frame_reads: number of reads mapping to first codon positions (i.e., in-frame reads) on the ORF  
- Expression-level: ribo-seq read count divided by the length of the ORF in nucleotides
- p_value: the p-value for binomial test for three nucleotide periodicity.
  
------------------------------------------------------------------------------
