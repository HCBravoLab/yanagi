# Yanagi: Transcript Segment Library Construction for RNA-Seq Quantification
## Source code based on the work presented in paper to appear in proceedings of WABI 2017

## Abstract

Analysis of differential alternative splicing from RNA-seq data is complicated by the fact that many RNA-seq reads map to multiple transcripts, and that annotated transcripts from a given gene are often a small subset of many possible complete transcripts for that gene. Here we describe Yanagi, a tool which segments a transcriptome into disjoint regions to create a segments library from a complete transcriptome annotation that preserves all of its consecutive regions of a given length L while distinguishing annotated alternative splicing events in the transcriptome. In this paper, we formalize this concept of transcriptome segmentation and propose an efficient algorithm for generating segment libraries based on a length parameter dependent on specific RNA-Seq library construction. The resulting segment sequences can be used with pseudo-alignment tools to quantify expression at the segment level. We characterize the segment libraries for the reference transcriptomes of Drosophila melanogaster and Homo sapiens. Finally, we demonstrate the utility of quantification using a segment library based on an analysis of differential exon skipping in Drosophila melanogaster and Homo sapiens. The notion of transcript segmentation as introduced here and implemented in Yanagi will open the door for the application of lightweight, ultra-fast pseudo-alignment algorithms in a wide variety of analyses of transcription variation.

## Usage

- Preprocessing:
  - Input: Transcriptome annotation (GTF file) and genome sequence (FASTA file).
  - Output: The processed transcriptome annotation based on disjoint genomic regions. The output includes two files 'disjoint_exons.txt' and 'txs_exons.txt' under the directory $GENOME/preprocessed/
  - How: Edit preprocessing_main.R to set the used genome. Then run it using R/Rscript
  
- Segmentation:
  - Input: The preprocessed transcriptome (both files obtained from the previous step) and the desired parameter value L
  - Output: The segments library as a FASTA file
  - How: Edit yanagi_main.py to set the used genome and value of L. Then run it using python (numpy package required) 
  
  ## Utilities
  
  The 'Utilities' directory contains some scripts used in the analysis presented in the paper.
  Most importantly, file 'segPairQuant.py' is used to quantify the segments pair, for the case of paired-end RNA-Seq reads.

![Alt text](workflow.pdf?raw=true "Yanagi-based Workflow")
