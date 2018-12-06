# Yanagi: Transcript Segment Library Construction for RNA-Seq Quantification
## Update July 10th, 2018: Check our recent manuscript for detailed formulation of segments and their usage in gene and alternative splicing analysis https://www.biorxiv.org/content/early/2018/07/08/364281
## Source code based on the work presented in paper to appear in proceedings of WABI 2017

## Abstract

Analysis of differential alternative splicing from RNA-seq data is complicated by the fact that many RNA-seq reads map to multiple transcripts, and that annotated transcripts from a given gene are often a small subset of many possible complete transcripts for that gene. Here we describe Yanagi, a tool which segments a transcriptome into disjoint regions to create a segments library from a complete transcriptome annotation that preserves all of its consecutive regions of a given length L while distinguishing annotated alternative splicing events in the transcriptome. In this paper, we formalize this concept of transcriptome segmentation and propose an efficient algorithm for generating segment libraries based on a length parameter dependent on specific RNA-Seq library construction. The resulting segment sequences can be used with pseudo-alignment tools to quantify expression at the segment level. We characterize the segment libraries for the reference transcriptomes of Drosophila melanogaster and Homo sapiens. Finally, we demonstrate the utility of quantification using a segment library based on an analysis of differential exon skipping in Drosophila melanogaster and Homo sapiens. The notion of transcript segmentation as introduced here and implemented in Yanagi will open the door for the application of lightweight, ultra-fast pseudo-alignment algorithms in a wide variety of analyses of transcription variation.

## Usage
----------------------------
**SET UP**
==========
----------------------------

**Requirements**
----------------

Yanagi has been developed and tested in Python 3.7 and R 3.5.
Yanagi uses the following modules:
Python:
- tqdm
R (Bioconductor):
- GenomicFeatures
- Biostrings

**Command and subcommand structure**
==============

Yanagi works with a command/subcommand structure:

```
yanagi.py subcommand options
```
where the subcommand can be one of these options:

- **preprocess**    : Preprocesses transcriptome annotation by breaking exons into disjoint exonic bins and find their transcript mapping.
- **segment**       : Generates a set of maximal L-disjoint segments from the preprocessed transcriptome annotation.
- **align**         : Pseudo aligns reads (single or paired-end) into the segments and obtain segment counts (single segment or segment pair counts).

----------------------------
**Annotation Preprocessing**
============================
----------------------------

Exons (and retained introns) in the transcriptome annotation can be overlapping within a gene (e.g. in 3'/5' splicing) or across genes. In order for Yanagi to guaranteeing L-disjointness property of the generated segments, a preprocessing step is needed to generate disjoint exonic bins.
Yanagi generate disjoint exonic bins and their transcripts mappings from an input annotation file (GTF format) and the genome sequence file (FASTA format).

### **Command and options** ###

To preprocess the transcriptome annotation subject to segmentation one has to run the following command in the following format:
```
python yanagi.py preprocess -gtf <gtf-file> -fa <fasta-file> -o <output-directory>
```
Note that throughout this tutorial, we will use the same directory <output-directory> as the working directory when needed in different commands.

### **Output files** ###

The preprocess operation outputs two files:

1.  disjoint_bins.tsv: A file with the structural and sequence information of each constructed disjoint exonic bin.

```
	chr	start	end	strand	seq
1	1	11869	11871	+	GTT
2	1	11872	11873	+	AA
3	1	11874	12009	+	CTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCT...
4	1	12010	12057	+	GTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAG
...
```

2.  txs2bin.tsv: A file with transcripts-to-bins information.

```
	chr	geneID	txID	bins	strand
1	1	ENSG00000223972	ENST00000456328	1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18	+
2	1	ENSG00000223972	ENST00000515242	2,3,4,5,6,8,9,12,13,14,15,16,17,18,19	+
3	1	ENSG00000223972	ENST00000518655	3,4,5,6,7,8,9,14,15,17,18	+
...
````


-----------------------
**Segments Generation**
=======================
-----------------------

This command executes the main operation preparing the segments library by Yanagi, to be used later for RNA-seq reads alignment.
Yanagi takes the preprocessed transcriptome as input to build segments graph, which is then parsed to generate minimal L-disjoint segments.

![alt text](https://github.com/mgunady/yanagi/yanagi_example.png)

**Fig 1.** The figure shows an illustrative example of transcriptome segmentation of one gene with three transcripts. The example shows the final segments generated by yanagi and how reads are aligned to them.

### **Command and options** ###

To segment the transcriptome one has to run the following command in the following format:
```
python yanagi.py segment -l <read-length> -wd <work-directory>
```
List of options available:

- **-l**  | **--max-overlap**: This (integer) parameter value controls the maximum overlap between any two (L-disjoint) generated segments. Yanagi allows for up to L-length overlap between segments to span possible junctions. A typical choice of L would equal to the expected read length.

- **-wd** | **--work-dir**: This is the work directory where the preprocessed annotation files exist (same output directory used in the preprocess subcommand). This directory must have two files ```disjoint_bins.tsv``` and ```txs2bins.tsv```.

- **-o**  | **--output-name**: (**Optional**) This is a name prefix used to name output files. If not provided, the default output files are named in the format ```segs_<L>```.

- **-ioe**| **--events-annotation**: (**Optional**) This is a list of .ioe files created by SUPPA. Used if downstream analysis is needed on alternative splicing events. More details in Segment-Based PSI Calculation section. 

### **Output files** ###

The segmentation operation outputs three files:

1. output_prefix.fa: A FASTA file of the segments library representing the transcriptome.

```
>SEG0000001
GCTAGATGCGGACACCTGGACCGCCGCGCCGAGGCTCCCGGCGCTCGCTGCTCCCGCGGCCCGCGCCATGCCCTCCT...
>SEG0000002
CCTGGACCGCCGCGCCGAGGCTCCCGGCGCTCGCTGCTCCCGCGGCCCGCGCCATGCCCTCCTACACGGT...
>SEG0000003
GGAATGACTTCGCCGACTTTGAGAAAATCTTTGTCAAGATCAGCAACACTATTTCTGAGCGGGTCATGAATCACTG...
>SEG0000004
GATCCGGCGCTGCACAGAGCTGCCCGAGAAGCTCCCGGTGACCACGGAGATGGTAGAGTGCAGCCTGGAG...
...
````

2. output_prefix.meta: A file of metadata describing the structure of each segment and how it was formed.

```
segID	chrom	geneID	txAnnIDs	binIDs	st	end	strand
SEG0000001	10	ENSG00000012779	ENST00000542434	57010,57011	45869661	45869774	+
SEG0000002	10	ENSG00000012779	ENST00000374391,ENST00000542434	57011,57012,57013,57014,57015,57016,57017	45869675	45920450	+
SEG0000003	10	ENSG00000012779	ENST00000374391,ENST00000483623,ENST00000542434	57016,57017	45919539	45920580	+
SEG0000004	10	ENSG00000012779	ENST00000483623	57017,57018	45920481	45923934	+
...
```

3. output_prefix.gtf: A GTF file describing both exonic bins and segments. This file is intended for visualization of segments. Check Section Visualization.

This GTF contains entries of two possible feature type (column 3): ```exonic_bin``` or ```segment```. Each exonic bin or segment has only one entry in the file. Entries that describe lists (like transcripts for exonic bins, or bins for segments) are placed as lists separated by '+'.

- Exonic_bin feature examples:
```
1	hg19_segs101	exonic_bin	11869	11871	.	+	.	gene_id "ENSG00000223972"; entry_id "1"; transcripts "ENST00000456328";
1	hg19_segs101	exonic_bin	11872	11873	.	+	.	gene_id "ENSG00000223972"; entry_id "2"; transcripts "ENST00000456328+ENST00000515242";
1	hg19_segs101	exonic_bin	11874	12009	.	+	.	gene_id "ENSG00000223972"; entry_id "3"; transcripts "ENST00000456328+ENST00000515242+ENST00000518655";
1	hg19_segs101	exonic_bin	12010	12057	.	+	.	gene_id "ENSG00000223972"; entry_id "4"; transcripts "ENST00000450305+ENST00000456328+ENST00000515242+ENST00000518655";
```

- Segment feature examples:
```
10	hg19_segs101	segment	45869661	45869774	.	+	.	gene_id "ENSG00000012779"; entry_id "SEG0000001"; exonic_bin_ids "57010+57011"; transcripts "ENST00000542434";
10	hg19_segs101	segment	45869675	45920450	.	+	.	gene_id "ENSG00000012779"; entry_id "SEG0000002"; exonic_bin_ids "57011+57012+57013+57014+57015+57016+57017"; transcripts "ENST00000374391+ENST00000542434";
10	hg19_segs101	segment	45919539	45920580	.	+	.	gene_id "ENSG00000012779"; entry_id "SEG0000003"; exonic_bin_ids "57016+57017"; transcripts "ENST00000374391+ENST00000483623+ENST00000542434";
10	hg19_segs101	segment	45920481	45923934	.	+	.	gene_id "ENSG00000012779"; entry_id "SEG0000004"; exonic_bin_ids "57017+57018"; transcripts "ENST00000483623";
```


-------------------------
**Alignment To Segments**
=========================
-------------------------

This command runs a given alignment command to align RNA-seq reads into the segments library as a reference. The main use of this command is to facilitate aligning paired-end reads to obtain segment-pair counts.

** Further details are soon to be provided.


-----------------
**Visualization**
=================
-----------------

Details to be added soon


---------------------------------
**Segment-based PSI Calculation**
=================================
---------------------------------

Details to be added soon


