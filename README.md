# Yanagi: Transcript Segment Library Construction for RNA-Seq Quantification
## Source code based on the work presented in Yanagi: Fast and interpretable segment-based alternative splicing and gene expression analysis (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2947-6)
## Update Nov 18th, 2019: Some major changes are pushed to improve the usability of the pipeline. And introduced the use of Yanagi-count as the preferred alignment tool based on RapMap's quasi mapping to provide segment counts.

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
* Python:
  - tqdm

* R (Bioconductor):
  - GenomicFeatures
  - Biostrings
 
----------------------------
**Download**
============================
----------------------------

Download yanagi by cloning the repository through the ``Clone or download`` button on the top right of this page. Or by running the clone command in git. Then change directory into the created directory where yanagi source is downloaded.
```
git clone https://github.com/HCBravoLab/yanagi.git
cd yanagi
```

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
- **psiCalc**       : Calculates PSI values of alternative splicing events based on their segment mappings.

**Note: This tutorial assumes that all commands are excuted from inside the directory where yanagi is downloaded (refer to the previous Download section).**

----------------------------
**Annotation Preprocessing**
============================
----------------------------

Exons (and retained introns) in the transcriptome annotation can be overlapping within a gene (e.g. in 3'/5' splicing) or across genes. In order for Yanagi to guaranteeing L-disjointness property of the generated segments, a preprocessing step is needed to generate disjoint exonic bins.
Yanagi generate disjoint exonic bins and their transcripts mappings from an input annotation file (GTF format) and the genome sequence file (FASTA format).

### **Command and options** ###

To preprocess the transcriptome annotation subject to segmentation one has to run the following command in the following format:
```
python yanagi.py preprocess -gtf <gtf-file> -fa <fasta-file> -o <work-directory>
```
Note that throughout this tutorial, we will use the same directory <output-directory> as the working directory when needed in different commands.

### **Output files** ###

The preprocess operation outputs two main files:

1.  disjoint_bins.tsv: A file with the structural and sequence information of each constructed disjoint exonic bin.

```
	chr	start	end	strand	seq
1	1	11869	11871	+	GTT
2	1	11872	11873	+	AA
3	1	11874	12009	+	CTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCT...
4	1	12010	12057	+	GTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAG
...
```

2.  txs2bins.tsv: A file with transcripts-to-bins information.

```
	chr	geneID	txID	bins	strand
1	1	ENSG00000223972	ENST00000456328	1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18	+
2	1	ENSG00000223972	ENST00000515242	2,3,4,5,6,8,9,12,13,14,15,16,17,18,19	+
3	1	ENSG00000223972	ENST00000518655	3,4,5,6,7,8,9,14,15,17,18	+
...
````

Another output file ``exons2bins.tsv`` is generated from that step. That extra file contains a mapping between the exons/introns annotated in the .gtf file and the disjoint exonic bins (reported in disjoint_bins.tsv file) that are used as the building blocks for the splice graph used inside of yanagi.

**Alternative Splicing Events Generation**
==========================================

If the downstream analysis involves studying alternative splicing events present in the transcriptome. Then this step is needed to prepared the annotation of those events (Skip this step otherwise).
Yanagi uses the same events definition and code used in [SUPPA](https://github.com/comprna/SUPPA)(eventGenerator command).

To generate the list of events given the GTF (unzipped) of the transcriptome one can run that command:
```
python eventGenerator.py -i <gtf-file> -o <output-directory-and-prefix> -f ioe -e <list-of-event-types-space-separated>
```
List of options available:
- **-e**  | **--event-type**: (only used for local AS events) space separated list of events to generate from the following list:

  - **SE**: Skipping exon (SE)
  - **SS**: Alternative 5' (A5) or 3' (A3) splice sites (generates both)
  - **MX**: Mutually Exclusive (MX) exons
  - **RI**: Retained intron (RI)
  - **FL**: Alternative First (AF) and Last (AL) exons (generates both)

Note that a description of each event type and definition can be found on [SUPPA](https://github.com/comprna/SUPPA)'s page.
The command generates a separate .ioe file of the list of events of each event type provided in the event-type option.
The shell script ```merge_ioe_files.sh``` can be edited for use to merge the separate .ioe files into one file, or to filter out events outside of the primary transcriptome assembly.

<a id="-segment"/>

-----------------------
**Segments Generation**
=======================
-----------------------

This command executes the main operation preparing the segments library by Yanagi, to be used later for RNA-seq reads alignment.
Yanagi takes the preprocessed transcriptome as input to build segments graph, which is then parsed to generate minimal L-disjoint segments.

![Yanagi's Segments Example](https://github.com/mgunady/yanagi/blob/master/yanagi_example.png)

**Fig 1.** The figure shows an illustrative example of transcriptome segmentation of one gene with three transcripts. The example shows the final segments generated by yanagi and how reads are aligned to them.

### **Command and options** ###

To segment the transcriptome one has to run the following command in the following format:
```
python yanagi.py segment -l <read-length> -wd <work-directory>
```
List of options available:

- **-l**  | **--max-overlap**: This (integer) parameter value controls the maximum overlap between any two (L-disjoint) generated segments. A typical choice of l would equal to the expected read length. Refer to Yanagi's [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2947-6) for more details on L-disjointness.

- **-wd** | **--work-dir**: This is the work directory where the preprocessed annotation files exist (same output directory used in the preprocess subcommand). This directory must have two files ```disjoint_bins.tsv``` and ```txs2bins.tsv```.

- **-o**  | **--output-name**: (**Optional**) This is a name prefix used to name output files. If not provided, the default output files are named in the format ```segs_<L>```.

- **-ioe**| **--events-annotation**: (**Optional**) This is a list of .ioe files annotating alternative splicing events present in the corresponding transcriptome. Used if downstream analysis is needed on alternative splicing events. More details in <a href="#-psicalc">`Segment-Based PSI Calculation section`</a>. 

### **Output files** ###

The segmentation operation outputs three files:

1. ``<output-name>.fa``: A FASTA file of the segments library representing the transcriptome.

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

2. ``<output-name>.fa.meta``: A file of metadata describing the structure of each segment and how it was formed.

```
segID	chrom	geneID	txAnnIDs	binIDs	st	end	strand
SEG0000001	10	ENSG00000012779	ENST00000542434	57010,57011	45869661	45869774	+
SEG0000002	10	ENSG00000012779	ENST00000374391,ENST00000542434	57011,57012,57013,57014,57015,57016,57017	45869675	45920450	+
SEG0000003	10	ENSG00000012779	ENST00000374391,ENST00000483623,ENST00000542434	57016,57017	45919539	45920580	+
SEG0000004	10	ENSG00000012779	ENST00000483623	57017,57018	45920481	45923934	+
...
```

3. ``<output-name>.gtf``: A GTF file describing both exonic bins and segments. This file is intended for visualization of segments. Check Section Visualization.

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

P.S. Section <a href="#-segsDownloads">`Ready-to-Download Segments Libraries`</a> provides pre-prepared segment libraries of some commonly used genomes.

-------------------------
**Alignment To Segments**
=========================
-------------------------

Since the update on Nov 2019, we strongly recommend to use [Yanagi-count](https://github.com/mgunady/Yanagi-count) as an alignment tool with segments support for better usability and support of yanagi's segments. 
Otherwise, use the following subcommand in yanagi to use any other transcriptome aligner as long as it reports output in SAM format.
-------------------------

This command runs a given alignment command to align RNA-seq reads into the segments library as a reference. The main use of this command is to facilitate aligning paired-end reads to obtain segment-pair counts.

In this tutorial, we assume the use of RapMap (https://github.com/COMBINE-lab/RapMap) to perform pseudo-alignment. However, other alignment tools can be used and the following commands are to be adjusted accordingly.
Note! The indexing step is done separately from yanagi's pipeline. As an example, to index the segments library using RapMap, follow this command format:

```
PATH/TO/RAPMAP quasiindex -t PATH/TO/segments.fa -i quasiindex/output/directory
```

Once the aligner's index is ready, one can run the alignment step using Yanagi's following command:

```
python yanagi.py align -ref <segments-meta> -o <output-filename> -cmd1 'PATH/TO/RAPMAP quasimap -i quasiindex/output/directory -r PATH/TO/FIRST_READS.fa' -cmd2 'PATH/TO/RAPMAP quasimap -i quasiindex/output/directory -r PATH/TO/SECOND_READS.fa'
```
**Alternatively**, use the provided shell script ```run_segAlign.sh``` by first editing the variables used in it, and choosing the right command for single or paired-end modes.

List of options available:

- **-ref**  | **--segs-ref-file**:  Specifies the segments reference metadata file (.fa.meta file generated by segment subcommand) to align reads against. Note that its corresponding FASTA file has to be indexed by the used aligner a priori.

- **-o**  | **--output-name**: This is a name of the output segment count file.

- **-cmd1**| **--align-command1**: This is the command used to run Rapmap's quasi mapping for the reads FASTA file <PATH/TO/FIRST_READS.fa> using segments indexed at <quasiindex/output/directory>.  

- **-cmd2**| **--align-command2**: (**Optional**) For the second-end reads (if paired-end reads). 

### **Output files** ###

The segmentation operation outputs three files:

<output-filename>.txt: A text (TSV file) contains segments counts (or segment-pairs counts if paired-end).

- Segments counts output example (Single-end reads):
```
segID	count	geneID	segLen	segStLoc
SEG0232653	7	ENSG00000000457	510	169822815
SEG0232655	11	ENSG00000000457	1039	169824007
SEG0232667	3	ENSG00000000460	1105	169631245
SEG0232671	12	ENSG00000000460	166	169764190
```

- Segment-Pair counts output example (Paired-end reads):
```
seg1ID	seg2ID	count	geneID	seg1Len	seg2Len	seg1StLoc	seg2StLoc	txs
SEG0232653	SEG0232653	4	ENSG00000000457	510	510	169822815	169822815	ENST00000367772
SEG0232653	SEG0232655	1	ENSG00000000457	510	1039	169822815	169824007	ENST00000367772
SEG0232655	SEG0232653	2	ENSG00000000457	1039	510	169824007	169822815	ENST00000367772
SEG0232655	SEG0232655	8	ENSG00000000457	1039	1039	169824007	169824007	ENST00000367772,ENST00000367771,ENST00000367770
```
Note: txs field for segment-pairs counts represent the intersecton of transcripts from both segments.

-----------------
**Visualization**
=================
-----------------

![Yanagi's Visualization Example](https://github.com/mgunady/yanagi/blob/master/geneViz.png)

To visualize segments of a specific gene, run the R script found in ```R/vizGeneSegments.R``` using Rstudio or the following Rscript command:
```
Rscript R/vizGeneSegments.R <geneID> <segments.fa.meta> <segmentsDir> <output-filename>
```

Support for visualizing segment counts will be added soon.


<a id="-psicalc"/>

---------------------------------
**Segment-based PSI Calculation**
=================================
---------------------------------

After samples are aligned to the segments using command ```align```, one can process the segments/segment-pairs counts obtained to perform alternative splicing analysis. Yanagi provides a command to calculat PSI values based on segments counts in each of the aligned samples.
This command calculates PSI values of alternative splicing events based on their segment mappings. 

### **Command and options** ###

To calculate PSI values one has to run the following command in the following format:
```
python yanagi.py psiCalc  -es <events-to=segments-mapping> -s <segments-meta> -i <segment-counts-directory> -o <output-directory> -opf <output-prefix>
```
List of options available:

- **-es**  | **--events2segs**: This is the path to the file mapping splicing events to segments. This file was prepared from yanagi command ```segment```. By passing the optional option ```-ioe``` with the events annotation file as input into the ```segment``` command, it outputs the mapping into ```.evs2segs``` file under the output directory specified in option ```-o```. Refer to section <a href="#-segment">`Segments Generation section`</a> for how to run ```segment``` command.

- **-s** | **--segs-meta**: This is the segments metadata file ```.fa.meta``` obtained from yanagi command ```segment```. Refer to section <a href="#-segment">`Segments Generation section`</a> for how to run ```segment``` command.

- **-i** | **---segCounts-dir**: The path to a directory with segments counts (or segment-pairs counts) obtained from yanagi command ```align```. The directory will have a .tsv file per sample.

- **-o** | **--out-dir**: The output directory.

- **-opf**  | **--output-prefix**: (**Optional**) This is a name prefix used to name output files. If not provided, the default output files will use the input counts filenames for each sample.

### **Output files** ###

The PSI calculation operation outputs a .psi file per sample. Each .psi file is of the following format:

```
eventID	incCount	exCount	PSI	incSegs	exSegs	incTxs	exTxs	incSegsLen	exSegsLen	incLen
ENSG00000177697;SE:11:836442-836769:836843-837250:+	6.150537634408602	0.6363636363636364	0.9061	SEG0009700,SEG0009707	SEG0009701	ENST00000322008,ENST00000397420,ENST00000397421,ENST00000524748,ENST00000526693,ENST00000527341,ENST00000528011,ENST00000530320,ENST00000530726	ENST00000529810	372	198	74
ENSG00000177697;SE:11:833026-834530:834591-836063:+	0.6874154262516915	0.41695501730103807	0.62188	SEG0009666,SEG0009671,SEG0009673,SEG0009684,SEG0009685,SEG0009686	SEG0009667,SEG0009668,SEG0009678,SEG0009679,SEG0009680	ENST00000322008,ENST00000531999	ENST00000397421,ENST00000526661,ENST00000529810,ENST00000530155	739	578	61
ENSG00000214063;SE:11:842915-847201:847300-850288:+	0.6484149855907781	0.4393939393939394	0.59552	SEG0010382,SEG0010384,SEG0010395,SEG0010396,SEG0010397	SEG0010383	ENST00000397397	ENST00000397411	694	198	99
...
```

<a id="-segsDownloads"/>

---------------------------------------
**Ready-to-download Segment Libraries**
=======================================
---------------------------------------

* Segments for human transcriptome Ensembl GRCh37:
	* For L=100: (<a href="https://zenodo.org/record/2646964/files/hg19_segs100.fa">`FASTA File`</a>, <a href="https://zenodo.org/record/2646964/files/hg19_segs100.fa.meta">`Meta File`</a>)
	* <a href="https://zenodo.org/record/2646964/files/hg19_5types_noAlt.ioe">`AS Events Annotation (SE, MX, RI, A3, A5)`</a>
	* <a href="https://zenodo.org/record/2646964/files/hg19_segs100_hg19_5types_noAlt.ioe.flex">`Events-to-segments mapping (SE, MX, RI, A3, A5)`</a>
	
* Segments for fruit fly transcriptome Ensembl BDGP6:
	* For L=100: (<a href="https://zenodo.org/record/2646964/files/Dm6_segs_100.fa">`FASTA File`</a>, <a href="https://zenodo.org/record/2646964/files/Dm6_segs_100.fa.meta">`Meta File`</a>)
	* <a href="https://zenodo.org/record/2646964/files/Dm6_5types_noAlt.ioe">`AS Events Annotation (SE, MX, RI, A3, A5)`</a>
	* <a href="https://zenodo.org/record/2646964/files/Dm6_segs_100_Dm6_5types_noAlt.ioe.flex">`Events-to-segments mapping (SE, MX, RI, A3, A5)`</a>


---------------------------------------
**License**
=======================================
---------------------------------------

Yanagi is released under the MIT license.
