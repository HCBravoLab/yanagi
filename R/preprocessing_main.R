
.libPaths("C:/Users/umiacs/Documents/R/win-library/3.5")

# source("https://bioconductor.org/biocLite.R")

library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)

# Read arguments
args = commandArgs(trailingOnly = TRUE)

print("Load GTF File")
txdb = makeTxDbFromGFF(args[1])	# GTF File
print("Load Fasta File")
genome = FaFile(args[2])	# FASTA File
outDir = args[3]	#output Directory
####################################################
####################################################

seqs = seqlevels(genome)[seqlevels(genome) %in% seqlevels(txdb)]
seqs = seqs[1:25] #No Alt

seqlevels(txdb) <- seqs

print("Finding Disjoint Exonic Bins")
startTime = proc.time()


source(file.path(getwd(), "R", "preprocess_transcriptome.R"))
preprocess_transcriptome(txdb, genome, outDir, do_tx = F)

print(proc.time() - startTime)
print("Done Preprocessing")