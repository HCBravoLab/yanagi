
# source("https://bioconductor.org/biocLite.R")
####################################################
####################################################
# outDir ="hg19\\preprocessed"
# 
# library(BSgenome.Hsapiens.UCSC.hg19)
# genome = BSgenome.Hsapiens.UCSC.hg19
# 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
####################################################
####################################################
# outDir ="mm10\\preprocessed"
# 
# library(BSgenome.Mmusculus.UCSC.mm10)
# genome = BSgenome.Mmusculus.UCSC.mm10
# 
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
####################################################
####################################################
# outDir ="dm6\\preprocessed"
# 
# library(BSgenome.Dmelanogaster.UCSC.dm6)
# genome = BSgenome.Dmelanogaster.UCSC.dm6
# 
# library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
####################################################
####################################################
# outDir ="hg38\\preprocessed"
# 
# library(BSgenome.Hsapiens.UCSC.hg38)
# genome = BSgenome.Hsapiens.UCSC.hg38
# 
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)

####################################################
####################################################
# outDir ="hg37\\preprocessed"
# 
# genome = FaFile("Homo_sapiens.GRCh37.71.dna.primary_assembly.fa") #readDNAStringSet("Homo_sapiens.GRCh37.71.dna.primary_assembly.fa")
# txdb = makeTxDbFromGFF("Homo_sapiens.GRCh37.71.gtf")
####################################################
####################################################
outDir ="dm3\\preprocessed"

genome = FaFile("Drosophila_melanogaster.BDGP5.70.dna.chromosome.2L.fa") #readDNAStringSet("Homo_sapiens.GRCh37.71.dna.primary_assembly.fa")
txdb = makeTxDbFromGFF("Drosophila_melanogaster.BDGP5.70.gtf")
####################################################
####################################################

seqs = seqlevels(genome)[seqlevels(genome) %in% seqlevels(txdb)]

seqlevels(txdb) <- seqs

startTime = proc.time()

source("preprocess_transcriptome.R")
preprocess_transcriptome(txdb, genome, outDir, do_tx = FALSE)

print(proc.time() - startTime)