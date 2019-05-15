#.libPaths("C:/Users/mgunady/Documents/R/win-library/3.5")

prepareImports <- function() {
	install.packages("BiocManager", repos='http://cran.us.r-project.org')
	library(BiocManager)

	list.of.packages <- c("GenomicFeatures", "Rsamtools", "RCurl")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

	if(length(new.packages)) BiocManager::install(new.packages)
}

#suppressPackageStartupMessages(prepareImports())


####################
## Load Libraries
####################

library(GenomicFeatures)
library(Rsamtools)

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("GenomeInfoDb", "GenomicFeatures", "Rsamtools"))

# Read arguments
args = commandArgs(trailingOnly = TRUE)

# Uncomment the next line and use the right paramters when running from R directly
###################################################################################
#args = c("../test_dm3/Drosophila_melanogaster.BDGP5.70.gtf", "../test_dm3/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa", "../test_dm3/")


print("Load GTF File")
txdb = makeTxDbFromGFF(args[1])	# GTF File
print("Load Fasta File")
genome = FaFile(args[2])	# FASTA File
outDir = args[3]	#output Directory
####################################################
####################################################

# do_tx: Whether to generate transcriptome.fa as well from the genome
preprocess_transcriptome <- function(txdb, genome, outDir, do_tx = FALSE) {
  # Prepare Intermediate data files
  ##################################
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  txtome_filename = file.path(outDir, "transcriptome.fa")
  exs_filename = file.path(outDir, "disjoint_bins.tsv")
  exs2bins_filename = file.path(outDir, "exons2bins.tsv")
  txs2exs_filename = file.path(outDir, "txs2bins.tsv")
  
  # Read TXDB data
  #################
  txs2exons = exonsBy(txdb, by ="tx")
  
  # Construct disjoint exons
  ###########################
  exons = unique(unlist(txs2exons))
  DEx = disjoin(exons)
  hits = as.data.frame(findOverlaps(txs2exons, DEx))
  txs2DExs_ranks = split(hits$subjectHits, f=hits$queryHits)
  txs2DExs = sapply(txs2DExs_ranks, function(exs) paste(exs, collapse=","))
  all_txIDs = names(txs2exons)
  names(txs2DExs) = all_txIDs
  
  # Report Disjoint Exons
  ########################
  exs_table = data.frame(chr = seqnames(DEx), 
                         start = start(DEx), end = end(DEx), strand=strand(DEx), 
                         seq = getSeq(genome, DEx))
  
  write.table(exs_table, file=exs_filename, quote = FALSE, col.names = NA, sep="\t")
  
  # Report Exons-to-Disjoint Exons
  #################################
  exhits = as.data.frame(findOverlaps(exons, DEx))
  exs2DEXs_table = data.frame(bin = exhits$subjectHits, 
                              exon = exhits$queryHits, exonInfo = as.character(exons[exhits$queryHits]))
  
  write.table(exs2DEXs_table, file=exs2bins_filename, quote = FALSE, col.names = NA, sep="\t")
  
  
  # Mapping Transcripts to Disjoint Exons
  ########################################
  if(length(txs2exons) != 0) {
    # Get TX2Gene
    txIDs_geneIDs = select(txdb, keys = all_txIDs, 
                           columns=c("GENEID", "TXCHROM", "TXNAME", "TXSTART", "TXEND", "TXSTRAND"), keytype="TXID")
    # Txs with no genes
    txs_undefgene = is.na(txIDs_geneIDs$GENEID)
    if(any(txs_undefgene)) {
      # UNDEF_GENE# for overlapping TXs
      undef_txs2exs = txs2exons[txs_undefgene]
      undef_txs_info = txIDs_geneIDs[txs_undefgene,]
      reduced_txs = GRanges(Rle(undef_txs_info$TXCHROM), 
                            IRanges(undef_txs_info$TXSTART, undef_txs_info$TXEND),
                            strand = Rle(strand(undef_txs_info$TXSTRAND)))
      # reduced_txs = GRangesList(lapply(undef_txs2exs, function(tx) 
      #   GRanges(seqnames(tx[1]), IRanges(min(start(tx)), max(end(tx))), strand(tx[1]))))
      undef_nums = subjectHits(findOverlaps(reduced_txs, reduce(unlist(reduced_txs))))
      txIDs_geneIDs$GENEID[txs_undefgene] <- paste("UNDEF_GENE", undef_nums, sep="_")
      
    }
	
    # Report Genes X Transcripts X DisjointExons
    txs2DExs_vec = as.vector(txs2DExs[match(txIDs_geneIDs$TXID, names(txs2DExs))])
    txs_table = data.frame(chr = txIDs_geneIDs$TXCHROM, geneID = txIDs_geneIDs$GENEID, 
                           txID = txIDs_geneIDs$TXNAME, strand = txIDs_geneIDs$TXSTRAND, bins = txs2DExs_vec)
    write.table(txs_table, file=txs2exs_filename, quote = FALSE, col.names = NA, sep="\t")
  }
  
  # Report Transcriptome
  #######################
  if(do_tx) {
    # Clear file
    close(file(txtome_filename, open="w"))
    
    print("Writing Transcriptome ...")
    tx_seqs = as.vector(extractTranscriptSeqs(genome, txs2exons))
    lines = paste(">", txIDs_geneIDs$TXNAME, " ", txIDs_geneIDs$TXCHROM, ":", txIDs_geneIDs$GENEID, ":", names(tx_seqs), "\n", 
                tx_seqs, sep="")
    
    f = file(txtome_filename, "w")
    writeLines(lines, f)
    close(f)
  }
}

####################################################
####################################################

#seqs = seqlevels(genome)[seqlevels(genome) %in% seqlevels(txdb)]
#seqs = seqs[1:25] #No Alt

#seqlevels(txdb) <- seqs

print("Finding Disjoint Exonic Bins")
startTime = proc.time()

preprocess_transcriptome(txdb, genome, outDir, do_tx = F)

print(proc.time() - startTime)
print("Done Preprocessing")