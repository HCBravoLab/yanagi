
# do_tx: Whether to generate transcriptome.fa as well from the genome
preprocess_transcriptome <- function(txdb, genome, outDir, do_tx = FALSE) {
  # Prepare Intermediate data files
  ##################################
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  txtome_filename = file.path(outDir, "transcriptome.fa")
  exs_filename = file.path(outDir, "disjoint_exons.txt")
  txs2exs_filename = file.path(outDir, "txs_exons.txt")
  
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
  
  write.table(exs_table, file=exs_filename, quote = FALSE, col.names = FALSE)
  
  # Mapping Transcripts to Disjoint Exons
  ########################################
  if(length(txs2exons) != 0) {
    # Get TX2Gene
    txIDs_geneIDs = select(txdb, keys = all_txIDs, 
                           columns=c("GENEID", "TXCHROM", "TXNAME", "TXSTART", "TXEND", "TXSTRAND"), keytype="TXID")
    # Txs with no genes
    txs_undefgene = is.na(txIDs_geneIDs$GENEID)
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
    
    # Report Genes X Transcripts X DisjointExons
    txs2DExs_vec = as.vector(txs2DExs[match(txIDs_geneIDs$TXID, names(txs2DExs))])
    txs_table = data.frame(chr = txIDs_geneIDs$TXCHROM, geneID = txIDs_geneIDs$GENEID, 
                           txID = txIDs_geneIDs$TXNAME, exs = txs2DExs_vec, strand = txIDs_geneIDs$TXSTRAND)
    write.table(txs_table, file=txs2exs_filename, quote = FALSE, col.names = FALSE)
  }
  
  # Report Transcriptome
  #######################
  if(do_tx) {
    # Clear file
    close(file(txtome_filename, open="w"))
    
    print("Writing Transcriptome ...")
    tx_seqs = as.vector(extractTranscriptSeqs(genome, txs2exons))
    lines = paste(">", txIDs_geneIDs$TXCHROM, ":", txIDs_geneIDs$GENEID, ":", names(tx_seqs), "\n", 
                tx_seqs, sep="")
    
    f = file(txtome_filename, "w")
    writeLines(lines, f)
    close(f)
  }
}