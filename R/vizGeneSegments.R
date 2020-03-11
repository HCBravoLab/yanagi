
library(grid)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(data.table)

#######################################
## Load Reference and Segments Info
#######################################

loadSegmentsMETA <- function(filename) {
  
  print("Loading Segments Metadata")
  df = as.data.frame(fread(filename, sep="\t", fill = TRUE, header=T))
  print(head(df))
  
  ##################################
  print("Preparing segs2exs/txs mapping")
  
  exs = strsplit(as.character(df$binIDs),',')
  names(exs) = df$segID
  
  txs = strsplit(as.character(df$txAnnIDs),',')
  names(txs) = df$segID
  #segs2exs = melt(exs)
  #colnames(segs2exs) = c("exonID", "segID")
  
  
  return(list("segs"=df, "segs2exs_list"=exs, "segs2txs_list"=txs))
}

loadDExs <- function(filename) {
  print("Loading text File")
  d = as.data.frame(fread(filename, sep="\t", fill = TRUE, header=T))
  colnames(d)[1] = "exID"
  print(head(d))
  #colnames(d) = c("exID", "chrm", "start", "end", "strand", "seq")
  d$lens = d$end - d$start + 1
  return(d)
}

loadTxs2Dxs <- function(filename) {
  print("Loading text File")
  d = as.data.frame(fread(filename, sep="\t", fill = TRUE, header=T))
  colnames(d)[c(1,4)] = c("txID", "txName")
  print(head(d))
  
  ##################################
  print("Preparing txs2exs mapping")
  exs = strsplit(as.character(d$bins),',')
  names(exs) = d$txName
  
  return(list("df"=d, "txs2exs_list"=exs))
}

########################################
## Load Counts
########################################

loadTXTPMs <- function(dir_prefix, samples, suffix) {
  df = NA
  for(sample in samples) {
    d = as.data.frame(fread(paste0(dir_prefix, sample, suffix), header = TRUE))
    if(is.na(df)) {
      df = d
      colnames(df) = c("txID", paste0("S", sample))
    } else {
      df[paste0("S", sample)] = d$tpm
    }
  }
  return(df)
}

loadSCs <- function(dir_prefix, samples, suffix) {
  df = NA
  for(sample in samples) {
    d = as.data.frame(fread(paste0(dir_prefix, sample, suffix), header = FALSE))
    colnames(d) = c("segID", "count")
    if(is.na(df)) {
      df = d
      colnames(df) = c("segID", paste0("SC", sample))
    } else {
      df[paste0("SC", sample)] = d$count
    }
  }
  if(grepl(":", df$segID[1])) {
    tokens = do.call('rbind', strsplit(as.character(df$segID),':',fixed=TRUE))
    df$geneID = tokens[,1]
    df$segID = tokens[,2]
  }
  return(df)
}

########################################

plotGene <- function(gene, sampleTable=NA,
                     segsFASTA, DExs, txs2DExs,
                     SCs=NA, txTPMs=NA, DET_true=c(),
                     segs=NA, segsIdxsRange=NA,
                     export_filename=NA,
                     landscape=FALSE,
                     L = 100,
                     view_SC='count',
                     DESegs = c(),
                     axis_size.y=15, axis_size.x=7) {
  # Filter for gene
  segsf = segsFASTA$segs[segsFASTA$segs$geneID == gene,]
  orderedSegs = segsf$segID[order(segsf$st, decreasing = segsf$strand[1]=="+")]
  
  ## Filter Segments
  if(!is.na(segs)) {
    segsf = segsf[segsf$segID %in% segs,]
  } else if(!is.na(segsIdxsRange)) {
    segsf = segsf[segsf$segID %in% orderedSegs[segsIdxsRange],]
  }
  orderedSegs = segsf$segID[order(segsf$st, decreasing = segsf$strand[1]=="+")]
  
  ## Obtain segs2exs table
  segs2exs_listf = segsFASTA$segs2exs_list[segsf$segID]
  segs2exs = melt(table(melt(segs2exs_listf)))
  colnames(segs2exs) = c("exons", "segs", "value")
  orderedExs = as.character(sort(unique(segs2exs$exons), decreasing = segsf$strand[1]!="+"))
  
  segsf$segID = factor(segsf$segID, levels=orderedSegs)
  segs2exs$segs = factor(segs2exs$segs, levels = orderedSegs)
  segs2exs$exons = factor(segs2exs$exons, levels = orderedExs)
  oddRows = orderedExs[c(TRUE, FALSE)]
  #segs2exs$value[segs2exs$exons %in% oddRows & segs2exs$value == 0] = 2
  
  ## Filter Txs
  txsdff = txs2DExs$df[txs2DExs$df$geneID == gene,]
  txsf = txsdff$txName
  txs2exs_listf = txs2DExs$txs2exs_list[as.vector(txsf)]
  
  #tab = table(melt(txs2exs_listf))
  #tab = tab[order(as.numeric(rownames(tab))), ]
  #txs2exs = melt(tab)
  txs2exs = melt(table(melt(txs2exs_listf)))
  
  colnames(txs2exs) = c("exs", "txs", "value")
  orderedTxs = sort(unique(txs2exs$txs))
  txs2exs$exs = factor(txs2exs$exs, levels = orderedExs)
  oddRows = orderedExs[c(TRUE, FALSE)]
  #txs2exs$value[txs2exs$exs %in% oddRows & txs2exs$value == 0] = 2
  
  ## Filter DExs
  DExsf = DExs[orderedExs, ]
  DExsf$exID = factor(DExsf$exID, levels = orderedExs)
  
  ## Filter counts
  if(!is.na(SCs)) {
    SCsf = SCs[match(orderedSegs, SCs$segID),]
    if(is.na(sampleTable)) {
      nsamples = dim(SCs)[2]-1
      sampleTable = data.frame(
        row.names = paste("SC", nsamples, sep=""),
        condition = rep("control", nsamples)
      )
    }
  }
  if(!is.na(txTPMs)) {
    txTPMsf = txTPMs[match(orderedTxs, txTPMs$txID),]
    if(is.na(sampleTable)) {
      nsamples = dim(txTPMs)[2]-1
      sampleTable = data.frame(
          row.names = paste("SC", nsamples, sep=""),
          condition = rep("control", nsamples)
      )
    }
  }
  
  ## Obtain segs2txs table
  segs2txs_listf = segsFASTA$segs2txs_list[as.vector(segsf$segID)]
  segs2txs = melt(table(melt(segs2txs_listf)))
  colnames(segs2txs) = c("txs", "segs", "value")
  segs2txs$segs = factor(segs2txs$segs, levels = orderedSegs)
  oddRows = orderedTxs[c(TRUE, FALSE)]
  # if(!is.na(txTPMs)) {
  #   segs2txs$TPMs = txTPMsf$txID segs2txs$txs
  # } else {
  #   
  # }
  #segs2txs$value[segs2txs$txs %in% oddRows & segs2txs$value == 0] = 2
  segs2txsNum = colSums(table(melt(segs2txs_listf)))
  
  #######################################
  ## Plots
  #######################################
  

  
  ####
  for(segID in names(segs2exs_listf)) { # iterate over the list of segments to be viewed
    seg = segsf[segsf$segID == segID, ] #segment info
    exs = segs2exs_listf[[segID]] #list of exonic bins in that segment
    if(length(exs)==1) { #if segment spans only one exonic bin
      DEx = DExsf[DExsf$exID == exs[1], ] # get that bin
      cover1 = (DEx$end - as.numeric(seg$st)+1)/DEx$lens #portion of the start offset of the segment in that bin
      cover2 = (as.numeric(seg$end) - DEx$start+1)/DEx$lens #portion of the end offset 
      segs2exs$value[segs2exs$exons == DEx$exID &  segs2exs$segs == segID] = ifelse(cover1==1 & cover2==1,1,0.5) #set value to 1 if full, 0.5 if partial
    } else { #if segment spans more than one bin
      DEx = DExsf[DExsf$exID == exs[1], ] # get first bin
      cover = (DEx$end - as.numeric(seg$st)+1)/DEx$lens
      segs2exs$value[segs2exs$exons == DEx$exID &  segs2exs$segs == segID] = ifelse(cover==1,1,0.5)
      DEx = DExsf[DExsf$exID == exs[length(exs)], ] # get last bin
      cover = (as.numeric(seg$end) - DEx$start+1)/DEx$lens
      segs2exs$value[segs2exs$exons == DEx$exID &  segs2exs$segs == segID] = ifelse(cover==1,1,0.5)
    }
    
  }
  ####
  
  exs2segs_p = segs2exs %>%
    mutate(value=ifelse(value==0,NA,value)) %>%
    mutate(shortID=factor(paste0("S.",str_sub(segs, -4)), levels=paste0("S.",str_sub(orderedSegs, -4)))) %>%
    ggplot(aes(x=shortID, y=as.factor(exons))) + geom_tile(aes(fill = as.factor(value), width=0.8, height=0.8)) +
    scale_fill_manual(values=c("salmon", "red", "grey75")) +
    scale_y_discrete(expand = c(0, 0)) +
    #scale_x_discrete(expand = c(0, 0)) +
    theme_minimal()  +
    guides(fill=FALSE) +
    theme(axis.text.x = element_blank()#text(colour = "black", size = 7)
          , axis.text.y = element_text(size=axis_size.y)
          , axis.title.x = element_blank()
          , axis.title.y = element_blank()
          , panel.grid.major = element_line(size=1, colour = "grey90")
    ) + coord_flip()
  
  if(length(DESegs)>0) {
    ids = which(orderedSegs %in% DESegs)
    exs2segs_p = exs2segs_p + geom_vline(xintercept = ids, color="blue")
  }
  
  ######################
  
  txs2exs_p = txs2exs %>%
    mutate(value=ifelse(value==0,NA,value)) %>%
    mutate(shortTxID=factor(paste0("T.",str_sub(txs, -4)), levels=paste0("T.",str_sub(orderedTxs, -4)))) %>%
    ggplot(aes(x=shortTxID, y=exs)) + geom_tile(aes(fill = as.factor(value), width=0.8, height=0.8)) +
    scale_fill_manual(values=c("darkgreen", "grey75")) +
    scale_y_discrete(expand = c(0, 0)) +
    #scale_x_discrete(xpand = c(0, 0)) +
    theme_minimal()  +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(colour = "black",  size = axis_size.x)
          , axis.text.y = element_text(size=axis_size.y)
          ,axis.title.x = element_blank()
          , axis.title.y = element_blank()
          , panel.grid.major = element_line(size=1, colour = "grey90")
          , panel.grid.major.y = element_line(size=1, colour = "grey70")
    ) + coord_flip()
  
  ######################
  
  segLen_p = segsf %>%
    mutate(loglen = log10(length)) %>%
    ggplot(aes(x=segID, y=loglen)) + 
    geom_bar(aes(x=segID, y=loglen), stat = "identity", fill = 'bisque') + 
    geom_abline(intercept = -log10(L), slope = 0, linetype="dashed", colour="red", size = 1) +
    geom_text(aes(x=segID, y=loglen, label=length, hjust= 1.5), position = position_dodge(width=1), size = 7, colour="white", fontface="bold") +
    scale_x_discrete(expand = c(0, 0)) +
    #scale_x_discrete(breaks = orderedSegs, expand = c(0, 0)) +
    #scale_y_log10() +
    theme_minimal() +
    #scale_y_reverse( expand = c(0, 0)) +
    theme(axis.title.x = element_blank()
          , axis.title.y = element_blank()
          , axis.text.x = element_blank()
          , axis.text.y = element_blank()
          , panel.grid = element_blank()
    ) + coord_flip()

  ######################  
  
  txs2segs_p = segs2txs %>%
    mutate(shortTxID=factor(paste0("T.",str_sub(txs, -4)), levels=paste0("T.",str_sub(orderedTxs, -4)))) %>%
    ggplot(aes(x=segs, y=shortTxID)) + geom_tile(aes(fill = as.factor(value), width=0.8, height=0.8)) +
    scale_fill_manual(values=c("grey90", "darkgreen", "grey75")) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal()  +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle = 90, size=9)
          , axis.text.y = element_text()
          , axis.title.x = element_blank()
          , axis.title.y = element_blank()
    ) + coord_flip()
  
  ######################
  
  txslens = melt(sapply(txs2exs_listf, function(exs) sum(DExsf$lens[DExsf$exID %in% exs])))
  colnames(txslens) = c("len")
  
  txLen_p = txslens %>%
    mutate(loglen = log10(len)) %>%
    ggplot(aes(x=rownames(txslens), y=loglen)) + 
    geom_bar(aes(x=rownames(txslens), y=loglen), stat = "identity", fill = 'bisque') + 
    geom_abline(intercept = -log10(L), slope = 0, linetype="dashed", colour="red", size = 1) +
    geom_text(aes(x=rownames(txslens), y=loglen, label=len, hjust= 1.5), position = position_dodge(width=1), size = 7, colour="white", fontface="bold") +
    scale_x_discrete(expand = c(0, 0)) +
    #scale_y_log10() +
    theme_minimal() +
    #scale_y_reverse( expand = c(0, 0)) +
    theme(axis.title.x = element_blank()
          , axis.title.y = element_blank()
          , axis.text.x = element_text(size=axis_size.x)
          , axis.text.y = element_blank()
          , panel.grid = element_blank()
    ) + coord_flip()
  
  #########################
  
  exLen_p = DExsf %>%
    mutate(loglen = log10(lens+1)) %>%
    mutate(exID = exID) %>%
    ggplot(aes(x=exID, y=loglen)) + 
    geom_bar(stat = "identity", fill = 'bisque') + 
    geom_abline(intercept = -log10(L), slope = 0, linetype="dashed", colour="red", size = 1) +
    geom_text(aes(label=lens, hjust= -.5), position = position_dodge(width=1), angle = 90, size = 7, colour="white", fontface="bold") +
    #scale_x_discrete(breaks = orderedExs, expand=c(0,0), position="top") +
    #scale_y_log10() +
    theme_minimal() +
    #coord_flip() +
    scale_y_reverse(  labels=function(l){paste0("        ", l)}) +
    theme(axis.title.y = element_blank()
          , axis.title.x = element_blank()
          , axis.text.y = element_text(size=axis_size.y)
          #, axis.text.x.top = element_blank()#text(angle = 90, vjust = .5)
          , axis.text.x = element_blank()#text(size=7)
          , panel.grid = element_blank()
    )
  
  ############################
  
  if(!is.na(SCs)) {
    sc_p_df = melt(SCsf, variable.name = "sample", value.name = "count") %>%
      mutate(condition = rep(sampleTable$condition, each=length(orderedSegs))) %>%
      mutate(rpk=count*L/(segsf$length)) %>% #*segs2txsNum
      #mutate(marks=markedSegs[segID]) %>%
      mutate(logcount=log2(count+1)) %>%
      mutate(shortID=factor(paste0("S.",str_sub(segID, -4)), levels=paste0("S.",str_sub(orderedSegs, -4))))
    print(view_SC)
    if(view_SC=='count') {
      sc_p_df$sc_y=sc_p_df$logcount
      print('count')
    } else {
      sc_p_df$sc_y=sc_p_df$rpk
      print('rpk')
    }
    sc_p = sc_p_df %>%
      mutate(DE=ifelse(segID %in% DESegs, "DE", "NDE")) %>%
      ggplot(aes(x=shortID, y=sc_y, color=condition)) +
      geom_point(position = position_dodge(width = 0.5), aes(shape=DE), size=2) +
      scale_shape_manual(values=c(16, 1))+
      #geom_text(aes(label=Var2),hjust=0, vjust=0, size=3, position = position_dodge(width=0.5)) +
      theme_minimal() +
      #scale_y_log10() +
      theme(axis.text.x = element_text(size=axis_size.x)
            , axis.text.y = element_blank()
            , axis.title.x = element_blank()
            , axis.title.y = element_blank(),#element_text(angle = 90),
            legend.position="right", legend.title = element_blank()
            #, panel.grid.major.x = element_line(size=1, colour = "grey75")
      ) + coord_flip()
    if(levels(sc_p_df$condition) == 1 & levels(sc_p_df$DE) == 1) {
      sc_p = sc_p + theme(legend.position = "none")
    }
  } else if(!is.na(txTPMs)) {
    sc_p = segLen_p
  } else {
    sc_p = ggplot()
  }
  
  ############################
  
  if(!is.na(txTPMs)) {
    txtpm_p_df = melt(txTPMsf, variable.name = "sample", value.name = "tpm") %>%
      mutate(condition = rep(sampleTable$condition, each=length(orderedTxs))) %>%
      mutate(logTPM=log2(tpm+1)) %>%
      mutate(DE=ifelse(txID %in% DET_true, "DE", "NDE")) %>%
      mutate(shortID=factor(paste0("T.",str_sub(txID, -4)), levels=paste0("T.",str_sub(orderedTxs, -4))))
    
    txtpm_p = txtpm_p_df %>%
      ggplot(aes(x=shortID, y=logTPM, color=condition)) +
      geom_point(position = position_dodge(width=0.5), aes(shape=DE), size=2)+#, size=2) +
      scale_shape_manual(values=c(16, 1))+
      #geom_text(aes(label=sampleID),hjust=0, vjust=0, size=3, position = position_dodge(width=0.5)) +
      #scale_x_discrete(breaks = paste0("T.",str_sub(orderedTxs, -4)), expand = c(0, 0)) +
      theme_minimal() +
      #scale_y_log10() +
      theme(axis.text.x = element_text(size=axis_size.x)
            , axis.text.y = element_blank()
            , axis.title.x = element_blank()
            , axis.title.y = element_blank(),#element_text(angle = 90),
            legend.position="right", legend.title = element_blank()
            #, panel.grid.major.x = element_line(size=1, colour = "grey75")
      ) + coord_flip()
    if(length(levels(txtpm_p_df$condition)) == 1 & length(unique(txtpm_p_df$DE)) == 1) {
      txtpm_p = txtpm_p + theme(legend.position = "none")
    }
  } else if(!is.na(SCs)) {
    txtpm_p = txLen_p
  } else {
    txtpm_p = ggplot()
  }

  
  ############################
  
  ###################################
  ## Plots layout
  ###################################
  
  if(!landscape) {
    
    pg1 <- ggplot_gtable(ggplot_build(ggplot()))
    pg2 <- ggplot_gtable(ggplot_build(txs2exs_p))
    pg3 <- ggplot_gtable(ggplot_build(txtpm_p))
    pg4 <- ggplot_gtable(ggplot_build(ggplot()))
    pg5 <- ggplot_gtable(ggplot_build(exs2segs_p))
    pg6 <- ggplot_gtable(ggplot_build(sc_p))
    pg7 <- ggplot_gtable(ggplot_build(ggplot()))
    pg8 <- ggplot_gtable(ggplot_build(exLen_p))
    pg9 <- ggplot_gtable(ggplot_build(ggplot()))
    
  } else {
    
  }
  
  maxWidth = unit.pmax(pg1$widths[2:3], pg4$widths[2:3], pg7$widths[2:3])
  pg1$widths[2:3] <- maxWidth
  pg4$widths[2:3] <- maxWidth
  pg7$widths[2:3] <- maxWidth
  maxWidth = unit.pmax(pg2$widths[2:3], pg5$widths[2:3], pg8$widths[2:3])
  pg2$widths[2:3] <- maxWidth
  pg5$widths[2:3] <- maxWidth
  pg8$widths[2:3] <- maxWidth
  maxWidth = unit.pmax(pg3$widths[2:3], pg6$widths[2:3], pg9$widths[2:3])
  pg3$widths[2:3] <- maxWidth
  pg6$widths[2:3] <- maxWidth
  pg9$widths[2:3] <- maxWidth
  
  maxHeight = unit.pmax(pg1$heights[2:3], pg2$heights[2:3], pg3$heights[2:3])
  pg1$heights[2:3] <- maxHeight
  pg2$heights[2:3] <- maxHeight
  pg3$heights[2:3] <- maxHeight
  maxHeight = unit.pmax(pg4$heights[2:3], pg5$heights[2:3], pg6$heights[2:3])
  pg4$heights[2:3] <- maxHeight
  pg5$heights[2:3] <- maxHeight
  pg6$heights[2:3] <- maxHeight
  maxHeight = unit.pmax(pg7$heights[2:3], pg8$heights[2:3], pg9$heights[2:3])
  pg7$heights[2:3] <- maxHeight
  pg8$heights[2:3] <- maxHeight
  pg9$heights[2:3] <- maxHeight
  
  if(!is.na(txTPMs) | !is.na(SCs)) {
    lay <- rbind(c(1,1,1,2),
                 c(3,3,3,4),
                 c(3,3,3,4),
                 c(3,3,3,4),
                 c(5,5,5,NA)
    )
    
    p = grid.arrange(pg2, pg3, 
                     pg5, pg6, 
                     pg8, layout_matrix = lay) 
  } else {
      
    lay <- rbind(c(1),
                 c(2),
                 c(2),
                 c(2),
                 c(3)
    )
    
    p = grid.arrange(pg2, 
                     pg5, 
                     pg8, layout_matrix = lay) 

  }
  

  if(!is.na(export_filename)) {
    print(paste0("Exporting Output to ", export_filename))
    ggsave(export_filename, plot = p, device=png(), width=40, height=25)
    dev.off()
  }
  #ggsave(plot = p, paste0(dir_prefix, gene, ".png"), width=20, height=10)
}


################################
################################

# Read arguments
args = commandArgs(trailingOnly = TRUE)
geneID = args[1]
segsMETA = args[2]
segsDir = args[3]
output_fname = args[4]

segsMeta = loadSegmentsMETA(segsMETA)
DExs = loadDExs(file.path(segsDir, "disjoint_bins.tsv"))
txs2DExs = loadTxs2Dxs(file.path(segsDir, "txs2bins.tsv"))




#samples = c("SRR3332174")
#txTPMs = loadTXTPMs("D:/yanagi/fruitfly/kallisto/", samples, ".tpm")
#SCs = loadTXTPMs("D:/yanagi/fruitfly/Dm6_segs100/psi/", samples, ".SPM")
#colnames(SCs)[1] = "segID"
#sampleTable = data.frame(
#  row.names = paste("SC", samples, sep=""),
#  condition = c(rep("control", 1), rep("case", 0))
#)
#gene = "FBgn0015609"

plotGene(geneID, L = 100,
         segsFASTA = segsMeta, DExs = DExs, txs2DExs = txs2DExs
		 , export_filename = output_fname
         #, txTPMs = txTPMs, sampleTable = sampleTable, SCs = SCs
)

