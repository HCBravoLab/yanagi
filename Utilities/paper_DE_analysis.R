library("ggplot2")
library(stringr)
library(reshape2)

library(gtable)
library(grid)

library(tidyverse)
library(limma)
library(ROCR)

loadSegmentLib <- function(libPath) {
  lines = readLines(libPath)
  headers = lines[c(TRUE, FALSE)]
  tokens = strsplit(headers," ",fixed=TRUE)
  segIDs = sapply(tokens, function(seg_tokens) substr(seg_tokens[1], 2, nchar(seg_tokens[1])))
  txs = sapply(tokens, function(seg_tokens) strsplit(substr(seg_tokens[3], 5, nchar(seg_tokens[3])), ","))
  names(txs) = segIDs
  return(txs)
}

getTx2Pair <- function(pairs, segs_txs) {
  txs_pair = lapply(pairs, function(pair) {
    segs = unlist(strsplit(pair, " "))
    txs = intersect(segs_txs[[segs[1]]], segs_txs[[segs[2]]])
  })
  names(txs_pair)=pairs
  txs_pair = cbind(melt(txs_pair), 1)
  colnames(txs_pair) = c("tx", "pair", "spliced")
  return(txs_pair)
}

loadTXTPM <- function(sample_ids, sample_genome) {
  TPM_df_list = list()
  for(id in sample_ids) {
    TPM = read.table(paste("simData/", sample_genome, "_sample", id, "_TPM.txt", sep=""), header=TRUE)
    TPM_df_list[[id]] = data.frame(TPM, sample=id)
  }
  TPM_df = do.call(rbind, TPM_df_list)
}

getSegCounts <- function(isample, segs, txs) {
  
  segs_txs_list = mapply(c,as.vector(segs),as.vector(txs),SIMPLIFY=FALSE)
  
  fdf = count_df %>%
    filter(sample==isample) %>%
    droplevels
  
  segCounts = sapply(segs_txs_list, function(seg_tx) {
    v = as.vector(unlist(seg_tx))
    seg = v[1]
    tx = v[2]
    segCount = fdf %>% filter(SEG1ID == seg | SEG2ID == seg) %>% filter(tx %in% strsplit(as.character(TXS), ","))
    if(nrow(segCount) > 0) {
      sum(segCount$count)
    } else {
      0
    }
  })
  return(segCounts)
}

LSV_TXTPM_SEGCounts <- function(txIDs, TXTPM, genome, L) {
  LSV = read.table(paste(genome, "/", genome, "_segs_", L, ".fa.segsWithSkippedExons.tsv", sep=""), header=TRUE)
  
  single_LSV = LSV %>%
    filter(NumSegs==2 & NumTxs==2) %>%
    filter(TxsWithEx %in% txIDs & TxsSkipping %in% txIDs) %>%
    droplevels
  
  eventIDs = paste("event", 1:(dim(single_LSV)[1]), sep="")
  print(length(eventIDs))
  
  data_points_list = list()
  for(id in sample_ids) {
    print(id)
    sample_TXTPM = TXTPM %>%
      filter(sample==id) %>%
      droplevels
    tpms = as.vector(sample_TXTPM$TPM)
    txs_vec = as.vector(sample_TXTPM$transcript_id)
    
    tI_TPM = tpms[match(as.vector(single_LSV$TxsWithEx), txs_vec)]
    tE_TPM = tpms[match(as.vector(single_LSV$TxsSkipping), txs_vec)]
    IE_TXTPM_ratio = log2(tI_TPM/tE_TPM)
    #whichDiff = abs(IE_TXTPM_ratio) >= 0
    #IE_TXTPM_ratio = IE_TXTPM_ratio[whichDiff]
    #diff_single_LSV = single_LSV[whichDiff,]
    
    
    segI_counts = getSegCounts(samples_desc[id], single_LSV$SegWithEx, single_LSV$TxsWithEx)
    segE_counts = getSegCounts(samples_desc[id], single_LSV$SegSkipping, single_LSV$TxsSkipping)
    IE_SEGCounts_ratio = log2(segI_counts/segE_counts)
    
    data_points_list[[id]]= data.frame(IE_TXTPM_ratio = IE_TXTPM_ratio, IE_SEGCounts_ratio = IE_SEGCounts_ratio, sample=id,
                                       segI = single_LSV$SegWithEx, segE = single_LSV$SegSkipping,
                                       txI = single_LSV$TxsWithEx, txE = single_LSV$TxsSkipping,
                                       tI_TPM=tI_TPM, tE_TPM=tE_TPM, 
                                       segI_counts=segI_counts, segE_counts=segE_counts, 
                                       eventIDs=eventIDs)
  }
  data_points = do.call(rbind, data_points_list) %>%
    mutate(label =ifelse(IE_TXTPM_ratio*IE_SEGCounts_ratio >= 0, "T", "F"))
  
  return(data_points)
}


getLMFit <- function(df, varI, varE) {
  matI=acast(df, eventIDs~sample, value.var=varI)#"tI_TPM")#"segI_cpk")#"segI_counts")
  matE=acast(df, eventIDs~sample, value.var=varE)#"tE_TPM")#"segE_cpk")#"segE_counts")
  mat = cbind(matI, matE)
  colnames(mat) = c(paste("I", 1:6, sep=""), paste("E", 1:6, sep=""))
  
  design = matrix(1, nrow=12, ncol=3)
  colnames(design) = c("intercept", "condition", "segType")
  design[,2] = c(c(0, 0, 0, 1, 1, 1), 
                 c(0, 0, 0, -1, -1, -1))
  design[,3] = c(rep(1, 6), rep(0, 6))
  
  exp <- voom(counts = mat, design = design, plot = T)
  fit <- lmFit(exp, design = design)
  fit <- eBayes(fit)
  ranks = topTable(fit, coef = "condition", adjust = 'BH', n = Inf)
  
  labels = df$DE_status[match(rownames(ranks), df$eventIDs)]
  labels[labels=="DE"] = 1
  labels[labels=="NDE"] = 0
  names(labels) = rownames(ranks)
  
  return(list("fit"=fit,
              "ranks"=ranks,
              "labels"=labels))
}


##############################################
##############################################

##############################################
##############################################

L = 10000
genome="hg37"
sample_genome="Hs"

experiment = paste(genome, "_segs", sep="")
samples = rep("", 6)
num_samples = length(samples)
for(i in 1:num_samples) {
  samples[i]=paste(sample_genome, "_sample_", i, sep="")
}
samples_desc = paste("S", 1:num_samples, sep="")

######################
## Load Segments Lib
######################
segs_txs = loadSegmentLib(paste(genome, "/", experiment, "_", L, ".fa", sep = ""))

######################
## Load Segment counts
######################
dir = "seg_counts/"
df = NULL
inter_segs = NULL
df_list = list()
for(i in 1:num_samples) {
  sample = samples[i]
  sample_desc = samples_desc[i]
  print(sample)
  d = read.table(paste(dir, "rapmap_", sample, "__", experiment, "_", L, "_seg_counts.tsv", sep=""), header = TRUE)
  pair = paste(d$SEG1ID, d$SEG2ID)
  df_list[[i]] = data.frame(d, sample=sample_desc, pair=pair)
  if(is.null(inter_segs)) {
    inter_segs = pair
  } else
    inter_segs = intersect(inter_segs, pair)
}

count_df = do.call(rbind, df_list)

count_df = count_df %>%
  mutate(eff_length = ifelse(SEGTYPES %in% c("J", "E"), SEG1LEN-101, SEG1LEN-101+SEG2LEN-101)) %>%
  mutate(cpk = count / eff_length * 10)

sample_ids = 1:6
TXTPM = loadTXTPM(sample_ids, sample_genome)

txIDs= unique(TXTPM$transcript_id)
data_points = LSV_TXTPM_SEGCounts(txIDs, TXTPM, genome, L)

#######
## DE
#######

DE_txIDs = read_lines(paste("simData/", sample_genome, "_DETxs.txt", sep=""))
data_points = data_points %>% mutate(DE_status=ifelse(txI %in% DE_txIDs & txE %in% DE_txIDs, "DE", "NDE"))

SRRTXTPM = read.table(paste("seg_counts/SRR_", sample_genome, "_kallisto_abundance_TXNAMES.tsv", sep=""), header=TRUE)
# SRRTXTPM = SRRTXTPM %>%
#   filter(TXNAME %in% txIDs) %>%
#   droplevels

SRRTI = SRRTXTPM$SRR.tpm[match(data_points$txI[data_points$sample==1], SRRTXTPM$TXNAME)]
SRRTE = SRRTXTPM$SRR.tpm[match(data_points$txE[data_points$sample==1], SRRTXTPM$TXNAME)]

f_data_points = data_points %>%
  mutate(SRRTI_TPM = rep(SRRTI, 6), SRRTE_TPM = rep(SRRTE, 6)) %>%
  filter(SRRTI_TPM+SRRTE_TPM > 1)

eventIDs = unique(f_data_points$eventIDs)
print(length(eventIDs))

segcounts_fit = getLMFit(f_data_points, "segI_counts", "segE_counts")
TXTPM_fit = getLMFit(f_data_points, "tI_TPM", "tE_TPM")

df = data.frame(eventIDs=eventIDs, segs_logFC=segcounts_fit$ranks$logFC[match(eventIDs, rownames(segcounts_fit$ranks))],
                segs_t=segcounts_fit$ranks$t[match(eventIDs, rownames(segcounts_fit$ranks))],
                tpm_logFC=TXTPM_fit$ranks$logFC[match(eventIDs, rownames(TXTPM_fit$ranks))], 
                tpm_t=TXTPM_fit$ranks$t[match(eventIDs, rownames(TXTPM_fit$ranks))],
                DE_status=segcounts_fit$labels[match(eventIDs, names(segcounts_fit$labels))],
                segs_sigma=segcounts_fit$fit$sigma[match(eventIDs, rownames(segcounts_fit$fit$coefficients))],
                tpm_sigma=TXTPM_fit$fit$sigma[match(eventIDs, rownames(TXTPM_fit$fit$coefficients))])


############
## Plots
############


ggplot(df, aes(segs_logFC, tpm_logFC, color=DE_status)) + geom_point() + geom_abline(slope = 1)

ggplot(df, aes(segs_t, tpm_t, color=DE_status)) + geom_point()  + geom_abline(slope = 1)

ggplot(df, aes(segs_sigma, tpm_sigma, color=DE_status)) + geom_point()  + geom_abline(slope = 1)


ggplot(rbind(data.frame(sigma=df$segs_sigma, type="segs"), data.frame(sigma=df$tpm_sigma, type="tpm")), aes(sigma, fill=type)) +
  geom_histogram(position=position_dodge(width=0.1))

lists = list(preds=list("segcounts" = abs(segcounts_fit$ranks$t),
                        "tpm" = abs(TXTPM_fit$ranks$t)),
             labels=list("segcounts" = segcounts_fit$labels,
                         "tpm" = TXTPM_fit$labels))
preds = prediction(lists$preds, lists$labels)
many.roc.perf = performance(preds, measure = "tpr", x.measure = "fpr")
plot(many.roc.perf, col=rainbow(10), colorize=TRUE)
abline(a=0, b= 1)

acc.perf = performance(preds, measure = "acc")
plot(acc.perf)

auc.perf <- performance(preds, measure = "auc")
print(auc.perf@y.values)


plot_df = data.frame(x=many.roc.perf@x.values[[1]], y=many.roc.perf@y.values[[1]], fit="Seg.Counts")
plot_df = rbind(plot_df, data.frame(x=many.roc.perf@x.values[[2]], y=many.roc.perf@y.values[[2]], fit="Tx.TPM"))

ggplot(data=plot_df, aes(x=x, y=y, color=fit)) +
  geom_line() +
  geom_abline(slope=1) +
  ggtitle(paste("(", genome, " genome)", sep="")) +
  labs(y = "True positive rate", x="False positive rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="top", legend.title = element_blank())

ggsave(paste(sample_genome, "_sim_ROC.pdf", sep=""), width=4, height=4)
