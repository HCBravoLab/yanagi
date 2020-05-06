from collections import defaultdict
from collections import Counter
from itertools import combinations
import numpy as np

from lib.utils import *
from lib.Segments import *
from lib.ReferenceLoader import *

from tqdm import tqdm
import time, os

# overriding function in Segments.py
def load_SegmentsLib(seg_file_meta):
    segs_dict = {}
    segIDs = []
    geneSegsDict = defaultdict(list)
    #geneTxsDict = defaultdict(set)
    with open(seg_file_meta) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            seg = Seg(line.strip())
            segs_dict[seg.ID] = seg
            segIDs.append(seg.ID)
            geneSegsDict[seg.geneID].append(seg.ID)
            #geneTxsDict[seg.geneID] |= set(seg.txs)
    return(segs_dict, segIDs, geneSegsDict)#, geneTxsDict)

# overriding function in ReferenceLoader.py
def load_Txs2Exs(inDir, exBins_dict):
    genes = defaultdict(list)
    geneIDSorted = set()
    txsDict = defaultdict(list)
    txIDs = []
    with open(os.path.join(inDir, 'txs2bins.tsv')) as file_hndl:
        for i, line in enumerate(file_hndl):
            if i == 0:   #Skip header
                continue
            key_id, chrome, geneID, txID, strand, exons = line.strip().split('\t')
            tx = TX(key_id, chrome, txID, exons, strand, exBins_dict)
            genes[geneID].append(txID)
            txsDict[txID] = tx
            geneIDSorted.add(chrome+":"+geneID)
            txIDs.append(txID)
    geneIDSorted = [genechr.split(":")[1] for genechr in sorted(geneIDSorted)]
    return(genes, geneIDSorted, txIDs, txsDict)

####################
## Aux Functions
####################

def calcEffLen(seg1ID, seg2ID, txs, segs_dict, DExons):
    seg1, seg2 = segs_dict[seg1ID], segs_dict[seg2ID]
    len1, len2 = seg1.length, seg2.length
    # Find overlap
    overlap = 0
    interDExs = set(seg1.exs) & set(seg2.exs)
    for e in interDExs:
        ex = DExons[e]
        st1, st2 = ex.start, ex.start
        end1, end2 = ex.end, ex.end
        if e == seg1.exs[0]:
            st1 = seg1.startLoc
        if e == seg2.exs[0]:
            st2 = seg2.startLoc
        if e == seg1.exs[-1]:
            end1 = seg1.endLoc
        if e == seg2.exs[-1]:
            end2 = seg2.endLoc
        overlap += max(min(end1, end2) - max(st1, st2), 0)
    if len1 > len2:
        len1 -= overlap
    else:
        len2 -= overlap
    lens_eff = [0]*len(txs)
    gaps = [0]*len(txs)
    # If zero overlap, find gap width per transcript
    if overlap == 0:
        gap = 0
        for i, tx in enumerate(txs):
            pos1 = seg1.posInTx[seg1.txs.index(tx)]
            pos2 = seg2.posInTx[seg2.txs.index(tx)]
            gap = pos2 - (pos1+seg1.length-1)
            if LEN_FR - gap < 2*LEN_R or len1 < LEN_R or len2 < LEN_R or len1 + len2 < LEN_FR - gap:
                len_eff = 1
            else:
                fr_st = max(len1 - (LEN_FR - (LEN_R+gap)), 0)
                fr_ed = min(len1 - LEN_R, len1 - (LEN_FR - (len2+gap)))
                len_eff = fr_ed - fr_st +1
                
            #if LEN_R > LEN_FR - (len2+gap):
            #    len_eff = - LEN_R + (LEN_FR-(LEN_R+gap)) + 1
            #else:
            #    len_eff = len2-LEN_R + 1

            #fr_st = len1 - (LEN_FR - (LEN_R+gap))
            #fr_end = min(len1 - LEN_R, len1 - ())
            #lens_eff[i] = max(1, lens_eff[i] - gap)
            gaps[i] = gap
            lens_eff[i] = len_eff
    # Calc eff length
    else:
        if len1 + len2 < LEN_FR or len1 < LEN_R or len2 < LEN_R:
            len_eff = 1
        else:
            len_eff = min(len1, len2) - LEN_R + 1
        lens_eff = [len_eff]*len(txs)

    #if seg1.geneID == "ENSG00000100097":
    #    print(seg1, seg2, len1, len2, overlap, gaps, lens_eff)
    return(lens_eff, overlap)

class SCVal:

    def __init__(self, segsMClass, count, segs_dict, DExons):
        self.SClassID = ""
        self.superTxs = []
        self.subTxs = []
        self.subSegKeys = []
        self.lens_eff = []
        self.genes = set()
        self.count = 0
        self.numSubTxs = 0

        lens_eff = []
        segPairs = segsMClass.split(',')        
        for segPair in segPairs:
            seg1ID, seg2ID = segPair.split("/")
            gene1, gene2 = segs_dict[seg1ID].geneID, segs_dict[seg2ID].geneID
            txs = list(set(segs_dict[seg1ID].txs) & set(segs_dict[seg2ID].txs))
            #txs = segTxs[seg1ID] & segTxs[seg2ID]
            segPairID = ""
            if seg1ID == seg2ID:
                segPairID = seg1ID
                len_eff = max(segs_dict[seg1ID].length - LEN_FR + 1, 1)
                lens_eff = [len_eff]*len(txs)
            else:
                segPairID = seg1ID + "_" + seg2ID
                lens_eff, overlap = calcEffLen(seg1ID, seg2ID, txs, segs_dict, DExons)
                #elif LEN_R + segs_dict[seg2ID].length < LEN_FR:
                #    fr_maxStart_in_seg1 = LEN_FR - segs_dict[seg2ID].length
                #len_eff = max( - LEN_FR + 1, 1)#min(LEN_FR, segs_dict[seg1ID].length) - LEN_R + 1
            for i, tx in enumerate(txs):
                if tx in self.superTxs:
                    self.lens_eff[self.superTxs.index(tx)] += lens_eff[i]
                else:
                    self.superTxs.append(tx)
                    self.lens_eff.append(lens_eff[i])
            
            self.subSegKeys.append(segPairID)
            self.subTxs.append(txs)
            self.genes.add(gene1)
            self.genes.add(gene2)
            self.numSubTxs += len(txs)
        self.count = count
        self.SClassID = ','.join(self.subSegKeys)


    def __str__(self):
        return "%s:%s:%d -- %s | %s" % (','.join(self.subSegKeys),
                                             ','.join(sorted(self.genes)), self.count,
                                             ','.join(self.superTxs), ','.join([str(l) for l in self.lens_eff]))
    def __repr__(self):
        return self.__str__()

def loadSegCounts(fname, segs_dict, DExons, segCounts,
                  genesClusters, geneClustDict, segCountsClusters,
                  isInputMulticlass=True):
    tot_counts = 0
    geneClustersCount = len(genesClusters)
    mode = 'SE' #SE for single-end, PE for paired-end
    len_eff = 0
    with open(fname) as f:
        secondCol = f.readline().split('\t')[1]
        if secondCol == 'SEG2ID' or isInputMulticlass:
            mode = 'PE'
        for i, line in enumerate(f):
            #if isInputMulticlass and i % 5000 == 0:
            #    print(i)
            tokens = line.strip().split('\t')
            if mode == 'SE':    #TODO Incomplete
                segID = tokens[0]
                seg = segs_dict[segID]
                len_eff = seg.length - LEN_R + 1
                segCounts[segID] = (seg.txs, len_eff, int(tokens[1]))
                tot_counts += int(float(tokens[1]))
            else:
                if isInputMulticlass:
                    #print(tokens[:10])
                    #print(line)
                    segsMClass, count = tokens[:2]
                else:
                    segID1, segID2, count = tokens[:3]
                    segsMClass = segID1+"/"+segID2
                if int(count) < 3 or (int(count) < 3 and isInputMulticlass): continue
                scval = SCVal(segsMClass, int(count), segs_dict, DExons)
                #if not isInputMulticlass and segID1 != segID2:
                #    scval.count = 1
                segPairID = scval.SClassID
                segCounts[segPairID] = scval
                #print(scval)
                [segCountsClusters[g].add(segPairID) for g in scval.genes]

                # Check if count is across genes
                genesClustIdxs = set([geneClustDict[g] for g in scval.genes])
                if len(genesClustIdxs) > 1:
                    newClust = set()
                    for clustIdx in genesClustIdxs:
                        clust = genesClusters[clustIdx]
                        newClust |= clust
                    oldClustIdxs = list(genesClustIdxs)
                    # Replace first clust in the list of clusts with the new clust
                    # and make everything in the new clust refer to it
                    newClustIdx = oldClustIdxs[0]   
                    genesClusters[newClustIdx] = newClust
                    geneClustersCount -= len(genesClustIdxs)-1
                    for oldclustIdx in oldClustIdxs[1:]:
                        for g in genesClusters[oldclustIdx]:
                            geneClustDict[g] = newClustIdx
                        genesClusters[oldclustIdx] = set()
                tot_counts += int(count)
    return(tot_counts, geneClustersCount)

def loadNewSegCounts(fname, segs_dict, genesClusters, geneClustDict):
    tot_counts = 0
    NSCs= defaultdict(list)
    with open(fname) as f:
        for line in f:
            tokens = line.strip().split("\t")
            seg1ID, seg2ID, count, pos1, pos2 = tokens[0:5]
            seg1 = segs_dict[seg1ID]
            seg2 = segs_dict[seg2ID]
            if seg1.geneID == seg2.geneID:
                NSCs[seg1.geneID].append((seg1ID, seg2ID, count, pos1, pos2))

###################################
## Discovery stuff
def findMissingCompatability(genesCluster, geneSegsDict, segs_dict,
                             segCountsKeys, segCounts, txs2SCsMat, HMatrix,
                             txs, thetas, TPMs):
    segIDs = []
    for geneID in genesCluster:
        segIDs.extend(geneSegsDict[geneID])

    badTxs = set()
    # Find segments with one tx
    for segID in segIDs:
        seg = segs_dict[segID]
        if len(seg.txs) == 1:
            coverageSum = 0
            for scKey in segCountsKeys:
                if segID in scKey and "," not in scKey:
                    scval = segCounts[scKey]
                    for i, tx in enumerate(scval.superTxs):
                        #if segID == "SEG0331316":
                        #    print(scKey, tx, scval)
                        coverageSum += (scval.count / scval.lens_eff[i]) if scval.lens_eff[i] > 1 else 0
            sumTPMs = 0
            for tx in seg.txs:
                j = txs.index(tx)
                sumTPMs += TPMs[j]
            print(segID, seg.txs, sumTPMs, coverageSum)
            if sumTPMs - coverageSum > 2:
                badTxs |= set(seg.txs)

    HMatrix_ = np.append(HMatrix, np.zeros((1, HMatrix.shape[1])), axis=0)
    txs2SCsMat_ = np.append(txs2SCsMat, np.zeros((1, txs2SCsMat.shape[1])), axis=0)
    print(HMatrix.shape, HMatrix_.shape, txs2SCsMat.shape, txs2SCsMat_.shape)
    #print(HMatrix, HMatrix_, txs2SCsMat, txs2SCsMat_)
    for tx in sorted(badTxs):
        j = txs.index(tx)
        sckeyIdxs = list(np.nonzero(txs2SCsMat[j, :])[0])
        for sckeyIdx in sckeyIdxs:
            print(tx, j, sckeyIdxs)

    


###################################
## Quant stuff
def initAbundances(txs):
    thetas0 = np.full(shape=len(txs), fill_value=1.0/len(txs))
    return thetas0

def getTxs2SegsMapping(segCountsKeys, segCounts, txs, txsLens):
    tx2SCMatrix = np.zeros(shape=(len(txs), len(segCountsKeys)))
    HMatrix = np.zeros(shape=(len(txs), len(segCountsKeys)))
    tot_counts = 0
    for i, scKey in enumerate(segCountsKeys):
        scval = segCounts[scKey]
        #print(scval)
        tot_counts += scval.count
        for tx in scval.superTxs:
            j = txs.index(tx)
            tx2SCMatrix[j, i] = scval.count
            HMatrix[j, i] = calcH(tx, scKey, scval)
    #print(tx2SCMatrix, HMatrix)
    #x = raw_input("getTxs2SegsMapping - Next: ")
    HMatrix = normalizeHMatrix(HMatrix, tx2SCMatrix, txsLens)
    #print(tx2SCMatrix, HMatrix)
    #x = raw_input("getTxs2SegsMapping - Next: ")
    return(tx2SCMatrix, HMatrix, tot_counts)

def calcH(tx, scKey, scval, model="U"):
    if model == "U":
        return scval.count #scval.lens_eff[scval.superTxs.index(tx)]
    elif model == "NU1":
        return scval.count
    elif model == "NU2":
        return scval.count*1.0/scval.numSubTxs
    return 0

def normalizeHMatrix(HMatrix, tx2SCMatrix, txsLens, model="U"):
    if model == "U":
        HMatrix /= np.maximum(txsLens[:,None] - LEN_FR + 1, 1) #txsLens[:,None]
    elif model == "NU1" or model == "NU2":
        rowSums =  np.sum(HMatrix, axis=1)
        HMatrix /= np.maximum(rowSums, 1)[:,None]
    return(HMatrix)

def runEM(txs, txs2SCsMat, HMatrix, tot_counts):
    ################
    ## EM Algorithm
    ################
    ## Termination condition parameters
    maxIters = 100000
    # Inits
    iters = 0
    maxDiff = 1
    thetas = initAbundances(txs)   # Abundances

    #print(txs2SCsMat, HMatrix)
    # EM loop
    while iters < maxIters and maxDiff > maxDiffThreshold:
        #print(">>>>>", iters, thetas)
        # E-step:
        gamma_ = HMatrix*thetas[:,None]
        colSums = np.sum(gamma_, axis=0)
        #print(gamma_, colSums)
        gamma = gamma_ / colSums
        #print(gamma)
        # M-step:
        new_thetas_ = np.sum(txs2SCsMat*gamma, axis=1)
        new_thetas = new_thetas_ / np.sum(new_thetas_)
        # Check Convergance
        maxDiff = np.max(np.abs(new_thetas - thetas)*tot_counts)
        #print(maxDiff)
        thetas = new_thetas
        iters += 1
        #x = raw_input("Next: ")

    return(thetas, iters)

####################
###################################
###################################

def quantTxsMain(segCountsKeys, segCounts, segs_dict,
                 txs, txsDict, debug=False):

    if not segCountsKeys:
        return (np.zeros(len(txs)), 0)

    txsLens = np.array([txsDict[tx].length for tx in txs])

    start = time.time()
    #print(txsLens)
    txs2SCsMat, HMatrix, tot_counts = getTxs2SegsMapping(segCountsKeys, segCounts,
                                                         txs, txsLens)

    thetas, iters = runEM(txs, txs2SCsMat, HMatrix, tot_counts)

    #input("sgs")
    #abunds = thetas*tot_counts
    #TPMs = abunds/np.maximum(txsLens - LEN_FR + 1.0, 1.0)
    if debug:
        print("EM: ")
        print(txs)
        print(list(thetas*tot_counts))
        #print(list(TPMs))

    end = time.time()
    em_time = end - start

    return(thetas, tot_counts)

###################################
###################################



####################

maxDiffThreshold = 0.001
def quantifyTxSample(genesClusters, geneClustDict, segCounts, segCountsClusters,
                     segs_dict, txsByGenes, txsDict, discover):
    tot = 0
    TPM_denom = 0
    resDict = {}
    for i in tqdm(range(len(genesClusters))):
        cluster = genesClusters[i]
        #if i % 1000==0: print(i)

        if not cluster:
            continue
        segCountsKeys = set()
        txs = set()
        for gene in cluster:
            segCountsKeys |= segCountsClusters[gene]
            txs |= set(txsByGenes[gene])
        sorted_txs = sorted(list(txs))
        segCountsKeys = sorted(list(segCountsKeys))
        #print(segCountsKeys, sorted_txs)
        thetas, gCounts = quantTxsMain(segCountsKeys, segCounts, segs_dict, sorted_txs, txsDict, debug=False)
        tot += gCounts
        
        for txID, theta in zip(sorted_txs, thetas):
            tlen = txsDict[txID].length
            tlen_eff = max(1.0, tlen-LEN_FR+1.0)
            cnt = theta*gCounts
            cnrm = cnt/tlen_eff
            resDict[txID] = [txID, tlen, tlen_eff, cnt, cnrm]
            TPM_denom += cnrm
        
    return(resDict, tot, TPM_denom)

LEN_R = 100
LEN_FR = 250
def quantifyTxExperiment(workDir, samplesFnames, outDir, segsMetaFname,
                         L_fr, L_read, discover=False):
    global LEN_FR, LEN_R
    LEN_FR, LEN_R = L_fr, L_read
    print("Load Segments Library")
    segs_dict, segIDs, segsByGenes = load_SegmentsLib(segsMetaFname)
    DExons = load_disjointExons(workDir)
    txsByGenes, genes, txIDs, txsDict = load_Txs2Exs(workDir, DExons)

    for i, sampleFname in enumerate(samplesFnames):
        sampleID = os.path.basename(sampleFname).split(".")[0]
        print("===============================")
        print("Processing sample#"+str(i+1)+":", sampleID)
        
        ## Init dataStructures
        genesClusters = [set([g]) for g in genes]
        geneClustDict = dict(zip(genes, range(len(genes))))
        segCounts = {}
        segCountsClusters = defaultdict(set)

        ## Load SegCounts
        print("...Load Main SegmentCounts")
        tot_counts = 0
        counts, geneClustersCount = loadSegCounts(sampleFname, segs_dict, DExons, segCounts,
                                                  genesClusters, geneClustDict, segCountsClusters,
                                                  isInputMulticlass=False)
        tot_counts += counts
        #print("...SCs = ", tot_counts)
        #print("...Number of gene Clusters =", geneClustersCount)

        print("...Load Multimapped SegmentCounts")
        counts, geneClustersCount = loadSegCounts(sampleFname+".multimappings", segs_dict, DExons, segCounts,
                                                  genesClusters, geneClustDict, segCountsClusters,
                                                  isInputMulticlass=True)
        tot_counts += counts
        #print("...Total SCs = ", tot_counts)
        print("...Total number of gene Clusters =", geneClustersCount, "from", tot_counts, "reads")

        if discover:
            print("...Load Novel SegmentCounts")
            counts = loadNewSegCounts(sampleFname+".njpcounts", segs_dict, genesClusters, geneClustDict)
            print("...Total Novel SCs = ", counts)

        print("...Estimating Abundances")
        res, totCnts, TPM_denom = quantifyTxSample(genesClusters, geneClustDict, segCounts, segCountsClusters,
                                                   segs_dict, txsByGenes, txsDict, discover)

        print("...Saving Output")
        outPath = os.path.join(outDir, sampleID+".tsv")
        with open(outPath, "w") as fout:
            fout.write('\t'.join(["target_id", "length", "eff_length", "est_counts", "tpm"])+'\n')
            for txID in txIDs:
                ans = res[txID]
                ans[-1] = ans[-1]*1000000.0/TPM_denom
                ans[:3] = [str(x) for x in ans[:3]]
                ans[3:] = [str("%.5f" % x) for x in ans[3:]]
                fout.write('\t'.join(ans)+'\n')


#import glob
if __name__ == '__main__':
    print("Hello There!")
