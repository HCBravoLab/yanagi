
from collections import defaultdict
from collections import Counter

from lib.utils import *
from lib.Segments import *

import os

def loadSegCounts(inFile, segsDict):
    segCounts = defaultdict(list)
    tot_counts = 0
    mode = 'SE' #SE for single-end, PE for paired-end
    with open(inFile) as f:
        secondCol = f.readline().split('\t')[1]
        if secondCol == 'SEG2ID':
            mode = 'PE'
            #segCounts = defaultdict(Counter)
        for line in f:
            tokens = line.strip().split('\t')
            if mode == 'SE':    #TODO Incomplete
                segCounts[tokens[0]].append(([], int(float(tokens[1]))))
                tot_counts += int(float(tokens[1]))
            else:
                seg1ID, seg2ID, count = tokens[:3]
                #txs = tokens[-1]
                txs = segsDict[seg1ID].txs & segsDict[seg2ID].txs
                txs = set2Str(txs)
                #txs = segTxs[seg1ID] & segTxs[seg2ID]
                if seg1ID == seg2ID:
                    segCounts[seg1ID].append((txs, 2*int(count)))
                else:
                    segCounts[seg1ID].append((txs, int(count)))
                    segCounts[seg2ID].append((txs, int(count)))
                tot_counts += 2*int(count)
    return(segCounts, tot_counts)

def writeSPM(outFile, segCounts, segIDs, segsDict):
    sc_n = dict([(segID, 0) for segID in segIDs])
    tot_sc_n = 0
    for segID in sorted(segCounts):
        sc = sum([c[1] for c in segCounts[segID]])
        c_n = sc/segsDict[segID].length*1.0
        sc_n[segID] = c_n
        tot_sc_n += c_n
    with open(outFile, 'w') as outF:
        outF.write('\t'.join(["target_id","tpm"])+"\n")
        for segID in sorted(sc_n):
            SPM = sc_n[segID] / tot_sc_n * 1e6
            #print(segID, SPM)
            outF.write('\t'.join([segID,str(SPM)])+"\n")

def sumSegsCount(segs, segCounts, txsStr):
    evTxs = str2Set(txsStr)
    total = 0
    for seg in segs:
        if seg == '':
            #print(txsStr)
            continue
        for (txs, count) in segCounts[seg]:
            if not txs or str2Set(txs).intersection(evTxs):
                total += count
    return total

def quantEvents(ioeSegFile, outFile, segCountsIdx):
    with open(ioeSegFile) as ioeFile, open(outFile, "w") as outF:
        ioeFile.readline()
        outF.write('\t'.join(['eventID', 'incCount', 'exCount', 'PSI',
                              'incSegs', 'exSegs', 'incTxs', 'exTxs',
                              'incSegsLen', 'exSegsLen', 'incLen'])+"\n")
        for line in ioeFile:
            seqname, geneID, eventID, incSegs, exSegs, incTxs, exTxs, incSegsLen, exSegsLen, incLen = line.strip().split('\t')
            incCount = sumSegsCount(str2Set(incSegs), segCountsIdx, incTxs)
            exCount = sumSegsCount(str2Set(exSegs), segCountsIdx, exTxs)
            incNC = incCount/max(int(incSegsLen), 1)
            exNC = exCount/max(int(exSegsLen), 1)
            PSI = incNC/(incNC+exNC) if incNC > 0 else 0
            PSI = int(PSI*100000.0)/100000.0
           
            outF.write('\t'.join([eventID, str(incCount), str(exCount), str(PSI),
                                  incSegs, exSegs, incTxs, exTxs,
                                  incSegsLen, exSegsLen, incLen])+"\n")

def quantifyEvents(segEventsFname, samplesFnames, outDir, out_prefix, segsMetaFname):
    print("Load Segments Library")
    segsDict, segIDs = load_SegmentsLib(segsMetaFname)

    print("Quantify Events")
    
    for sampleFname in samplesFnames:
        sampleID = os.path.basename(sampleFname).split(".")[0]
        outFile = os.path.join(outDir, out_prefix+"_"+sampleID+".psi")
        print(outFile)

        segCountsIdx, tot_counts = loadSegCounts(sampleFname, segsDict)
        print("segCounts Loaded")
        quantEvents(segEventsFname, outFile, segCountsIdx)
        print("Calculate SPM")
        writeSPM(os.path.join(outDir, sampleID+".SPM"), segCountsIdx, segIDs, segsDict)

#import glob
if __name__ == '__main__':
    samplesFnames = glob.glob(os.path.join("../../../yanagi/simUnannTxs/segs/c1/", "*.tsv"))
    segIDs = []
    segLens = {}
    with open("../../output/hg19_segs100_dropped.fa.meta") as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            tokens = line.strip().split("\t")
            segIDs.append(tokens[0])
            segLens[tokens[0]] = int(tokens[-1])
    quantifyEvents("../../output/hg19_segs100_dropped_hg37_7types_dropped.ioe.flex",
                   samplesFnames, "../../../yanagi/simUnannTxs/segs/c1/psi/", "hg37_7types", segIDs, segLens)
