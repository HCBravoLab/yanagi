from collections import defaultdict
from collections import Counter

import itertools
import subprocess
import threading

import os, sys, time
import glob

from lib.Segments import *
from lib.utils import *

class alignThread (threading.Thread):
   def __init__(self, threadID, cmd):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.cmd = cmd
      self.aligns = defaultdict(list)
   def run(self):
      print("Starting " + self.threadID)
      self.aligns = runAlignment(self.cmd)

def writeSegPairCounts(outputCountsFilename, segPairs_counts, newsegPairs_counts, segPairs_txs, segsDict):
    with open(outputCountsFilename, "w") as f:
        f.write("seg1ID\tseg2ID\tcount\tgeneID\tseg1Len\tseg2Len\tseg1StLoc\tseg2StLoc\ttxs\n")
        
        for segPair in sorted(segPairs_counts):
            count = segPairs_counts[segPair]
            segs = [segsDict[segID] for segID in segPair.split("_")]
            line = "\t".join([segs[0].ID, segs[1].ID, str(count), 
                              segs[0].geneID,
                              str(segs[0].length), str(segs[1].length),
                              str(segs[0].startLoc), str(segs[1].startLoc), ','.join(segPairs_txs[segPair])])
            f.write(line + "\n")
    with open(outputCountsFilename+".newJuncs", "w") as f:
        f.write("seg1ID\tseg2ID\tcount\tgeneID\tseg1Len\tseg2Len\tseg1StLoc\tseg2StLoc\n")
        
        for segPair in sorted(newsegPairs_counts):
            count = newsegPairs_counts[segPair]
            segs = [segsDict[segID] for segID in segPair.split("_")]
            line = "\t".join([segs[0].ID, segs[1].ID, str(count), 
                              segs[0].geneID, segs[1].geneID,
                              str(segs[0].length), str(segs[1].length),
                              str(segs[0].startLoc), str(segs[1].startLoc)])
            f.write(line + "\n")

def writeSegCounts(outputCountsFilename, seg_counts, segsDict):
    with open(outputCountsFilename, "w") as f:
        f.write("segID\tcount\tgeneID\tsegLen\tsegStLoc\n")
        for segID in sorted(seg_counts):
            seg = segsDict[segID]
            count = seg_counts[segID]
            line = "\t".join([seg.ID, str(count), seg.geneID,
                              str(seg.length), str(seg.startLoc)])
            f.write(line + "\n")
    

######################################
######################################

def readSAM(fin):
    readAligns = defaultdict(list)
    #c = 0
    if sys.version_info[0] == 3: # use for python3
         p = fin
    else:
         p = iter(fin.readline, b'')
    for line in p:
        #print("======================", line)
        tokens = line.strip().split()
        if line.startswith('@'): #If header line:
            continue
        #print(tokens)
        readID, flags, segID, pos = tokens[:4]
        if readID[-2]=='/': # take off the suffixes /1,/2 from readID
            readID = readID[:-2]
        if segID == '*': # end unmapped
            continue
        #readLen = len(tokens[9])
        readAligns[readID].append((segID, int(flags)))
        #c = c+1
        #if c % 1000000 == 0:
        #    print("Seen", str(c), "alignments")
    return(readAligns)

def validOrientation(flags1, flags2):
    flag1Bin = bin(flags1)[2:]
    flag2Bin = bin(flags2)[2:]
    if(len(flag1Bin) < 5):
        flags1Orientation = 0
    else:
        flags1Orientation = 16*int(flag1Bin[-5])
    if(len(flag2Bin) < 5):
        flags2Orientation = 0
    else:
        flags2Orientation = 16*int(flag2Bin[-5])
    validOrientationTrue = ((flags1Orientation + flags2Orientation) == 16)
    return validOrientationTrue
        
def processPairAligns(aligns1, aligns2, segsDict):
    segPairs_counts = Counter()
    newsegPairs_counts = Counter()
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    
    segPairsTxs = {}
    for readID in aligns1:
        a1 = aligns1[readID]
        a2 = aligns2[readID]
        if len(a1) == 0 or len(a2) == 0:    # only one-end mapped
            continue
        mapped = False
        passed_txs = False
        for align in itertools.product(a1, a2):
            segID1, segID2 = align[0][0], align[1][0]
            pairkey = segID1 + "_" + segID2
            if pairkey in segPairsTxs:
                txs = segPairsTxs[pairkey]
            else:
                seg1 = segsDict[segID1]
                seg2 = segsDict[segID2]
                txs = (seg1.txs & seg2.txs)
                segPairsTxs[pairkey] = txs
            validOri = validOrientation(align[0][1], align[1][1])
            if len(txs) < 1 and validOri: # alignment between segs with no tx in common
              newsegPairs_counts[pairkey] += 1    # Counted as a novel junction
            elif not validOri:  # not valid alignment
              passed_txs = True
            else:
              segPairs_counts[pairkey] += 1
              mapped = True
              break
        if mapped:  # Correctly counted (Mapped read)
            readsMapped[0] += 1
        elif passed_txs:    # Passed the txs check, but not valid alignment (Unmapped read)
            readsMapped[2] += 1
        else:   # Didn't find any segments-pair with common txs (Unmapped read)
            readsMapped[1] += 1
    print("Mapped Reads:", readsMapped[0], "Unmapped Reads:", readsMapped[1]+readsMapped[2], \
          "Unmapped Due Txs:", readsMapped[1], "Unmapped Due Direction:", readsMapped[2])
    return(segPairs_counts, newsegPairs_counts, segPairsTxs)

def processSingleAligns(aligns):
    seg_counts = Counter()
    readsMapped_cnt, multimapped_cnts = 0, 0
    for readID in aligns:
        als = aligns[readID]
        if len(als) > 1:
            multimapped_cnts += 1
        for al in als:
            flags = al[1]
            if flags == 0 or flags == 16:
                seg_counts[al[0]] += 1
    print("Mapped Reads:", readsMapped_cnt, "Multimapped Reads:", multimapped_cnts)
    return(seg_counts)

def runAlignment(cmd):
    print("Running Alignment Command...", cmd)
    start_t = time.time()
    align_proc = subprocess.Popen(cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=sys.stdout, universal_newlines=True)
    with align_proc.stdout:
        aligns = readSAM(align_proc.stdout)
    print("Done")
    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)
    return(aligns)


#############################
######## Main Logic #########
#############################

def alignAndCount(segmentReferenceFilename, outf, cmd1, cmd2=None):

    print("Loading Segments Lib...")
    start_t = time.time()
    segsDict, segIDs = load_SegmentsLib(segmentReferenceFilename)
    print("Done!")
    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)
    ##print(process.memory_info().rss)

    # alignment threads
    thread1 = alignThread("AlignThread-1", cmd1.split())
    thread1.start()

    if cmd2:
        thread2 = alignThread("AlignThread-2", cmd2.split())
        thread2.start()

    thread1.join()
    if cmd2:
        thread2.join()

    print("Processing Alignments...")
    start_t = time.time()
    if cmd2:    # Paired-end
        segPairs_counts, newsegPairs_counts, segPairs_txs = processPairAligns(thread1.aligns, thread2.aligns, segsDict)
    else:       # Single-end
        seg_counts = processSingleAligns(thread1.aligns)

    print("Done!")
    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)
    #print(process.memory_info().rss)

    print("Writing Segments Counts...")
    #print(seg_counts)
    start_t = time.time()
    if cmd2:
        writeSegPairCounts(outf, segPairs_counts, newsegPairs_counts, segPairs_txs, segsDict)
    else:
        writeSegCounts(outf, seg_counts, segsDict)

    elapsed = time.time() - start_t
    print("Elapsed Time: ", elapsed)

if __name__ == '__main__':
    createSegments(101, "C:/yanagi_new/output/", "hg37_segs101.fa")

