
from lib.utils import *

from collections import defaultdict

def loadSegmentsIndex(segFile, DExons):
    segRanges = defaultdict(lambda : defaultdict(set))
    segTxs = {}
    segLens = {}
    with open(segFile) as sFile:
        for lc, line in enumerate(sFile):
            if lc == 0: #skip header
                continue
            tokens = line.strip().split("\t")
            segID, chrm, geneID, txAnIDs, binIDs, start, end, strand, length = tokens[:9]
            segTxs[segID] = str2Set(txAnIDs)
            segLens[segID] = int(length)
            exons = [DExons[int(ex)] for ex in binIDs.split(',')]
            if strand == "-":
                if int(start) == exons[0].start:
                    start = str(exons[0].end)
                if int(end) == exons[-1].end:
                    end = str(exons[-1].start)
            for i, exon in enumerate(exons):
                if strand == "+":
                    st = int(start) if i == 0 else exon.start
                    ed = int(end) if i == len(exons)-1 else exon.end
                else:
                    st = int(start) if i == 0 else exon.end
                    ed = int(end) if i == len(exons)-1 else exon.start
                segRanges[geneID][st].add(segID)
                segRanges[geneID][ed].add(segID)

    segRangesSortedKeys = {}
    for geneID in segRanges:
        segRangesSortedKeys[geneID] = sorted(segRanges[geneID].keys())
    return (segRanges, segRangesSortedKeys, segTxs, segLens)


# A grange in granges list is a tuple: (start,end)
def getSegsForRanges(geneID, granges, txs, index, sortedIdxKeys, segTxs, mode):
    #print(granges)
    segs = set()
    geneIdx = index[geneID]
    sortedKeys = sortedIdxKeys[geneID]
    for i, grange in enumerate(granges):
        segs |= (geneIdx[grange[0]] & geneIdx[grange[1]])
        if grange[2] == "E":
                #print(sortedKeys)
            i = sortedKeys.index(grange[0])+1
            while i < len(sortedKeys) and sortedKeys[i] < grange[1]:
                segs |= geneIdx[sortedKeys[i]]
                i = i+1
    
    strictSegs = set()
    for seg in segs:
        #use intersection for .flex, subset for strict
        if mode == "strict":
            if segTxs[seg].subset(txs):
                strictSegs.add(seg)
        else:
            if segTxs[seg].intersection(txs):
                strictSegs.add(seg)
    return (segs, strictSegs)

def event2Ranges(eventID):
    tokens = eventID.strip().replace(';', ':').replace(':-', '^').replace('-', ':').replace('^', ':-').split(':')
    geneID, etype, seqname = tokens[:3]
    strand = tokens[-1]
    ## Event Ranges
    locs = [int(loc) for loc in tokens[3:-1]]
    switchTxs = (etype == "AL" and strand == "+")
    ## Handle strand-dependent eventTypes
    if strand == "-":
        if etype == "A5":
            etype = "A3"
        elif etype == "A3":
            etype = "A5"
        elif etype == "AF":
            etype = "AL"
            switchTxs = True
        elif etype == "AL":
            etype = "AF"
        elif etype == "MX":
            buf = list(locs)
            locs[:4] = buf[4:]
            locs[4:] = buf[:4]
            switchTxs = True
    #print(tokens)
    if etype == 'SE':   #Skipped Exon
        return (geneID, [
            [[locs[0], locs[1], "J"], [locs[1], locs[2], "E"], [locs[2], locs[3], "J"]],   #Inclusion
            [[locs[0], locs[3], "J"]]    #Exclusion
            ], switchTxs, locs[2]-locs[1])       
    elif etype == 'MX': # Mutually Exclusive Exons
        return (geneID, [
            [[locs[0], locs[1], "J"], [locs[1], locs[2], "E"], [locs[2], locs[3], "J"]],   #Inclusion
            [[locs[4], locs[5], "J"], [locs[5], locs[6], "E"], [locs[6], locs[7], "J"]]    #Exclusion
            ], switchTxs, locs[2]-locs[1])
    elif etype == 'A5':  # Alt. 5' splice-site
        return (geneID, [
            [[locs[2], locs[2]+1, "J"], [locs[2], locs[0], "E"], [locs[0], locs[1], "J"]],   #Inclusion
            [[locs[2], locs[3], "J"]]    #Exclusion
            ], switchTxs, locs[0]-locs[2])
    elif etype == 'A3':  # Alt. 3' splice-site
        return (geneID, [
            [[locs[0], locs[1], "J"], [locs[1], locs[3], "E"], [locs[3]-1, locs[3], "J"]],   #Inclusion
            [[locs[0], locs[3], "J"]]    #Exclusion
            ], switchTxs, locs[3]-locs[1])
    elif etype == 'RI':   #Retained Intron
        return (geneID, [
            [[locs[1], locs[1]+1, "J"], [locs[1]+1, locs[2]-1, "E"], [locs[2]-1, locs[2], "J"]
             #[locs[0], locs[1]], [locs[2], locs[3]]
             ],   #Inclusion
            [[locs[1], locs[2], "J"]
             #[locs[0], locs[1]], [locs[2], locs[3]]
             ]    #Exclusion
            ], switchTxs, locs[2]-locs[1]) 
    elif etype == 'AF':   #Alt. First Exon
        return (geneID, [
            [[locs[0], locs[1], "E"], [locs[1], locs[2], "J"]],   #Inclusion
            [[locs[3], locs[4], "E"], [locs[4], locs[5], "J"]]    #Exclusion
            ], switchTxs, 0)
    elif etype == 'AL':   #Alt. Last Exon
        return (geneID, [
            [[locs[3], locs[4], "J"], [locs[4], locs[5], "E"]],   #Inclusion
            [[locs[0], locs[1], "J"], [locs[1], locs[2], "E"]]    #Exclusion
            ], switchTxs, 0)
    else:
        return (geneID, [], switchTxs, 0)

def getCumSegLens(segs, segLens):
    totLen = sum([segLens[segID] for segID in segs])
    return (totLen)

gevent = ""

def getSegsForIOEFile(ioeF, outF, segRanges, sortedIdxKeys, segTxs, segLens, mode):
    if mode != "strict":
        mode = "flex"
    with open(ioeF) as f, open(outF, 'w') as outFile, open(outF+'.'+mode+".evs2segs", 'w') as outFS, open(outF+'.'+mode+'.segs', 'w') as outFSS:
        header = '\t'.join(['seqname', 'geneID', 'eventID', 'incSegs', 'exSegs', 'incTxs', 'exTxs',
                            'incSegLen', 'exSegLen', 'incLen'])+'\n'
        outFile.write(header)
        outFS.write(header)
        outFSS.write(f.readline())        
        for i, ioeLine in enumerate(f):
            #print(i)
            #print(ioeLine)
            seqname, geneID, eventID, inc_txs, tot_txs = ioeLine.strip().split('\t')
            global gevent
            gevent = eventID
            geneID, ranges, switchTxs, incLen = event2Ranges(eventID)
            incTxs = str2Set(inc_txs)
            exTxs = str2Set(tot_txs) - incTxs
            if switchTxs:
                buf = incTxs
                incTxs = exTxs
                exTxs = buf
            try:
                incSegs, strictIncSegs = getSegsForRanges(geneID, ranges[0], incTxs, segRanges, sortedIdxKeys, segTxs, mode)
                exSegs, strictExSegs = getSegsForRanges(geneID, ranges[1], exTxs, segRanges, sortedIdxKeys, segTxs, mode)
                exSegs -= incSegs
            except Exception as e:
                print("Event "+eventID+" can't be mapped to segments. Annotation doesn't match") 
                outFile.write('\t'.join([seqname, geneID, eventID,
                                         set2Str(set()), set2Str(set()),
                                         set2Str(incTxs), set2Str(exTxs),
                                         str(0),
                                         str(0),
                                         str(0)]) + "\n")
                outFS.write('\t'.join([seqname, geneID, eventID,
                                         set2Str(set()), set2Str(set()),
                                         set2Str(incTxs), set2Str(exTxs),
                                         str(0),
                                         str(0),
                                         str(0)]) + "\n")
                outFSS.write('\t'.join([seqname, geneID, eventID,
                                        set2Str(set()), set2Str(set())]) + "\n")
                print(e)
                continue

            if not strictIncSegs or not strictExSegs:
                print(eventID, "has empty segments set in either inclusion or exclusion")
            
            outFile.write('\t'.join([seqname, geneID, eventID,
                                     set2Str(incSegs), set2Str(exSegs),
                                     set2Str(incTxs), set2Str(exTxs),
                                     str(getCumSegLens(incSegs, segLens)),
                                     str(getCumSegLens(exSegs, segLens)),
                                     str(incLen)]) + "\n")
            outFS.write('\t'.join([seqname, geneID, eventID,
                                   set2Str(strictIncSegs), set2Str(strictExSegs),
                                   set2Str(incTxs), set2Str(exTxs),
                                   str(getCumSegLens(strictIncSegs, segLens)),
                                   str(getCumSegLens(strictExSegs, segLens)),
                                   str(incLen)]) + "\n")
            outFSS.write('\t'.join([seqname, geneID, eventID,
                                    set2Str(strictIncSegs), set2Str(strictIncSegs | strictExSegs)]) + "\n")
        print(i)


def generateEventsSegsIOE(segFile, DExons, eventsFile, outFname, mode):
    segRanges, sortedIdxKeys, segTxs, segLens = loadSegmentsIndex(segFile, DExons)
    getSegsForIOEFile(eventsFile, outFname, segRanges, sortedIdxKeys, segTxs, segLens, mode)


#from ReferenceLoader import *
if __name__ == '__main__':
    DExons = load_disjointExons("../../output")
    generateEventsSegsIOE("../../output/hg19_segs100.fa.meta", DExons, "../../output/hg37_5types_noAlt.ioe",
                                  "../../output/hg19_segs100_hg37_5types_noAlt.ioe", "flex")
