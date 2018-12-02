
from lib.utils import *

from collections import defaultdict

def loadSegmentsIndex(segFile, DExons):
    segRanges = defaultdict(lambda : defaultdict(set))
    segTxs = {}
    segLens = {}
    with open(segFile) as sFile:
        for lc, line in enumerate(sFile):
            if lc % 2 != 0: #read headers only
                segLens[segID] = len(line.strip())
                continue
            segID, header, txs, stype = line[1:].strip().split()
            segTxs[segID] = str2Set(txs.split(':')[1])
            seqname, geneID, start, exs, end, strand = header.split(':')
            exons = [DExons[int(ex)] for ex in exs[1:-1].split(',')]
            for i, exon in enumerate(exons):
                # Add seg to exon_st loc
                if i == 0:
                    segRanges[geneID][int(start)].add(segID)
                else:
                    if strand == "+":
                        segRanges[geneID][exon.start].add(segID)
                    else:
                        segRanges[geneID][exon.end].add(segID)  
                if i == len(exons)-1:
                    segRanges[geneID][int(end)].add(segID)
                else:
                    if strand == "+":
                        segRanges[geneID][exon.end].add(segID)
                    else:
                        segRanges[geneID][exon.start].add(segID)
    return (segRanges, segTxs, segLens)


# A grange in granges list is a tuple: (start,end)
def getSegsForRanges(geneID, granges, txs, index, segTxs):
    #print(granges)
    segs = set()
    geneIdx = index[geneID]
    #print(geneIdx)
    for i, grange in enumerate(granges):
        segs |= (geneIdx[grange[0]] & geneIdx[grange[1]])
    strictSegs = set()
    for seg in segs:
        #use intersection for .flex, subset for strict
        if txs.intersection(segTxs[seg]):#segTxs[seg].issubset(txs):
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
            [[locs[0], locs[1]], [locs[1], locs[2]], [locs[2], locs[3]]],   #Inclusion
            [[locs[0], locs[3]]]    #Exclusion
            ], switchTxs, locs[2]-locs[1])       
    elif etype == 'MX': # Mutually Exclusive Exons
        return (geneID, [
            [[locs[0], locs[1]], [locs[1], locs[2]], [locs[2], locs[3]]],   #Inclusion
            [[locs[4], locs[5]], [locs[5], locs[6]], [locs[6], locs[7]]]    #Exclusion
            ], switchTxs, locs[2]-locs[1])
    elif etype == 'A5':  # Alt. 5' splice-site
        return (geneID, [
            [[locs[2], locs[2]+1], [locs[2], locs[0]], [locs[0], locs[1]]],   #Inclusion
            [[locs[2], locs[3]]]    #Exclusion
            ], switchTxs, locs[0]-locs[2])
    elif etype == 'A3':  # Alt. 3' splice-site
        return (geneID, [
            [[locs[0], locs[1]], [locs[1], locs[3]], [locs[3]-1, locs[3]]],   #Inclusion
            [[locs[0], locs[3]]]    #Exclusion
            ], switchTxs, locs[3]-locs[1])
    elif etype == 'RI':   #Retained Intron
        return (geneID, [
            [[locs[1], locs[1]+1], [locs[1]+1, locs[2]-1], [locs[2]-1, locs[2]]
             #[locs[0], locs[1]], [locs[2], locs[3]]
             ],   #Inclusion
            [[locs[1], locs[2]]
             #[locs[0], locs[1]], [locs[2], locs[3]]
             ]    #Exclusion
            ], switchTxs, locs[2]-locs[1]) 
    elif etype == 'AF':   #Alt. First Exon
        return (geneID, [
            [[locs[0], locs[1]], [locs[1], locs[2]]],   #Inclusion
            [[locs[3], locs[4]], [locs[4], locs[5]]]    #Exclusion
            ], switchTxs, 0)
    elif etype == 'AL':   #Alt. Last Exon
        return (geneID, [
            [[locs[3], locs[4]], [locs[4], locs[5]]],   #Inclusion
            [[locs[0], locs[1]], [locs[1], locs[2]]]    #Exclusion
            ], switchTxs, 0)
    else:
        return (geneID, [], switchTxs, 0)

def getCumSegLens(segs, segLens):
    totLen = sum([segLens[segID] for segID in segs])
    return (totLen)

def getSegsForIOEFile(ioeF, outF, segRanges, segTxs, segLens):
    with open(ioeF) as f, open(outF, 'w') as outFile, open(outF+'.flex', 'w') as outFS:
        header = '\t'.join(['seqname', 'geneID', 'eventID', 'incSegs', 'exSegs', 'incTxs', 'exTxs',
                            'incSegLen', 'exSegLen', 'incLen'])+'\n'
        outFile.write(header)
        outFS.write(header)
        f.readline()
        for i, ioeLine in enumerate(f):
            #print(i)
            seqname, geneID, eventID, inc_txs, tot_txs = ioeLine.strip().split('\t')
            geneID, ranges, switchTxs, incLen = event2Ranges(eventID)
            incTxs = str2Set(inc_txs)
            exTxs = str2Set(tot_txs) - incTxs
            if switchTxs:
                buf = incTxs
                incTxs = exTxs
                exTxs = buf
            incSegs, strictIncSegs = getSegsForRanges(geneID, ranges[0], incTxs, segRanges, segTxs)
            exSegs, strictExSegs = getSegsForRanges(geneID, ranges[1], exTxs, segRanges, segTxs)
            exSegs -= incSegs
            
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


def generateEventsSegsIOE(segFile, DExons, eventsFile, outFname):

    segRanges, segTxs, segLens = loadSegmentsIndex(segFile, DExons)
    getSegsForIOEFile(eventsFile, outFname, segRanges, segTxs, segLens)
