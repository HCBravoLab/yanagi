
def str2Set(mystr):
    return set(mystr.split(','))

def set2Str(myset):
    return ','.join(sorted(myset))

def txs2segs(segsDict):
    txsPerGene = Counter()
    txs2segs = defaultdict(set)
    for segID in segsDict:
        seg = segsDict[segID]
        for txID in seg.txs:
            if not txID in txs2segs:
                txsPerGene[seg.geneID] += 1
            txs2segs[txID].add(segID)
    return(txs2segs, txsPerGene)

#pos: 1-indexed, returns 1-indexed
def posInTx2Exs(txID, pos, txsDict, exBins_dict):
    tx = txsDict[txID]
    offset = 0
    for ex in tx.exs:
        exB = exBins_dict[ex]
        if offset+exB.width >= pos:
            return(ex, pos-offset)
        offset += exB.width

#pos: 1-indexed
def posInEx2Segs(exID, pos, exs2segs_dict, segs_dict, exBins_dict):
    segIDs = exs2segs_dict[exID]
    resSegIDs = []
    for segID in segIDs:
        seg = segs_dict[segID]
        if exID == seg.exs[0]:
            ex = exBins_dict[exID]
            if seg.startLoc <= ex.start+pos-1:
                resSegIDs.append(segID)
        elif exID == seg.exs[-1]:
            ex = exBins_dict[exID]
            if seg.endLoc >= ex.start+pos-1:
                resSegIDs.append(segID)
        else:
            resSegIDs.append(segID)
    return(resSegIDs)

def calcDistInExs(fromEx, fromPos, toEx, toPos, exs, exBins_dict):
    if fromEx == toEx:
        return(toPos-fromPos+1)
    dist = 0
    fromIdx = exs.index(fromEx)
    dist += exBins_dict[fromEx].width - fromPos + 1
    toIdx = exs.index(toEx)
    dist += toPos
    if toIdx == fromIdx+1:
        return(dist)
    for i in range(fromIdx+1, toIdx):
        dist += exBins_dict[exs[i]].width
    return(dist)

def rangeInTx2Segs(txID, frompos, topos, txsDict, exBins_dict, segs_dict, exs2segs_dict):
    fromEx, i = posInTx2Exs(txID, frompos, txsDict, exBins_dict)
    fromSegs = posInEx2Segs(fromEx, i, exs2segs_dict, segs_dict, exBins_dict)
    
    toEx, j = posInTx2Exs(txID, topos, txsDict, exBins_dict)
    toSegs = posInEx2Segs(toEx, j, exs2segs_dict, segs_dict, exBins_dict)

    interSegs = set(fromSegs).intersection(set(toSegs))
    resSegs = []
    for segID in interSegs:
        if calcDistInExs(fromEx, i, toEx, j, segs_dict[segID].exs, exBins_dict) == topos-frompos+1:
            seg = segs_dict[segID]
            if txID in seg.txs:
                resSegs.append(segID)
    #if txID == "ENST00000341065":
        #print(fromEx, i, fromSegs, toEx, j, toSegs)
    return(resSegs, (fromEx, toEx))

def sortAndMergePairCounts(SC_pairs):
    SCPair_dict = {}
    for segID1, segID2, SC in SC_pairs:
        key = segID1+"_"+segID2 if segID1 <= segID2 else segID2+"_"+segID1
        if key in SCPair_dict:
            SCPair_dict[key] += int(SC)
        else:
            SCPair_dict[key] = int(SC)
    return(SCPair_dict)
