
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
