
import lib.utils
from collections import defaultdict

class Seg:
    ID = ""
    header = ""
    txs = []

    chrom = ""
    length = 0
    geneID = ""
    startLoc = 0
    endLoc = 0
    exs = []
    strand = "+"

    posInTx = []

    def __init__(self, header):
        self.header = header

        tokens = header.strip().split("\t")
        segID, chrom, geneID, txIDs, binIDs, st, end, strand, length = tokens[:9]
        self.ID = segID
        self.chrom = chrom
        self.geneID = geneID
        self.startLoc = int(st)
        self.endLoc = int(end)
        self.strand = strand
        self.length = int(length)
        
        self.txs = txIDs.split(",")
        self.exs = [int(ex) for ex in binIDs.split(",")]

        self.posInTx = [int(p) for p in tokens[9].split(",")]

    def __str__(self):
        return "%s:%d" % (self.ID, self.length)
    def __repr__(self):
        return self.__str__()


def load_SegmentsLib(seg_file_meta):
    segs_dict = {}
    segIDs = []
    with open(seg_file_meta) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            seg = Seg(line.strip())
            segs_dict[seg.ID] = seg
            segIDs.append(seg.ID)
    return(segs_dict, segIDs)

def load_segmentsSeq(seg_file):
    segs_dict = {}
    with open(seg_file) as f:
        for i, line in enumerate(f):
            if i % 2 == 0:
                segID = line.strip()[1:]
            else:
                segs_dict[segID] = line.strip()
    return(segs_dict)

def build_Exs2Segs(segs_dict):
    print("Building Segs2Exs index")
    exs2segs = defaultdict(list)
    for segID in segs_dict:
        seg = segs_dict[segID]
        for ex in seg.exs:
            exs2segs[ex].append(segID)
    return(exs2segs)

def calcGCContent(segsSeq_dict):
    segsGC = {}
    for segID in segsSeq_dict:
        seq = segsSeq_dict[segID].upper()
        GC_ratio = (seq.count("C")+seq.count("G"))/(len(seq)*1.0)
        segsGC[segID] = GC_ratio
    return(segsGC)

def relativePosInTx(segID, txID, segs_dict, txsDict, exBins_dict):
    seg = segs_dict[segID]
    if txID not in seg.txs:
        return -1
    tx = txsDict[txID]
    ex0 = seg.exs[0]
    pos = 0
    for ex in tx.exs:
        if ex == ex0:
            pos += seg.startLoc - exBins_dict[ex0].start
            break
        else:
            pos += exBins_dict[ex].width
    return pos, tx.length, seg.length
