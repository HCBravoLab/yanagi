from collections import defaultdict
import itertools

class ExSkippingEvent:
    txs = [[], []]

    def __init__(self):
        self.txs = [[], []]

    def addTxs(self, tx1, tx2):
        if tx1 not in self.txs[0]:
            self.txs[0].append(tx1)
        if tx2 not in self.txs[1]:
            self.txs[1].append(tx2)
        

def getConseqEx(ex, exs_list, offset, exceptEx):
    conseqExs = []
    for exs in exs_list:
        if ex in exs:
            i = exs.index(ex)
            i_next = i+offset
            if i_next >= 0 and i_next < len(exs) and exs[i_next] != exceptEx:
                conseqExs.append(exs[i_next])
                return(conseqExs)
            elif i_next < 0 or i_next >= len(exs):
                conseqExs.append("DUMMY_EX")
                return(conseqExs)
    return(conseqExs)

def checkSkippedExon(tx1, tx2, txs, gene, txs2exs, skipEvents2txs):
    exs1 = txs2exs[tx1]
    exs2 = txs2exs[tx2]
    otherExs = [txs2exs[tx] for tx  in txs if tx not in [tx1, tx2]]
    ex_diff = [ex for ex in exs1[1:-1] if not ex in exs2[1:-1]]
    for ex in ex_diff:
        i = exs1.index(ex)
        e1, e2 = exs1[i-1], exs1[i+1]   # e1-ex-e2 sequence in tx1
        # e1-e2 seq in tx2
        cond1 = e1 in exs2 and e2 in exs2 and exs2.index(e1) == (exs2.index(e2) - 1)
        #there is no other ej-ex in any other tx
        cond2 = len(getConseqEx(ex, otherExs, -1, e1)) == 0
        # there is no other ex-ek in any other tx
        cond3 = len(getConseqEx(ex, otherExs, 1, e2)) == 0
        if cond1 and cond2 and cond3:
            event_exs = [str(ex),str(e1),str(e2)]
            event = skipEvents2txs[gene+'\t'+'\t'.join(event_exs)]
            event.gene = gene
            event.addTxs(tx1, tx2)

#Main#1
def generateExSkippingEvents(txdb_dir):
    print "Loading TXDB Data...",
    genes2txs = defaultdict(list)
    txs2exs = {}
    with open(txdb_dir+"txs_exons.txt") as txsf:
        for line in txsf:
            gene, tx, exs, strand = line.strip().split(" ")[2:]
            genes2txs[gene].append(tx)
            txs2exs[tx] = sorted([int(x) for x in exs.split(",")], reverse=strand=="-")
    print "Done!"

    print "Checking for skipped exons..."
    skipEvents2txs = defaultdict(ExSkippingEvent)
    geneCount = len(genes2txs)
    for i, gene in enumerate(sorted(genes2txs.iterkeys())):
        if i % 500 == 0:
            print str(i*100.0/geneCount)+"%"
        txs = genes2txs[gene]
        for txsPair in list(itertools.combinations(txs, 2)):
            checkSkippedExon(txsPair[0], txsPair[1], txs, gene, txs2exs, skipEvents2txs)
            checkSkippedExon(txsPair[1], txsPair[0], txs, gene, txs2exs, skipEvents2txs)
    print "Done!"

    print len(skipEvents2txs)
    print "Writing skipped exons events...",
    with open(txdb_dir+"skippedExons.tsv", "w") as fout:
        for event in sorted(skipEvents2txs.iterkeys()):
            txs = skipEvents2txs[event].txs
            fout.write("\t".join([event, ','.join(txs[0]), ','.join(txs[1])])+"\n")
    print "Done!"

##################
##################

class Seg:
    ID = ""
    name = ""
    header = ""
    txs = set()

    length = 0
    segtype = "E"
    geneID = ""
    startLoc = 0
    strand = "+"

    def __init__(self, header):
        self.header = header

        ID, name, txs, segtype = header.split(" ")
        self.ID = ID
        self.name = name
        
        tokens = name.split(":")
        self.geneID = tokens[1]
        self.startLoc = long(tokens[2])
        self.strand = tokens[-1]
        
        self.segtype = segtype.split(":")[-1]
        self.txs = set(txs.split(":")[1].split(","))

    def __str__(self):
        return "%s:<%s>:%d" % (self.ID, self.name, self.length)
    def __repr__(self):
        return self.__str__()
    
def load_SegmentsLib(seg_file):
    segs_dict = {}
    genes2segs = defaultdict(list)
    with open(seg_file) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i % 2 == 0:
                seg = Seg(line[1:])
                segs_dict[seg.ID] = seg
                genes2segs[seg.geneID].append(seg.ID)
            else:
                seg.length = len(line)
    return(segs_dict, genes2segs)

def segHasExs(seg, exs):
    return ",".join(exs) in seg.name

def segsForSkippedExs(event_line, segs_dict, genes2segs, filterTxs=None):
    gene, skippedEx, preEx, postEx, txsInc, txsExc = event_line.strip().split("\t")
    if filterTxs:
        onlyOneI, indexI = txsListContainsOnlyOneTx(txsInc.split(","), filterTxs)
        onlyOneE, indexE = txsListContainsOnlyOneTx(txsExc.split(","), filterTxs)
        if (not onlyOneI or not onlyOneE): #not txsListContainsAnyTx((txs1+","+txs2).split(","), filterTxs):
            return(None)
        diff = (indexE-indexI)
    else:
        diff = 0
    segs = genes2segs[gene]
    res = [[], [], diff]
    for segID in segs:
        seg = segs_dict[segID]
        if segHasExs(seg, [preEx,skippedEx, postEx]):
            res[0].append(segID)
        elif segHasExs(seg, [preEx, postEx]):
            res[1].append(segID)
    return(res)

#Main#2
def generateExSkippingSegs2Events(segs_filename, txdb_dir, run_label, filterTxs=None):
    print "Loading Segments Data...",
    segs_dict, genes2segs = load_SegmentsLib(segs_filename)
    print "Done!"
    
    print "Getting Segmentns for Skipped Exons...",
    if filterTxs:
        segs_filename += "."+run_label
    with open(txdb_dir+"skippedExons.tsv") as f, \
         open(segs_filename+".segsWithSkippedExons.tsv", "w") as fout:
        fout.write("Gene\tEx\tExBefore\tExAfter\tTxsWithEx\tTxsSkipping\tSegWithEx\tSegSkipping\tNumTxs\tNumSegs\tInc_Exc_ratio\n")
        for line in f:
            segs = segsForSkippedExs(line, segs_dict, genes2segs, filterTxs)
            if not segs: continue
            num_txs = sum([len(l.split(",")) for l in line.strip().split('\t')[4:6]])
            num_segs = len(segs[0])+len(segs[1])
            fout.write("\t".join([line.strip(), ",".join(segs[0]), ",".join(segs[1]), str(num_txs), str(num_segs), str(segs[2])])+"\n")
            if len(segs[0]) == 0 or len(segs[1]) == 0:
                print "Error, segments missing for splice event:"
                print line, segs[0], segs[1]
    print "Done!"


def loadFilterTxs(txs_filename):
    with open(txs_filename) as f:
        txs = f.read().splitlines() 
    return(txs)

def txsListContainsAnyTx(txsList, filterTxs):
    txs = [tx for tx in txsList if tx in filterTxs]
    return len(txs) > 0

def txsListContainsOnlyOneTx(txsList, filterTxs):
    txs = [tx for tx in txsList if tx in filterTxs]
    i = filterTxs.index(txs[0]) if txs else -1
    return (len(txs) == 1, i)

txdb_dir = "hg37/preprocessed/"#"hg37/preprocessed/"#"unit_tests/"#"hg37/preprocessed/"
segs_filename = "hg37/hg37_segs_10000.fa"#"hg37/hg37_segs_108.fa"
#generateExSkippingEvents(txdb_dir)
run_label = "DETxs.logfiltered0"
filterTxs = []#loadFilterTxs("simData/Hs_"+run_label+".txt")
print filterTxs[:10]
generateExSkippingSegs2Events(segs_filename, txdb_dir, run_label)
