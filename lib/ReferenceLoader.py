from collections import defaultdict
import os

# A disjoint exonic bin. This class holds disjoint exonic bin data prepared
# in the preprocessing step
class DExon:
    chrome = ""
    exonID = 0
    start = 0
    end = 0
    width = 0
    strand = ""
    seq = ""

    def __init__(self, chrome, exonID, start, end, strand, seq):
        self.chrome = chrome
        self.exonID = int(exonID)
        self.start = int(start)
        self.end = int(end)
        self.width = int(self.end-self.start+1)
        self.strand = strand
        self.seq = seq

    def __str__(self):
        return "<%d : [%d:%d:%d:%s]>" % (self.exonID, self.start, self.end, self.width, self.strand)
    def __repr__(self):
        return self.__str__()

def load_disjointExons(inDir):
    exons = dict()
    with open(os.path.join(inDir, 'disjoint_bins.tsv')) as file_hndl:
        for i, line in enumerate(file_hndl):
            if i == 0:   #Skip header
                continue
            exonID, chrome, start, end, strand, seq = line.strip().split('\t')
            exon = DExon(chrome, exonID, start, end, strand, seq)
            exons[exon.exonID] = exon
    return exons

#######################
#######################

# Transcript class

class TX:
    key_id = 0
    chrome = ""
    txID = ""
    exons = ""

    def __init__(self, key_id, chrome, txID, exons):
        self.key_id = key_id
        self.chrome = chrome
        self.txID = txID
        self.exons = exons

def load_Txs2Exs(inDir):
    numTxs = 0
    genes = defaultdict(list)
    geneIDSorted = set()
    with open(os.path.join(inDir, 'txs2bins.tsv')) as file_hndl:
        for i, line in enumerate(file_hndl):
            if i == 0:   #Skip header
                continue
            key_id, chrome, geneID, txID, strand, exons = line.strip().split('\t')
            genes[geneID].append(TX(key_id, chrome, txID, exons))
            geneIDSorted.add(chrome+":"+geneID)
            numTxs += 1
    geneIDSorted = [genechr.split(":")[1] for genechr in sorted(geneIDSorted)]
    return(genes, geneIDSorted, numTxs)
