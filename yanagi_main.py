
from collections import defaultdict
from collections import Counter

from ReferenceLoader import *
from SegGraph import *

#from TransReconstructor import *

#################################
######## SegmentsGraph ##########
#################################

# A final segment that can be built from aggregating a sequence of consecutive segments
class Segment:
    seg_key = SG_Key()  # has fields: exs, start
    end = 0
    seq = ""

    txIDs = set()
    segtype = "E"

    nodeIDs = []    # SG_Nodes merged by the maximality property

    def __init__(self):
        self.seg_key = SG_Key()
        self.end = 0
        self.seq = ""

        self.txIDs = set()
        self.segtype = "E"

        self.nodeIDs = []

    def appendNode(self, node):
        if self.seg_key.pos == 0:
            self.seg_key.pos = node.key.pos
            self.txIDs = node.getTxIDs()
            self.segtype = node.ntype
        else:
            self.txIDs &= node.getTxIDs()
            
        self.end = node.end
        self.seg_key.addExs(node.key.exs)
        self.nodeIDs.append(node.keyStr())


    def finalizeSegment(self, chrome, geneID, exon_dict):
        exs = self.seg_key.getExs()
        strand = DExon.getExon(exon_dict, exs[0]).strand
        exs = sorted(exs, reverse=(strand == "-"))
        
        seg_start = self.seg_key.pos
        self.seq = ""
        for i, exonID in enumerate(exs):
            ex = DExon.getExon(exon_dict, exonID)
            if i == 0:  # First Exon
                from_pos = seg_start - ex.start
                # Set segment strand
                strand = ex.strand
            else:
                from_pos = 0
                if strand != ex.strand: #Sanity check for segment strand
                    print("ERROR, Segment covers exons of different strands")
            if i == len(exs)-1: # Last Exon
                to_pos = self.end - ex.start + 1
            else:
                to_pos = ex.width
            
            #print("%s > %d:%d %s" % (ex.key, from_pos, to_pos, ex.seq[from_pos:to_pos]))
            self.seq += ex.seq[from_pos:to_pos]

        exs = ','.join(map(str, exs))
        txs = ','.join(sorted(self.txIDs))

        global segsCount
        segsCount += 1
        segID = "SEG%07d" % (segsCount)
        identifier = ">%s %s:%s:%d:(%s):%d:%s TXs:%s segtype:%s" % (segID, chrome, geneID,
                                                    seg_start, exs, self.end, strand, txs, self.segtype)
        #identifier = ">%s" % (segID)
        meta = '\t'.join([segID, chrome, geneID, txs, str(seg_start), str(self.end)])+'\n'
        return segID, meta, identifier + '\n' + self.seq + '\n'

#############################
######## Main Logic #########
#############################

# Main Logic
#kernprof-script.py -l -v SegMaker_main.py
#@profile
def createSG(k, inDir, outname):
    print("Loading preprocessed data...",)
    start_time = time.time()
    # Load input
    DExons = load_disjointExons(inDir)
    txs2exons, geneIDSorted = load_Txs2Exs(inDir)
    #geneIDs = txs2exons.keys()
    print("ET: ", time.time() - start_time)

    output_file = open(outname, "w")
    outf_meta = open(outname+".meta", "w")
    outf_extra = open(outname+".extra", "w")
    # Statistics files
    spl_count_file = open(outname+".gene_splices", "w")
    num_segs_len_txs_file = open(outname+".tx_segs_count", "w")
    
    chr_start_time = time.time()
    chrome = ""

    tot = 0

    SG_total_time = 0
    Segs_total_time = 0
    total_reductions = 0
    # Iterate over genes
    for geneID in geneIDSorted:
        start_time = time.time()
        start_time = time.time()
        #Init SegmentGraph
        SG = defaultdict(dict)
        startNodes = set()
        redundantNodes = set()

        txs = txs2exons[geneID]
        
        DUMMY_START_NODE = SG_Node(SG_Key(), 0, 0)
        # Iterate over transcripts of the gene
        for tx in txs:
            chrome = tx.chrome
            exs = [int(x) for x in tx.exons.strip().split(',')]
            tx_len = sum([DExon.getExon(DExons, ex).width for ex in exs])
            if tx_len < k:
                tot = tot+1
            current_exon = DExon.getExon(DExons, exs[0])
            exs = sorted(exs, reverse=(current_exon.strand == "-"))
            
            current_exon = DExon.getExon(DExons, exs[0])
            stream = current_exon.stream()
            loc = current_exon.start
            end_tx = DExon.getExon(DExons, exs[-1]).end
            prev = DUMMY_START_NODE
            seg_width = k
            # main segmentMaker loop
            while True:
                seg_key, seg_width, end_pos, used_exs, used_widths, ntype = DExon.getSegExs(DExons, exs, loc, k)
                exs_count = len(used_exs)#exs_count = seg_key.exsLen()
                seg_exs = seg_key.exsKey()
                new_loc = loc + 1
                seg_width, end_pos, new_loc = DExon.refineSegExs(loc, used_exs, used_widths, seg_width, end_pos, new_loc, k)
                node = None
                if loc in SG[seg_exs]:
                    node = SG[seg_exs][loc]
                else:
                    node = SG_Node(seg_key, seg_width, end_pos)
                    # Set ntype: e.g. E, J, ...
                    node.ntype = ntype
                node.addTxID(tx.txID)
                prev.my_next.add(node.key)
                SG[seg_exs][loc] = node
                if prev.color != node.color:
                    for next_key in prev.my_next:
                        n = SG[next_key.exsKey()][next_key.pos]
                        n.isStart = True
                        SG[next_key.exsKey()][next_key.pos] = n
                        startNodes.add(next_key)
                if shreded:
                    node.isStart = True
                    startNodes.add(node.key)
                if seg_width < k and not(prev == DUMMY_START_NODE):
                    redundantNodes.add(node.key)
                prev = node
                loc = new_loc
                if seg_width < k:
                    break
                if loc > current_exon.end or loc < current_exon.start:
                    exs.pop(0)
                    if len(exs) == 0:
                        break
                    else:
                        current_exon = DExon.getExon(DExons, exs[0])
                        stream = current_exon.stream()
                        loc = current_exon.start
        #print("Done creating SG")
        SG_total_time += time.time() - start_time
        start_time = time.time()
        output, out_meta, output_extra, txs_num_segs, reductions = generateSegments(startNodes, redundantNodes, chrome, geneID, SG, DExons, k)
        total_reductions += reductions #- len(txs)*(k-1)
        output_file.write(output)
        outf_meta.write(out_meta)
        outf_extra.write(output_extra)
        # Statistics
        #splices_count = count_splices(geneID, txs)
        #spl_count_file.write(geneID+"\t"+str(splices_count)+"\n")
        #for tx, c in txs_num_segs.items():
        #    num_segs_len_txs_file.write(geneID + "\t" + tx + "\t" + str(c) + "\n")
        Segs_total_time += time.time() - start_time
    output_file.close()
    outf_meta.close()
    outf_extra.close()
    spl_count_file.close()
    num_segs_len_txs_file.close()
    print("Creating SG ET:", SG_total_time)
    print("Creating Segments ET:", Segs_total_time)
    print("Total Reduction:", str(total_reductions/1000000.0), "Gbp")
    print("Total of ", tot, "Transcripts shorter than", k)

def generateSegments(startNodes, redundantNodes, chrome, geneID, SG, DExons, k):
    reductions = 0
    output = ""
    out_meta = ""
    out_extra = ""
    # Statistics
    txs_num_segs = Counter()
    for key in sorted(startNodes, key=lambda x:(x.exs_str, x.pos)):
        current = SG[key.exsKey()][key.pos]
        if key in redundantNodes:
            continue
        #reductions += (len(current.getTxIDs())-1)*current.width
        #reductions += (len(current.getTxIDs())-1)*(k-1)
        segment = Segment()
        done = False
        while not done:
            segment.appendNode(current)

            has_next = len(current.my_next) > 0
            done = not has_next         
            
            if has_next:
                next_key = next(iter(current.my_next))
                next_node = SG[next_key.exsKey()][next_key.pos]
                done = next_node.isStart
                current = next_node
        segID, meta, seg_str = segment.finalizeSegment(chrome, geneID, DExons)
        reductions += (len(segment.txIDs)-1)*len(segment.seq)
        output = output + seg_str
        out_meta += meta
        out_extra += segID + " " + '|'.join(segment.nodeIDs) + "\n"
        for tx in segment.txIDs:
            txs_num_segs[tx] += 1 
    return(output, out_meta, out_extra, txs_num_segs, reductions)

def count_splices(geneID, txs2exons):
    exs_dict = defaultdict(list)
    for tx in txs2exons:
        exs = [int(x) for x in tx.exons.strip().split(',')]
        for ex in exs:
            exs_dict[ex].append(tx.txID)
    splices_count = 0
    exs = sorted(exs_dict.keys())
    for i in range(0,len(exs)-1):
        diff = set(exs_dict[exs[i]]) - set(exs_dict[exs[i+1]])
##        if geneID == "1":
##            print exs[i], set(exs_dict[exs[i]]), set(exs_dict[exs[i+1]]), diff
        splices_count += len(diff)
    return(splices_count)

import time
import sys
#from pathlib import Path

#Ls = [40, 108, 1000, 10000]#sys.argv[2]
exper = "dm3"
Ls = [108]
shreded=True
for L in Ls:
    print("L=" + str(L))
    inDir = exper+'\\preprocessed'
    shred = "shreded_" if shreded else ""
    outfile = exper+'\\'+shred+exper+'_segs_'+str(L)+'.fa'#sys.argv[3]
    print(outfile)

    total_start = time.time()
    segsCount = 0
    createSG(L, inDir, outfile)
    elapsed = time.time() - total_start
    print("Elapsed Time: ", elapsed)
