
from collections import defaultdict
from collections import Counter

from lib.SegGraph import *
from lib.ReferenceLoader import *
from lib.GTFGenerator import *
from lib.Seg2EventMapper import generateEventsSegsIOE

from tqdm import tqdm
import time, os

############################################
## Useful Functions while creating SegGraph
############################################

# Get the list of exons spanning a segment of width _width_
# starting from location _start_
def getSegExs(exBins_dict, exs, start, width):
    rem_width = width
    coverage_in_exon = 0
    used_exs = []
    used_widths = []
    #from_first_ex = 0
    #from_last_ex = 0
    seg_exs = SG_Key()
    seg_exs.pos = start
    i = 0
    while rem_width > 0 and i < len(exs):
        ex = exBins_dict[exs[i]]
        start = ex.start if i != 0 else start
        coverage_in_exon = min(rem_width, ex.end-start+1)
        rem_width -= coverage_in_exon
        seg_exs.addEx(ex.exonID)
        i += 1
        used_exs.append(ex)
        used_widths.append(coverage_in_exon)
    end_pos = start+coverage_in_exon-1
    seg_width = width - rem_width

    # Node Type
    ntype = "E"
    ex1 = exBins_dict[exs[0]]
    if ex1.start != seg_exs.pos:
        ntype = "J"
    
    return(seg_exs, seg_width, end_pos, used_exs, used_widths, ntype)

#################################

# Function handles very short exonic bins shorter than L

def refineSegExs(loc, used_exs, used_widths, seg_width, end_pos, new_loc, L):
    exs_count = len(used_exs)
    inner_coverage = sum(used_widths[1:-1])
    widths_1 = sum(used_widths[:-1])
    if exs_count == 1:
        if seg_width == L:   # Current exon can fit the whole L-lengthed segment
            seg_width = used_exs[0].width
            end_pos = used_exs[0].end
            new_loc = used_exs[0].end + 1 - (L-1)
    else:
        last_ex = used_exs[-1]
        from_last_ex = min(last_ex.width, L -(1 + inner_coverage))
        seg_width = widths_1 + from_last_ex
        end_pos = last_ex.start - 1 + from_last_ex
        if from_last_ex < last_ex.width:   # covered part of the last exon
            new_loc = used_exs[1].start
        else:   # covered the whole of the last exon
            rem = L - (inner_coverage + last_ex.width + 1)
            new_loc = used_exs[0].end + 1 - rem
    return seg_width, end_pos, new_loc


#################################
######## SegmentsGraph ##########
#################################

# A final segment that can be built from aggregating a sequence of consecutive segments
class SegContig:

    # Static variable
    SEGS_COUNT = 0
    
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


    def finalizeSegment(self, chrome, geneID, exBins_dict):
        exs = self.seg_key.getExs()
        strand = exBins_dict[exs[0]].strand
        exs = sorted(exs, reverse=(strand == "-"))
        
        seg_start = self.seg_key.pos
        self.seq = ""
        for i, exonID in enumerate(exs):
            ex = exBins_dict[exonID]
            if i == 0:  # First Exon
                from_pos = seg_start - ex.start
                # Set segment strand
                strand = ex.strand
            else:
                from_pos = 0
                if strand != ex.strand: #Sanity check for segment strand
                    print("ERROR, Segment covers exons of different strands in Gene "+geneID)
            if i == len(exs)-1: # Last Exon
                to_pos = self.end - ex.start + 1
            else:
                to_pos = ex.width
            
            self.seq += ex.seq[from_pos:to_pos]

        exs = ','.join(map(str, exs))
        txs = ','.join(sorted(self.txIDs))

        SegContig.SEGS_COUNT += 1
        segID = "SEG%07d" % (SegContig.SEGS_COUNT)
        identifier = ">%s %s:%s:%d:(%s):%d:%s TXs:%s segtype:%s" % (segID, chrome, geneID,
                                                    seg_start, exs, self.end, strand, txs, self.segtype)
        identifier = ">%s" % (segID)
        meta = '\t'.join([segID, chrome, geneID, txs, exs, str(seg_start), str(self.end), strand])+'\n'
        return segID, meta, identifier + '\n' + self.seq + '\n'

def parseSegGraph(SG, startNodes, redundantNodes,
                  chrome, geneID, DExons):
    output = ""
    out_meta = ""
    out_extra = ""
    for key in sorted(startNodes, key=lambda x:(x.exs_str, x.pos)):
        current = SG[key.exsKey()][key.pos]
        if key in redundantNodes:
            continue
        segment = SegContig()
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
        output = output + seg_str
        out_meta += meta
        out_extra += segID + " " + '|'.join(segment.nodeIDs) + "\n"
    return(output, out_meta, out_extra)


#############################
######## Main Logic #########
#############################

# Main Logic
#kernprof-script.py -l -v SegGraphCreator.py
#@profile
def createSegments(L, inDir, outname, eventsFiles=None, shreded=False):
    print("Loading preprocessed Annotation...",)
    start_time = time.time()
    # Load input
    DExons = load_disjointExons(inDir)
    txs2exons, geneIDSorted, numTxs = load_Txs2Exs(inDir)
    print("ET: ", time.time() - start_time)

    fulloutname = os.path.join(inDir, outname+".fa")
    output_file = open(fulloutname, "w")
    outf_meta = open(fulloutname+".meta", "w")
    outf_meta.write('\t'.join(["segID", "chrom", "geneID", "txAnnIDs",
                               "binIDs", "st", "end", "strand"])+'\n')  #Meta header
    
    output_gtf = open(fulloutname+".gtf", "w")
    writeExonicBinsToGTF(output_gtf, outname, DExons, txs2exons, geneIDSorted)
    
    chr_start_time = time.time()
    chrome = ""

    SegContig.SEGS_COUNT = 0

    num_shortTxs = 0
    SG_total_time = 0
    Segs_total_time = 0
    # Iterate over genes
    for it in tqdm(range(len(geneIDSorted))):
        geneID = geneIDSorted[it]
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
            tx_len = sum([DExons[ex].width for ex in exs])
            if tx_len < L:
                num_shortTxs += 1
            current_exon = DExons[exs[0]]
            exs = sorted(exs, reverse=(current_exon.strand == "-"))
            
            current_exon = DExons[exs[0]]
            loc = current_exon.start
            end_tx = DExons[exs[-1]].end
            prev = DUMMY_START_NODE
            seg_width = L
            # main segmentMaker loop
            while True:
                seg_key, seg_width, end_pos, used_exs, used_widths, ntype = getSegExs(DExons, exs, loc, L)
                exs_count = len(used_exs)#exs_count = seg_key.exsLen()
                seg_exs = seg_key.exsKey()
                new_loc = loc + 1
                seg_width, end_pos, new_loc = refineSegExs(loc, used_exs, used_widths,
                                                           seg_width, end_pos, new_loc, L)
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
                if seg_width < L and not(prev == DUMMY_START_NODE):
                    redundantNodes.add(node.key)
                prev = node
                loc = new_loc
                if seg_width < L:
                    break
                if loc > current_exon.end or loc < current_exon.start:
                    exs.pop(0)
                    if len(exs) == 0:
                        break
                    else:
                        current_exon = DExons[exs[0]]
                        loc = current_exon.start
        #print("Done creating SG")
        SG_total_time += time.time() - start_time
        start_time = time.time()
        # Generate Segments in a gene by parseing the generated segments graph
        output, out_meta, output_extra = parseSegGraph(SG, startNodes, redundantNodes,
                                                       chrome, geneID, DExons)
        output_file.write(output)
        outf_meta.write(out_meta)
        writeSegmentsToGTF(output_gtf, outname, out_meta.split("\n"))
        
        Segs_total_time += time.time() - start_time
    output_file.close()
    outf_meta.close()
    print("Creating SG ET:", SG_total_time)
    print("Creating Segments ET:", Segs_total_time)
    print("Processed", len(geneIDSorted), "Genes,", numTxs, "Transcripts,", SegContig.SEGS_COUNT, "Segments")
    print("Total of", num_shortTxs, "Transcripts shorter than", L)
    if eventsFiles:
        print("Mapping Segments to Alternative Splicing Events")
        for eventsFile in eventsFiles:
            print("Processing File", eventsFile)
            generateEventsSegsIOE(fulloutname, DExons, eventsFile,
                                  fulloutname[:-3]+"_"+os.path.basename(eventsFile))
            print("Done")
        

if __name__ == '__main__':
    createSegments(101, "C:/yanagi_new/output/", "hg37_segs101.fa")
