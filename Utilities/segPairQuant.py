import time
from collections import defaultdict
from collections import Counter
import itertools
import sys

from InputLoader import *

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
        

class ReadAlignments:
    readID = ""
    readLen = 0
    aligns = [[], []]

    def __init__(self, readID, readLen):
        self.readID = readID
        self.readLen = readLen
        self.aligns = [[], []]

    def addAlignment(self, seg, flags, end_idx):
        self.aligns[end_idx].append((seg, flags))

    def __str__(self):
        return "%s:%d:[%s,%d]" % (self.readID, self.readLen, self.aligns[0], self.aligns[1])
    def __repr__(self):
        return self.__str__()

def load_SegmentsLib(seg_file):
    segs_dict = {}
    with open(seg_file) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i % 2 == 0:
                seg = Seg(line[1:])
                segs_dict[seg.ID] = seg
            else:
                seg.length = len(line)
    return(segs_dict)


class TX:
    key_id = 0
    chrome = ""
    geneID = ""
    txID = ""
    exons = []
    ex_widths = []

    def __init__(self, key_id, chrome, geneID, txID, exons, DExons):
        self.key_id = key_id
        self.chrome = chrome
        self.geneID = geneID
        self.txID = txID
        exs = [int(x) for x in exons.strip().split(',')]
        current_exon = DExon.getExon(DExons, exs[0])
        exs = sorted(exs, reverse=(current_exon.strand == "-"))
        self.exons = exs
        for exID in self.exons:
            ex = DExon.getExon(DExons, exID)
            self.ex_widths.append(ex.width)

    def calcLength(self, fromEx, fromLoc, toEx, toLoc, DExons):
        if fromEx == toEx:
            return(toLoc-fromLoc)

        ex_idx = self.exons.index(fromEx)
        length = self.ex_widths[ex_idx] - fromLoc
        for ex_idx in range(self.exons.index(fromEx)+1, self.exons.index(toEx)):
            length += self.ex_widths[ex_idx]
        length += toLoc
        return(length)
        

def load_Txs2Exs(inDir, DExons):
    txs2exs = {}
    with open(os.path.join(inDir, 'txs_exons.txt')) as file_hndl:
        for i, line in enumerate(file_hndl):
            key_id, chrome, geneID, txID, exons = line.strip().split(' ')
            txs2exs[txID] = TX(key_id, chrome, geneID, txID, exons, DExons)
    return(txs2exs)

def read_alignment(f, read_aligns, readIDs, end_idx):
    line = f.readline()
    if not line:
        return ("", False)
    tokens = line.split("\t")
    readID = tokens[0][:-2]
    segID = tokens[2]

    flags = int(tokens[1])

    if segID == "*":
        return(readID, False)

    # Append to readIDs list
    if not readIDs[end_idx] or readID != readIDs[end_idx][-1]:
        readIDs[end_idx].append(readID)
    
    pos = int(tokens[3])
    readLen = len(tokens[9])
    if readID in read_aligns:
        #if flags <= 16:
        read_aligns[readID].addAlignment(segID, flags, end_idx)
    else:
        align = ReadAlignments(readID, readLen)
        #if flags <= 16:
        align.addAlignment(segID, flags, end_idx)
        read_aligns[readID] = align
    return(readID, True)
    

def validLength(fl_params, length):
    f_mean, f_sd, sd_factor = fl_params
    return True #((length < (f_mean + sd_factor * f_sd)) and (length > (f_mean - sd_factor * f_sd)))

def validOrientation(flags1, flags2):
    return flags1+flags2 == 16

#@profile
def processReads(readIDs, read_aligns, segs_dict, segPairs_counts, segPairs_txs):
    readsMapped = [0, 0, 0] # [mapped, notmapped_due_txs, notmapped_due_ori]
    c=0
    for readID in readIDs:
        aligns = read_aligns.pop(readID)
        mapped = False
        passed_txs = False
        c += len(aligns.aligns[0]+aligns.aligns[1])
        for align in itertools.product(aligns.aligns[0], aligns.aligns[1]):
            segID1, segID2 = align[0][0], align[1][0]
            pairkey = segID1 + "_" + segID2
            if pairkey in segPairs_txs:
                txs = segPairs_txs[pairkey]
            else:
                seg1 = segs_dict[segID1]
                seg2 = segs_dict[segID2]
                txs = (seg1.txs & seg2.txs)
                segPairs_txs[pairkey] = txs
            if len(txs) < 1: continue
            if not validOrientation(align[0][1], align[1][1]):
                passed_txs = True
                #print align[0][1], align[1][1]
                continue
            segPairs_counts[pairkey] += 1
            mapped = True
            break
        if mapped:
            readsMapped[0] += 1
        elif passed_txs:
            readsMapped[2] += 1
        else:
            readsMapped[1] += 1
    print (c*1.0)/len(readIDs)
    return readsMapped

#@profile
def main(sampleID, rootDIR, process_count):
    index_ref = "dm3_segs_10000"
    #bam_R1 = rootDIR+"/"+sampleID+"__"+index_ref+"_R1/paligns_noh.sam"
    #bam_R2 = rootDIR+"/"+sampleID+"__"+index_ref+"_R2/paligns_noh.sam"
    bam_R1 = rootDIR+"/"+sampleID+"_1_noh.sam"
    bam_R2 = rootDIR+"/"+sampleID+"_2_noh.sam"    
    seg_file = index_ref+".fa"
    #seg_inDir = 'hg37/preprocessed'
    out_counts = rootDIR+"/seg_counts/rapmap_"+sampleID+"__"+index_ref+"_seg_counts.tsv"

    #bam_R1="test_readAnalysis/o2_1.sam"
    #bam_R2="test_readAnalysis/o2_2.sam"
    #out_counts="test_readAnalysis/o2_counts.txt"

    #bam_R1="output/HPGL0445__hg38_segs_108_R1/o2_1.sam"
    #bam_R2="output/HPGL0445__hg38_segs_108_R2/o2_2.sam"
    #out_counts="output/seg_counts/o2_counts.txt"

    #process_count = 1000000
  
    print "Loading Segments Lib...",
    start_t = time.time()
    #DExons = load_disjointExons(seg_inDir)
    #txs2exs = load_Txs2Exs(seg_inDir, DExons)
    segs_dict = load_SegmentsLib(seg_file)
    print "Done!"
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Estimate Fragment Lengths..."
    start_t = time.time()
    fl_params = [181, 62.01, 1]

    print "Reading Read alignments..."
    start_t = time.time()
    read_aligns = {}
    segPairs_counts = Counter()
    segPairs_txs = {}
    #TODO clean up based on sorted readIDs
    readsMapped = [0,0,0]
    with open(bam_R1) as f1, open(bam_R2) as f2:
        done = [False, False]
        readIDs = [[], []]
        i = 0
        while not all(done):
            if not done[0]:
                readID, aligned = read_alignment(f1, read_aligns, readIDs, 0)
                done[0] = not readID
            if not done[1]:
                readID, aligned = read_alignment(f2, read_aligns, readIDs, 1)
                done[1] = not readID

            # Process some records
            if len(readIDs[0]) > process_count and len(readIDs[1]) > process_count:
                process_readIDs = readIDs[0][:process_count]
                readIDs[0] = readIDs[0][process_count:]
                readIDs[1] = readIDs[1][process_count:]
                print i, len(read_aligns)
                start_t_2 = time.time()
                rMapped = processReads(process_readIDs, read_aligns, segs_dict, segPairs_counts, segPairs_txs)
                readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
                elapsed = time.time() - start_t_2
                print "Processed:", len(process_readIDs), len(read_aligns), rMapped, elapsed
            i+=1
    print i, len(read_aligns)
    rMapped = processReads(readIDs[0], read_aligns, segs_dict, segPairs_counts, segPairs_txs)
    readsMapped = [readsMapped[0]+rMapped[0], readsMapped[1]+rMapped[1], readsMapped[2]+rMapped[2]]
    print "Processed:", len(process_readIDs), len(read_aligns), rMapped
    print "Done!"
    print "Mapped Reads:", readsMapped[0], "Unmapped Reads:", readsMapped[1]+readsMapped[2], \
          "Unmapped Due Txs:", readsMapped[1], "Unmapped Due Direction:", readsMapped[2]
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed

    print "Writing Segments Counts..."
    start_t = time.time()
    with open(out_counts, "w") as f:
        f.write("SEG1ID\tSEG2ID\tcount\tSEGTYPES\tGENE\tSEG1LEN\tSEG2LEN\tSEG1SLoc\tSEG2SLoc\tTXS\n")
        
        for segPair in sorted(segPairs_counts.iterkeys()):
            count = segPairs_counts[segPair]
            segs = [segs_dict[segID] for segID in segPair.split("_")]
            types = segs[0].segtype
            types += segs[1].segtype if segs[0].ID != segs[1].ID else ""
            line = "\t".join([segs[0].ID, segs[1].ID, str(count), types,
                              segs[0].geneID,
                              str(segs[0].length), str(segs[1].length),
                              str(segs[0].startLoc), str(segs[1].startLoc), ','.join(segPairs_txs[segPair])])
            f.write(line + "\n")
    elapsed = time.time() - start_t
    print "Elapsed Time: ", elapsed


#log_f = open("pairends_log.txt", "w")
#log_f2 = open("pairends_log2.txt", "w")
#log_f3 = open("pairends_log3.txt", "w")

sampleID = sys.argv[1]#"HPGL0365"
rootDir = sys.argv[2]
process_count = int(sys.argv[3])
main(sampleID, rootDir, process_count)

#log_f.close()
#log_f2.close()
#log_f3.close()
