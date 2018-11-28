
import lib.utils

class Seg:
    ID = ""
    name = ""
    header = ""
    txs = set()

    chrm = ""
    length = 0
    segtype = "E"
    geneID = ""
    startLoc = 0
    endLoc = 0
    exs = []
    strand = "+"

    seq = ''

    def __init__(self, header):
        if header[0] == '>':
            header = header[1:]
        self.header = header

        ID, name, txs, segtype = header.split(" ")
        self.ID = ID
        self.name = name
        
        tokens = name.split(":")
        self.chrm = tokens[0]
        self.geneID = tokens[1]
        self.startLoc = int(tokens[2])
        self.exs = [int(ex) for ex in tokens[3][1:-1].split(",")]
        self.endLoc = int(tokens[4])
        self.strand = tokens[-1]
        
        self.segtype = segtype.split(":")[-1]
        self.txs = set(txs.split(":")[1].split(","))

        self.seq = ''

    def __str__(self):
        return "%s:<%s>:%d" % (self.ID, self.name, self.length)
    def __repr__(self):
        return self.__str__()

def load_SegmentsLib(seg_file, loadSeqs=False):
    print("Loading Segments Library")
    segs_dict = {}
    with open(seg_file) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i % 2 == 0:
                seg = Seg(line[1:])
                segs_dict[seg.ID] = seg
            else:
                seg.length = len(line)
                if loadSeqs:
                    seg.seq = line
    return(segs_dict)

