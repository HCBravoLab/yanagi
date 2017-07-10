import numpy as np
from collections import defaultdict
from sets import Set

import os
from SegGraph import *

# A disjoint Exon. This class holds disjoint exon data prepared
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
        self.start = long(start)
        self.end = long(end)
        self.width = long(self.end-self.start+1)
        self.strand = strand
        self.seq = seq

    def __str__(self):
        return "<%d : [%d:%d:%d:%s]>" % (self.exonID, self.start, self.end, self.width, self.strand)
    def __repr__(self):
        return self.__str__()

    def stream(self):
        return 1 if self.strand == "+" else -1;

    @staticmethod
    def getExon(my_dict, exonID):
        return(my_dict[exonID])

    # Get the list of exons spanning a segment of width _width_
    # starting from location _start_
    @staticmethod
    def getSegExs(my_dict, exs, start, width):
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
            ex = DExon.getExon(my_dict, exs[i])
            start = ex.start if i != 0 else start
            coverage_in_exon = min(rem_width, ex.end-start+1)
            rem_width -= coverage_in_exon
            seg_exs.addEx(ex.exonID)
            i += 1
            used_exs.append(ex)
            used_widths.append(coverage_in_exon)
##            if from_first_ex == 0:
##                from_first_ex = coverage_in_exon
##            else:
##                from_last_ex = coverage_in_exon
        end_pos = start+coverage_in_exon-1
        seg_width = width - rem_width
        #sum_inner_widths = seg_width - from_first_ex - from_last_ex
        #sum_widths_1 = seg_width - from_last_ex

        # Node Type
        ntype = "E"
        ex1 = DExon.getExon(my_dict, exs[0])
        if ex1.start != seg_exs.pos:
            ntype = "J"
        
        return(seg_exs, seg_width, end_pos, used_exs, used_widths, ntype)#sum_inner_widths, sum_widths_1)

    @staticmethod
    def refineSegExs(loc, used_exs, used_widths, seg_width, end_pos, new_loc, k):
        exs_count = len(used_exs)
        inner_coverage = sum(used_widths[1:-1])
        widths_1 = sum(used_widths[:-1])
        if exs_count == 1:
            if seg_width == k:   # Current exon can fit the whole k-length segment
                seg_width = used_exs[0].width
                end_pos = used_exs[0].end
                new_loc = used_exs[0].end + 1 - (k-1)
        else:
            last_ex = used_exs[-1]
            from_last_ex = min(last_ex.width, k -(1 + inner_coverage))
            seg_width = widths_1 + from_last_ex
            end_pos = last_ex.start - 1 + from_last_ex
            if from_last_ex < last_ex.width:   # covered part of the last exon
                new_loc = used_exs[1].start
            else:   # covered the whole of the last exon
                rem = k - (inner_coverage + last_ex.width + 1)
                new_loc = used_exs[0].end + 1 - rem
        return seg_width, end_pos, new_loc
            
        
def load_disjointExons(inDir):
    exons = dict()
    with open(os.path.join(inDir, 'disjoint_exons.txt')) as file_hndl:
        for line in file_hndl:
            exonID, chrome, start, end, strand, seq = line.strip().split(' ')
            exon = DExon(chrome, exonID, start, end, strand, seq)
            exons[exon.exonID] = exon
    return exons

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
    genes = defaultdict(list)
    with open(os.path.join(inDir, 'txs_exons.txt')) as file_hndl:
        for i, line in enumerate(file_hndl):
            key_id, chrome, geneID, txID, exons, strand = line.strip().split(' ')
            genes[geneID].append(TX(key_id, chrome, txID, exons))
    return(genes)
