
from collections import defaultdict

from lib.ReferenceLoader import DExon
from lib.utils import *

def writeExonicBinsToGTF(fout, fname, DExons, txs2exons, geneIDs):
    exs2txs = defaultdict(list)
    exs2genes = defaultdict(set)
    for gene in geneIDs:
        txs = txs2exons[gene]
        for tx in txs:
            [(exs2txs[int(x)].append(tx.txID), exs2genes[int(x)].add(gene)) for x in tx.exons.strip().split(',')]

    for exID in sorted(DExons):
        ex = DExons[exID]
        line = '\t'.join([ex.chrome, fname, "exonic_bin",
                          str(ex.start), str(ex.end),
                          ".", ex.strand, ".",
                          "gene_id \""+'+'.join(sorted(exs2genes[exID]))+"\"; "+
                          "entry_id \""+str(exID)+"\"; "+
                          "transcripts \""+'+'.join(sorted(exs2txs[exID]))+"\";"])+"\n"
        fout.write(line)
    print(len(DExons.keys()), "exonic bins")      


def writeSegmentsToGTF(fout, fname, segsHeaders):
    for header in segsHeaders:
        if len(header)==0:
            continue
        segID, chrome, geneID, txs, exs, seg_start, seg_end, strand = header.split("\t")
        line = '\t'.join([chrome, fname, "segment",
                          str(min(int(seg_start), int(seg_end))),
                          str(max(int(seg_start), int(seg_end))),
                          ".", strand, ".",
                          "gene_id \""+geneID+"\"; "+
                          "entry_id \""+segID+"\"; "+
                          "exonic_bin_ids \""+exs.replace(",", "+")+"\"; "+
                          "transcripts \""+txs.replace(",", "+")+"\";"])+"\n"
        fout.write(line)
