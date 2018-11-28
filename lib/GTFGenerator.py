
from collections import defaultdict

from lib.ReferenceLoader import DExon
from lib.Segments import Seg

def writeExonicBinsToGTF(fout, fname, DExons, txs2exons, geneIDs):
    exs2txs = defaultdict(list)
    exs2genes = defaultdict(set)
    for gene in geneIDs:
        txs = txs2exons[gene]
        for tx in txs:
            [(exs2txs[int(x)].append(tx.txID), exs2genes[int(x)].add(gene)) for x in tx.exons.strip().split(',')]

    for exID in sorted(DExons):
        ex = DExons[exID]
        line = '\t'.join([ex.chrome, fname, "exonic_bins",
                          str(ex.start), str(ex.end),
                          ".", ex.strand, ".",
                          "gene_id \""+'+'.join(sorted(exs2genes[exID]))+"\"; "+
                          "entry_id \""+str(exID)+"\"; "+
                          "transcripts \""+'+'.join(sorted(exs2txs[exID]))+"\";"])+"\n"
        fout.write(line)
    print(len(DExons.keys()), "exonic bins")      


def writeSegmentsToGTF(fout, fname, segsHeaders):
    for header in segsHeaders:
        seg = Seg(header)
        line = '\t'.join([seg.chrm, fname, "segment",
                          str(min(seg.startLoc, seg.endLoc)),
                          str(max(seg.startLoc, seg.endLoc)),
                          ".", seg.strand, ".",
                          "gene_id \""+seg.geneID+"\"; "+
                          "entry_id \""+seg.ID+"\"; "+
                          "exonic_bin_ids \""+'+'.join([str(ex) for ex in sorted(seg.exs)])+"\"; "+
                          "transcripts \""+'+'.join(sorted(seg.txs))+"\";"])+"\n"
        fout.write(line)
