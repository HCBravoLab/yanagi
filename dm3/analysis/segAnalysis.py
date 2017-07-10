#######################
## Sequences Lengths ##
#######################
def seq_lengths(exper):
    infile = "..\\"+exper+".fa"
    outfile = exper+"_lens.txt"
    with open(infile, "r") as inf, open(outfile, "w") as outf:
        for i, line in enumerate(inf):
            if i % 2 == 1:  # odd-numbered lines: seqence lines
                outf.write(txID +"\t"+str(len(line)) + "\n")
            else:
                txID = line.strip().split(":")[2]
    print (i+1)/2
    outf.close()
################################
## Filtered Sequences Lengths ##
################################
from collections import defaultdict
def gene_segs_count(exper):
    infile = "..\\"+exper+".fa"
    outfile = exper+".fa.gene_segs_count"
    genes = defaultdict(list)

    with open(infile, "r") as inf, open(outfile, "w") as outf:
        while True:
            header = inf.readline()
            seq = inf.readline()
            if not seq: break  # EOF
            tokens = header.split(":")
            genes[tokens[1]].append(len(seq))

        for geneID in sorted(genes.keys()):
            lens = genes[geneID]
            outf.write(geneID + "\t" + str(len(lens))+"\n")
        
def f_seq_lengths(exper, THRESHOLD):
    infile = "..\\"+exper+".fa"
    outfile = exper+"_lens_f.txt"
    outfile2 = exper+"_lens_fout.txt"
    genes = defaultdict(list)
    with open(infile, "r") as inf, open(outfile, "w") as outf, open(outfile2, "w") as outf2:
        while True:
            header = inf.readline()
            seq = inf.readline()
            if not seq: break  # EOF
            tokens = header.split(":")
            genes[tokens[1]].append(len(seq))

        for geneID in sorted(genes.keys()):
            lens = genes[geneID]
            if len(lens) > THRESHOLD:
                [outf.write(str(i)+"\n") for i in lens]
            else:
                [outf2.write(str(i)+"\n") for i in lens]

##########################
## Filtered Gene Counts ##
##########################
def x_tx_gene_segs_count(exper, X):
    outfile4 = exper+".fa.gene_segs_count_f"+str(X)
    indexes = []
    with open("transcriptome.fa.gene_segs_count", "r") as inf:
        for i, line in enumerate(inf):
            count = line.strip().split("\t")[1]
            if count == x:
                indexes.append(i)

    with open(outfile4, "w") as outf, open(exper+".fa.gene_segs_count", "r") as inf:
        for i, count in enumerate(inf):
            if i in indexes:
                outf.write(count)
    outf.close()

############################
## #Segs/#splices/TX_lens ##
############################
from collections import Counter
def tx_lens_segs_count(exper):
    tx_lens = {}
    with open("transcriptome_lens.txt") as f:
        for line in f:
            tokens = line.strip().split("\t")
            tx_lens[tokens[0]] = tokens[1]
    tx_segs_count = {}
    with open("..\\"+exper+".fa.tx_segs_count") as f:
        for line in f:
            tokens = line.strip().split("\t")
            tx_segs_count[tokens[1]] = tokens[2]
    with open(exper+".fa.tx_lens_segs_count", "w") as f:
        for tx in sorted(tx_lens):
            f.write(tx+"\t"+tx_lens[tx]+"\t"+tx_segs_count[tx]+"\n")

def gene_splices_segs_count(exper):
    gene_splices = {}
    with open("..\\"+exper+".fa.gene_splices") as f:
        for line in f:
            tokens = line.strip().split("\t")
            gene_splices[tokens[0]] = tokens[1]
    gene_segs_count = {}
    with open(exper+".fa.gene_segs_count") as f:
        for line in f:
            tokens = line.strip().split("\t")
            gene_segs_count[tokens[0]] = tokens[1]
    with open(exper+".fa.gene_splices_segs_count", "w") as f:
        for geneID in sorted(gene_splices):
            f.write(geneID+"\t"+gene_splices[geneID]+"\t"+gene_segs_count[geneID]+"\n")

exper = "transcriptome"
print exper
gene_segs_count(exper)
seq_lengths(exper)
f_seq_lengths(exper, 1)

ks = [40, 100, 108, 250, 1000, 10000]
for k in ks:
    exper = "dm3_segs_"+str(k) #"transcriptome"
    print exper
    
    seq_lengths(exper)
    f_seq_lengths(exper, 1)
    
    tx_lens_segs_count(exper)

    gene_segs_count(exper)
    gene_splices_segs_count(exper)

