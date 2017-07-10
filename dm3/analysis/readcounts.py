from collections import Counter
import sys

def report_count(count, genes, read_len, outf):
    num_genes = len(genes)
    outf.write(str(count)+"\t"+str(num_genes))
    if read_len != 101:
        outf.write("\t"+str(read_len))
    outf.write("\n")
    counts[count] += 1
    counts_per_gene[count-num_genes+1] += 1
    genes_counts[num_genes] += 1
    

infile = sys.argv[1] #"po.txt"
outfile = sys.argv[2] #"read_counts.txt"
outf = open(outfile, "w")
outf.write("Count\tGenesCount\tReadLen\n")

## Statistics for multimapping counts
statfile = outfile+".stat"
statf = open(statfile, "w")
counts = Counter()
counts_per_gene = Counter()
genes_counts = Counter()

count = 0
read_len = 0
unmapped_count = 0
genes = set()
readName =  None
i = 0
for line in open(infile):
    tokens = line.split("\t")
    read = tokens[0]
    if tokens[2] == '*':
        unmapped_count += 1
        continue
    gene = tokens[2].split(":")[1]
    read_len = len(tokens[9])
    if readName and readName != read:
        report_count(count, genes, read_len, outf)
        
        count = 1
        genes.clear()
        genes.add(gene)
    else:
        count += 1
        genes.add(gene)
    readName = read
    if i % 1000000 == 0:
        print i
    i += 1

report_count(count, genes, read_len, outf)

outf.close()

statf.write("Counts Frequencies:\n")
statf.write(','.join(map(str, counts.keys()))+"\n")
statf.write(','.join(map(str, counts.values()))+"\n")

statf.write("Counts/Gene Frequencies:\n")
statf.write(','.join(map(str, counts_per_gene.keys()))+"\n")
statf.write(','.join(map(str, counts_per_gene.values()))+"\n")

statf.write("GenesCounts Frequencies:\n")
statf.write(','.join(map(str, genes_counts.keys()))+"\n")
statf.write(','.join(map(str, genes_counts.values()))+"\n")

statf.write("Mappings\tUnmapped\n"+str(i)+"\t"+str(unmapped_count)+"\n")

statf.close()
