import sys

chrs = {}

for line in sys.stdin:
    if line.startswith("#"):
        if line.startswith("##contig="):
            data = line.rstrip().split("=")
            name = data[2].split(",")[0]
            length = int(data[-1][0:-1])
            chrs[name] = length
        print(line, end = '')
    else:
        data = line.split("\t")
        chrom = data[0]
        pos = int(data[1])
        if pos >= chrs[chrom]:
            print("Malformed SV discarded (POS column larger than CONTIG length)", file = sys.stderr)
            print(line, file = sys.stderr, end = '')
        elif pos == 0:
            print("Malformed SV discarded (POS is 0)", file = sys.stderr)
            print(line, file = sys.stderr, end = '')
        else:
            print(line, end = '')
        