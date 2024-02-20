#! /usr/env/python
import subprocess

tmpFasta = snakemake.output.fasta
svafDb = snakemake.params.dbPrefix
svafTypes = snakemake.params.svafTypes

# Read rm.bed and save SVA-F in a temp fasta file
with open(snakemake.input[0]) as inputBed, \
        open(tmpFasta, "w") as fasta:
    fastaLines = []
    for line in inputBed.readlines():
        line = line.rstrip().split("\t")
        name = line[3]
        seqId = line[6]
        sequence = line[17]

        fastaLines.append(f">{seqId}\n{sequence}\n")
    fasta.writelines(fastaLines)

# Create blastn from MAST2 repeat
# Run blastn with fasta and MAST2 db
print(f"svaf db: {svafDb}")
subprocess.run(["makeblastdb", "-in", svafTypes,
               "-out", svafDb, "-dbtype", "nucl"])
subprocess.run(["blastn", "-task", "blastn", "-evalue", "1e-20", "-db", svafDb,
               "-query", tmpFasta, "-outfmt", "6", "-out", 
               snakemake.output.blast])
# Create dict from blastn
blastResults = {}
with open(snakemake.output.blast) as blast:
    for line in blast.readlines():
        line = line.split("\t")
        seqId = line[0]
        if seqId not in blastResults:
            blastResults[seqId] = line[1]

# Read rm.bed again and if seqId is in dict, replace SVA-F with SVA-F1
# and save it to output
with open(snakemake.input[0]) as inputBed, \
    open(snakemake.output.bed, "w") as o:
    lines = []
    for line in inputBed.readlines():
        line = line.rstrip().split("\t")
        name = line[3]
        seqId = line[6]
        if seqId in blastResults:
            line[3] = blastResults[seqId]
        lines.append("\t".join(line) + "\n")
    o.writelines(lines)

