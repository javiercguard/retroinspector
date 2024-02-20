#! /usr/bin/env python

# Should be run with python >= 3.6
import sys, gzip, subprocess
from math import floor

def getInfoValue (infoField, fieldName):
        if ";" + infoField + "=": # infoField is not the first field
            # thus its preceeded by a semicolon
            # this avoids getting CIPOSEND when looking for END, etc.
            # or END_SOMETHING, etc.
            result = infoField[infoField.index(";" + fieldName + "=") + len(fieldName) + 2:]
            # the "+ 2" is for the "=" and the ";"
        else: # its the first field
            #the first summand might as well be 0
            result = infoField[infoField.index(fieldName + "=") + len(fieldName) + 1:]
            # the "+ 1" is for the "="
        result = result[0:result.index(";")] 
        # This works even if there is not trailing semicolon
        # because python is just like that
        return result

def revCom (seq):
    return seq.translate("".maketrans("ATCG", "TAGC"))[::-1]

inputVcf = sys.argv[1]
outputDir = sys.argv[2]
species = sys.argv[3]
threads = int(sys.argv[4])

fileName = inputVcf[inputVcf.rfind("/") + 1:]
fileName = fileName[:fileName.rfind(".vcf")] # wont contain a trailing period
outFasta = f"{outputDir}/{fileName}.fa"
rmOut = f"{outputDir}/{fileName}.fa.out"
insBed = f"{outputDir}/{fileName}.rm.unproc.bed"

threads = 16

insertions = {} # dict of dicts

insertionLines = []

with open(inputVcf) if inputVcf.endswith(".vcf") else gzip.open(inputVcf, "rt") as fileIterator, \
    open(outFasta, "w") as fastaHandler:

    fastaContent = []
    insN = 0
    
    lines = fileIterator.readlines()

    for i in range(len(lines)):

        line = lines[i]

        if line[0] == "#": # We skip the header
            continue

        lineFields = line.rstrip().split("\t")

        chrom, pos, svId, ref, alt, qual, filterResult, info = lineFields[0:8]
        genotypes = lineFields[8:]
        if alt == "<INS>": # we cant work with these
            continue
        svType = getInfoValue(info, "SVTYPE")
        if "INS" not in svType:
            continue
        
        insId = "seq" + str(insN)

        insertions[insId] = {
            "chrom": chrom,
            "pos": str(int(pos) - 1),
            "end": str(pos), #!
            "lineN": i,
            "repeats": [],
        }
        fastaContent.append(f">{insId}\n{alt}\n")

        insN += 1

    fastaHandler.writelines(fastaContent)
    insertionLines = lines

print(len(list(insertions.keys())))

# sys.exit()

repeatMaskerCommand = ["RepeatMasker", "-e", "rmblast", "-species", species, "-no_is",
    "-pa", str(floor(threads / 4)), # divided by 4 because RepMasker (allegedly) multiplies it by 4
    outFasta]

print("Running RepeatMasker")

print(subprocess.list2cmdline(repeatMaskerCommand))
subprocess.run(repeatMaskerCommand)

print("Finished RepeatMasker")

# Get RepeatMasker results

with open(rmOut) as f:
    for i in range(4): # Skip the first 3 lines, RepeatMakser header
        f.readline()
    for line in f.readlines():

        values = list(filter(lambda x: x != "", line.rstrip().split(" ")))

        data = values[0:11] # repeat masker output line split by whitespace
        insType = data[10]

        fields = ["score", "divergencePercentage", "delPercentage", "insPercentage",
        "seqId", "begin", "end", "left", "strain", "repeat"]

        result = {}
        for i in range(len(fields)):
            result[fields[i]] = data[i]
        for x in ["score", "begin", "end"]:
            result[x] = int(result[x])
        result["left"] = int(result["left"][1:len(result["left"])-1])
        for x in ["divergencePercentage", "delPercentage", "insPercentage"]:
            result[x] = float(result[x])
        typeData = insType.split("/")
        if len(typeData) == 2:
            result["class"], result["subclass"] = typeData
        elif len(typeData) > 2:
            result["class"] = "/".join(typeData)
            result["subclass"] = "-"
        else:
            result["class"] = typeData[0]
            result["subclass"] = "-"
        insertions[result["seqId"]]["repeats"].append(result)

print(len(insertions))
for key, x in insertions.copy().items(): # Remove Insertions dicts without repeats
    # print(x)
    if not len(x["repeats"]):
        del insertions[key]
    else:
        x["repeats"] = max(x["repeats"], key = lambda x: x["score"])
print(len(insertions))

with open(insBed, "w") as t:
    lines = []
    indexes = range(12) # number of fields in genotype number
    for seqId, ins in insertions.items():
        fields = insertionLines[ins["lineN"]].rstrip().split("\t")
        usefulInfo = [fields[4], fields[7], *fields[9:]] # The ALT, INFO and sample columns
        for i in range(2, len(usefulInfo)): # remove alt allele from sample columns
            usefulInfo[i] = ":".join(
                list(
                    map(
                        lambda x,y: x if y != 9 else ".", 
                        usefulInfo[i].split(":"), 
                        indexes)
                    )
            )

        line = "\t".join([
            ins["chrom"], ins["pos"], 
            ins["end"],
            ins["repeats"]["repeat"],
            "1", # score is now a meaningless value
            "*",
            seqId,
            str(ins["repeats"]["score"]), ins["repeats"]["class"], ins["repeats"]["subclass"],
            str(ins["repeats"]["begin"]), str(ins["repeats"]["end"]),
            str(ins["repeats"]["left"]), ins["repeats"]["strain"], 
            str(ins["repeats"]["divergencePercentage"]),
            str(ins["repeats"]["delPercentage"]), 
            str(ins["repeats"]["insPercentage"]), 
            # ins["info"], ins["alt"], 
            *usefulInfo,
            ]) + "\n"
        lines.append(line)

    t.writelines(lines)

