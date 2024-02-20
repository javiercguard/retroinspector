#! /usr/env python3

# To be used with python >= 3.7
# Bed file must be tabix-indexed
# VCF should be a SURVIVOR output: this relies on the SVEND field

# from Bio import pairwise2 as pw
from statistics import median
from math import floor
from typing import Type
from pysam import TabixFile
from collections import namedtuple
import sys, gzip, re

vcf = sys.argv[1]
bed = sys.argv[2]
o = sys.argv[3]

def getField(string, field = "END"):
    result = re.search(r"(?:^|;)" + field + "=([^;]+)", string)
    if result:
        return result.group(1)
    else:
        return None # result = None anyways, but this is clearer

bedRecord = namedtuple("bedRecord", # Not a general description of BED records, just BED4 + meName + meFamily
    ["chrom", "start", "end", "name","family"])
def processBedRecord(string):
    s = string.split("\t")

    return bedRecord(s[0], int(s[1]) + 1, int(s[2]), s[4], s[5])

def processDEL (sv, bedPath):
    bed = TabixFile(bedPath, index = bedPath + ".csi")
    
    fields = sv.split("\t")
    chrom = fields[0]
    try:
        pos = int(fields[1])
    except:
        print(sv)
        sys.exit()

    info = fields[7]
    end = int(getField(info))
    size = abs(int(getField(info, "SVLEN")))
    try:
        candidates = [processBedRecord(str(x)) for x in bed.fetch(chrom, pos, end)]
    except ValueError as e:
        print(e, file = sys.stderr)
        return None

    if not len(candidates):
        return None

    best = candidates[0]
    score = size * 2
    for c in candidates:
        currentScore = abs(c.end - end) + abs(c.start - pos)
        if currentScore < score:
            best = c
            score = currentScore

    if not ((best.end - best.start) >= size * 0.85 and (best.end - best.start) <= size * 1.15):
        return None
    else:
        info += f";ME={best.name};MEFAM={best.family};MECOORDS={best.chrom}:{best.start}-{best.end}"
        return "\t".join(fields[0:7] + [info] + fields[8:])

with open(vcf) if vcf.endswith("vcf") else gzip.open(vcf, "rt") as ivcf, \
    open(o, "w") as ovcf:

    from joblib import Parallel, delayed

    outLines = []

    inLines = ivcf.readlines()

    for i in range(len(inLines)): # we first take care of the header
        line = inLines[i]
        if line[0:2] != "##":
            break
        outLines.append(line)
    outLines.append('##INFO=<ID=ME,Number=1,Type=String,Description="Mobile element name per RepeatMasker\'s output">\n')
    outLines.append('##INFO=<ID=MEFAM,Number=1,Type=String,Description="Family of mobile element name per RepeatMasker\'s output">\n')
    outLines.append('##INFO=<ID=MECOORDS,Number=1,Type=String,Description="Coordinates of deleted mobile element sequence name per RepeatMasker\'s output">\n')
    outLines.append(inLines[i]) # Column title row
    i += 1

    inLines = [x for x in inLines[i:] if "SVTYPE=DEL" in x]

    moo = Parallel(n_jobs = 15)(delayed(processDEL)(x, bed) for x in inLines)

    outLines += moo

    outLines = [x for x in outLines if x]

    ovcf.writelines(outLines)
