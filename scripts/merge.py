#! /usr/env python3
from statistics import median
from math import floor, ceil
from pprint import pp
from pysam import VariantFile
from collections import namedtuple
from fuzzywuzzy import fuzz as fz
import sys
import argparse
import traceback

signature = "merge.py"

aln = namedtuple("aln", "score")


def getGtQv(var, sampleN, hasGT, hasQV):
    gt = (None, None)
    qv = -1
    if not hasGT:
        return gt, qv
    if sampleN == 1:
        gt = var.samples[0]["GT"]
    elif sampleN > 1:
        if hasQV:
            genotypes = [(var.samples[sample]["GT"],
                          float(var.samples[sample]["QV"]) if isNumber(var.samples[sample]["QV"]) else -1)
                         for sample in var.samples]
            genotypes.sort(key=lambda x: x[1], reverse=True)
            gt = genotypes[0][0]
            qv = genotypes[0][1]
        else:
            gt = var.samples[0]["GT"]
    return gt, qv


def isNumber(x):
    try:
        float(x)
        return True
    except:
        return False


def globAlign(seq1, seq2):
    """
    Globally align 2 sequences without gap penalties.
    str, str -> [namedtuple])
    """
    r = aln(score=fz.ratio(seq1, seq2))
    return r


def replaceChar(str, what, idx):
    """
    Puts what at index idx of string str. str, int, str -> str
    """
    return str[0:idx] + what + str[idx + 1:]


def toGenotype(inputDict):
    """
    dict -> str
    """
    fields = ("GT", "PSV", "LN", "DR", "ST", "QV",
              "TY", "ID", "RAL", "AAL", "CO", "SC")
    result = ""
    for i in fields:
        result += str(inputDict[i]) + ":"
    return result[0:-1]  # remove trailing colon


def getInfoValue(sv, what="SVLEN", idx=0, toInt=True):
    result = ""
    try:
        result = sv.info[what][idx]
    except TypeError:
        result = sv.info[what]
    if toInt:
        try:
            result = int(result)
        except ValueError:
            print("The following record had no length, this may cause trouble!")
            print(sv)
            result = None
    return result


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-samples", dest="samples", nargs="+",
                    required=True, help="Samples for sample columns")
parser.add_argument("-vcf", dest="files", nargs="+",
                    required=True, help="Indexed VCF files")
parser.add_argument("-o", dest="out", required=True, help="Output VCF file")
parser.add_argument("-d", dest="dist", default=0, type=int,
                    help="Maximum dist between starting poitns for SVs")
parser.add_argument("--multi-sample-input", "-msi", dest="singleSample", default=True,
                    action="store_true", help="The input files correspond to one sample each")

args = parser.parse_args()

sampleList = args.samples
vcfs = args.files
outFile = args.out
distance = args.dist
singleSample = args.singleSample
considerSeq = 1

print(vcfs)

vcfs = [VariantFile(vcf) for vcf in vcfs]

consumed = set()  # For SVs that have been matched
# emptyGeno = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN" # Not found in sample
# formatStr = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO"

out = VariantFile(outFile, "w")
for i in sampleList:
    out.header.add_sample(i)
for i in vcfs[0].header.records:
    if i.type == "CONTIG":
        out.header.add_record(i)

out.header.add_meta("INFO", items=[("ID", "CIEND"), ("Number", 2), (
    "Type", "String"), ("Description", "PE confidence interval around END")])
out.header.add_meta("INFO", items=[("ID", "CIPOS"), ("Number", 2), (
    "Type", "String"), ("Description", "PE confidence interval around POS")])
out.header.add_meta("INFO", items=[("ID", "CHR2"), ("Number", 1), (
    "Type", "String"), ("Description", "Chrom. for END coord. for translocations")])
out.header.add_meta("INFO", items=[("ID", "END"), ("Number", 1), (
    "Type", "Integer"), ("Description", "End position of the structural variant")])
out.header.add_meta("INFO", items=[("ID", "MAPQ"), ("Number", 1), (
    "Type", "Integer"), ("Description", "Median mapping quality of paired-ends")])
out.header.add_meta("INFO", items=[(
    "ID", "RE"), ("Number", 1), ("Type", "Integer"), ("Description", "read support")])
out.header.add_meta("INFO", items=[(
    "ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "Length of the SV")])
out.header.add_meta("INFO", items=[("ID", "MAXDIFF"), ("Number", 1), ("Type", "Integer"), (
    "Description", "Max. difference in start coordinates between samples.")])
out.header.add_meta("INFO", items=[("ID", "SVMETHOD"), ("Number", 1), (
    "Type", "String"), ("Description", "Method for generating this merged VCF file")])
out.header.add_meta("INFO", items=[(
    "ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "Type of the SV")])
out.header.add_meta("INFO", items=[("ID", "SUPP_VEC"), ("Number", 1), (
    "Type", "String"), ("Description", "Vector of supporting samples")])
out.header.add_meta("INFO", items=[("ID", "SUPP"), ("Number", 1), (
    "Type", "Integer"), ("Description", "Number of samples supporting the variant")])
out.header.add_meta("INFO", items=[("ID", "STRANDS"), ("Number", 1), (
    "Type", "String"), ("Description", "Indicating the direction of the reads with resp")])

out.header.add_meta("FORMAT", items=[
                    ("ID", "GT"), ("Number", 1), ("Type", "String"), ("Description", "Genotype")])
out.header.add_meta("FORMAT", items=[("ID", "PSV"), ("Number", 1), (
    "Type", "String"), ("Description", "Previous support vector")])
out.header.add_meta("FORMAT", items=[(
    "ID", "LN"), ("Number", 1), ("Type", "Integer"), ("Description", "predicted length")])
out.header.add_meta("FORMAT", items=[("ID", "DR"), ("Number", 2), ("Type", "Integer"), (
    "Description", "#supporting reference,variant reads in that order")])
out.header.add_meta("FORMAT", items=[(
    "ID", "ST"), ("Number", 1), ("Type", "String"), ("Description", "Strand of SVs")])
out.header.add_meta("FORMAT", items=[("ID", "QV"), ("Number", 1), ("Type", "String"), (
    "Description", "Quality values: if not defined a . otherwise the r")]) # string to keep compatibility with SURVIVOR's output
out.header.add_meta("FORMAT", items=[
                    ("ID", "TY"), ("Number", 1), ("Type", "String"), ("Description", "Types")])
out.header.add_meta("FORMAT", items=[(
    "ID", "ID"), ("Number", 1), ("Type", "String"), ("Description", "Variant ID from input")])
out.header.add_meta("FORMAT", items=[("ID", "RAL"), ("Number", 1), (
    "Type", "String"), ("Description", "Reference allele sequence reported from input")])
out.header.add_meta("FORMAT", items=[("ID", "AAL"), ("Number", 1), ("Type", "String"), (
    "Description", "Alternative allele sequence reported from input")])
out.header.add_meta("FORMAT", items=[(
    "ID", "CO"), ("Number", 1), ("Type", "String"), ("Description", "Coordinates")])
out.header.add_meta("FORMAT", items=[("ID", "SC"), ("Number", 1), (
    "Type", "Integer"), ("Description", "Levenshtein distance ratio to first sample")])

out.header.add_meta("ALT", items=[("ID", "INS"), ("Description", "Insertion")])

out.header.add_meta("METHOD", value=signature)

sampleParams = {x: {"supportKey": "", "method": ""} for x in sampleList}

for i in range(len(sampleList)):
    vcf = vcfs[i]
    sample = sampleList[i]

    exampleRecord = ""  # Will be a VariantRecord instance

    for i in vcf:
        exampleRecord = i
        break  # We take the first one as an example
    if f"##METHOD={signature}" in str(vcf.header):
        sampleParams[sample]["method"] = signature
    if "RE" in exampleRecord.info.keys():
        sampleParams[sample]["supportKey"] = "RE"
    if "PE" in exampleRecord.info.keys():
        sampleParams[sample]["supportKey"] = "PE"
    if "SUPPORT" in exampleRecord.info.keys() and sampleParams[sample] != signature:
        sampleParams[sample]["supportKey"] = "SUPPORT"

emptyGenotype = {
    "GT": (0, 0),
    "PSV": ".",
    "LN": None,
    "DR": (None, None),
    "ST": ".",
    "QV": ".",
    "TY": ".",
    "ID": ".",
    "RAL": ".",
    "AAL": ".",
    "CO": ".",
    "SC": None,
}
genotypeFields = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO:SC".split(":")

# poscheck = 21871943
n = 0
for i in range(len(vcfs) - 1):

    vcf = vcfs[i]
    sample = sampleList[i]
    others = vcfs[i + 1:]

    # print(sample)
    sampleN = len(list(vcf.header.samples))
    hasGT = "GT" in vcf.header.formats.keys()
    hasQV = "QV" in vcf.header.formats.keys()

    for var in vcf.fetch():

        uid = f"{sample}_{var.id}"
        # if poscheck-10 <= var.pos <= poscheck+10:
        #         print(var.chrom, var.pos)

        if var.info["SVTYPE"] != "INS":
            continue

        if uid in consumed:
            continue

        consumed.add(uid)

        readSupport = 0
        if sampleParams[sample]["method"] == signature:
            readSupport = max([
                int(x) if isNumber(x) else -1
                for x in [geno["DR"][1] for geno in var.samples.values()]
            ])
        elif sampleParams[sample]["supportKey"]:
            readSupport = var.info[sampleParams[sample]["supportKey"]]
        # See output header for description.

        altField = var.alts[0].replace(":", "_") if var.alts else var.alts
        if not altField:
            print(var)
            sys.exit("INS without ALT?")

        try:
            gt, qv = (None, None), None
            if singleSample:
                gt, qv = getGtQv(var, sampleN, hasGT, hasQV)
            varFormat = {
                "GT": gt,
                "PSV": var.info['SUPP_VEC']
                if "SUPP_VEC" in var.info.keys() else "NAN",
                "LN": getInfoValue(var),
                "DR": (None, readSupport),
                "ST": "NAN",
                "QV": str(qv),
                "TY": "NaN",
                "ID": var.id,
                "RAL": var.ref,
                "AAL": altField,
                "CO": f"{var.chrom}_{var.pos}_{var.stop}",
                "SC": 100,
                'start': var.pos,
                'end': var.stop,
            }
        except Exception as e:
            print(traceback.format_exc())
            print(var.alts)
            sys.exit(
                f"Error when creating genotype column at sample {sample} for var\n{var}."
            )
        result = {
            "info":
            dict([
                ("SUPP", 1),
                ("SUPP_VEC", replaceChar("0" * len(sampleList), "1", i)),
                ("SVLEN", getInfoValue(var)),
                ("SVTYPE", var.info["SVTYPE"]),
                ("SVMETHOD", "merge.py"),
            ]),
            "samples": {
                sample: varFormat
            }  # length of samples
        }

        svLen = getInfoValue(var)
        svLenToCompare = ceil(svLen * 0.15)
        svScoreToCompare = 60

        for j in range(len(vcfs)):
            sampleToCompare = sampleList[j]
            if j == i:  # no intra-VCF merging
                continue
            other = vcfs[j]
            candidates = []  # candidates are reset for every VCF
            for candidate in other.fetch(var.chrom, max(var.pos - distance, 0),
                                         var.pos + distance):
                if candidate.pos not in range(max(var.pos - distance, 0),
                                              var.pos + distance):
                    continue
                uid2 = sampleToCompare + "_" + candidate.id
                if uid2 in consumed:
                    continue
                candidates.append(candidate)
            candidates = list(
                filter(
                    lambda x:
                        x.info["SVTYPE"] == var.info["SVTYPE"] and
                    abs(svLen - getInfoValue(x, toInt=1)) <= svLenToCompare,
                    candidates)
            )
            # if var.pos == poscheck:
            #     print(var.chrom, var.pos)
            #     print(f"# Initial candidates: {len(candidates)}")
            if not candidates:
                result["samples"][sampleToCompare] = emptyGenotype
                continue  # To next VCF
            try:
                candidates = sorted(
                    candidates,
                    key=lambda x: abs(
                        getInfoValue(x, toInt=1) - var.info["SVLEN"]))  # !
            except:
                print(var.info["SVLEN"])
                sys.exit(
                    f"Could not sort candidates in sample {sample} for variant\n{var}"
                )
            if considerSeq:
                alns = [{
                    "candidate": x,
                    "aln": globAlign(altField, x.alts[0])
                } for x in candidates]
                chosen = max(alns, key=lambda x: x["aln"].score)
                cand = chosen["candidate"]
                # if var.pos == poscheck:
                #     print("chosen")
                #     print(chosen)
                # print(chosen['aln'].score)
                if chosen["aln"].score < svScoreToCompare:
                    result["samples"][sampleToCompare] = emptyGenotype
                    continue  # To next VCF
            else:
                cand = candidates[0]
                chosen = {"cand": cand, "aln": 0}
            # if var.pos == poscheck:
            #         print("we have chosen")
            #         pp(chosen)
            uid2 = sampleToCompare + "_" + cand.id
            consumed.add(uid2)
            result["info"]["SUPP"] += 1
            result["info"]["SUPP_VEC"] = replaceChar(
                result["info"]["SUPP_VEC"], "1", j)
            readSupport = 0
            if sampleParams[sampleToCompare]["method"] == signature:
                readSupport = max([
                    int(x) if isNumber(x) else -1
                    for x in [geno["DR"][1] for geno in cand.samples.values()]
                ])
            elif sampleParams[sampleToCompare]["supportKey"]:
                readSupport = cand.info[sampleParams[sampleToCompare]
                                        ["supportKey"]]

            gt, qv = (None, None), None
            if singleSample:
                gt, qv = getGtQv(cand, len(list(cand.samples)), "GT" in cand.header.formats.keys(), "QV" in cand.header.formats.keys())
            result["samples"][sampleToCompare] = {
                "GT":
                gt,
                "PSV":
                cand.info['SUPP_VEC']
                if "SUPP_VEC" in cand.info.keys() else "NAN",
                "LN":
                getInfoValue(cand),
                "DR": (None, readSupport),
                "ST":
                "NAN",
                "QV":
                str(qv),
                "TY":
                "NaN",
                "ID":
                cand.id,
                "RAL":
                cand.ref,
                "AAL":
                cand.alts[0].replace(":", "_"),
                "CO":
                f"{cand.chrom}_{cand.pos}_{cand.stop}",
                "SC":
                chosen["aln"].score
                if isinstance(chosen["aln"], tuple) else chosen["aln"],
                'start':
                cand.pos,
                'end':
                cand.stop,
            }

        startToUse = median([
            x["start"] for x in result["samples"].values()
            if "start" in x.keys()
        ])
        endToUse = median([
            x["end"] for x in result["samples"].values() if "end" in x.keys()
        ])
        lenToUse = floor(
            median([x["LN"] for x in result["samples"].values() if x["LN"]]))
        starts = [
            x["start"] for x in result["samples"].values()
            if "start" in x.keys()
        ]
        maxDiffInCoords = max(starts) - min(starts)
        result["info"]["MAXDIFF"] = maxDiffInCoords
        result["info"]["SVLEN"] = lenToUse
        record = out.new_record(
            contig=var.chrom,
            start=startToUse - 1,
            stop=endToUse,
            alleles=(var.ref, var.alts[0]),
            filter="PASS",
            id=str(n),
        )
        # if var.pos == 172571922:
        #     pp(result)
        n += 1
        for field, value in result["info"].items():
            try:
                record.info[field] = value
            except:
                sys.exit(
                    f"Error when setting INFO field {field} for sample {sample} at variant\n{var}"
                )
        for sampleName, sampleGeno in result["samples"].items():
            for field in genotypeFields:
                try:
                    record.samples[sampleName][field] = sampleGeno[field]
                except:
                    print(traceback.format_exc())
                    sys.exit(
                        f"Error when setting {field} with value {sampleGeno[field]} at column for sample {sampleName} for var\n{var}"
                    )

        out.write(record)

print("Finished all but last VCF")

# last vcf
sample = sampleList[-1]
print("last sample:", sample)
for var in vcfs[-1].fetch():  # Add the unmatched vars for the last one
    i = len(vcfs)
    # For index in lists, etc.
    uid = f"{sample}_{var.id}"

    if var.info["SVTYPE"] != "INS":
        continue

    if uid in consumed:
        continue

    consumed.add(uid)

    readSupport = 0
    if sampleParams[sample]["method"] == signature:
        readSupport = max([
            int(x)
            if isNumber(x) else -1
            for x in [geno["DR"][1] for geno in var.samples.values()]
        ])
    elif sampleParams[sample]["supportKey"]:
        readSupport = var.info[sampleParams[sample]["supportKey"]]
    try:
        gt, qv = (None, None), None
        if singleSample:
            gt, qv = getGtQv(var, sampleN, hasGT, hasQV)
    except:
        print(traceback.format_exc())
        sys.exit(
            f"Encountered error when processing genotypes for (last) sample {sample} at var\n{var}")
    varFormat = {
        "GT": gt,
        "PSV": var.info['SUPP_VEC'] if "SUPP_VEC" in var.info.keys() else "NAN",
        "LN": getInfoValue(var),
        "DR": (None, readSupport),
        "ST": "NAN",
        "QV": str(qv),
        "TY": "NaN",
        "ID": var.id,
        "RAL": var.ref,
        "AAL": var.alts[0].replace(":", "_"),
        "CO": f"{var.chrom}_{var.pos}_{var.stop}",
        "SC": globAlign(var.alts[0], var.alts[0]).score,
    }
    result = {
        "info": dict([
            ("SUPP", 1),
            ("SUPP_VEC", replaceChar("0" * len(sampleList), "1", i)),
            ("SVLEN", getInfoValue(var)),
            ("SVTYPE", var.info["SVTYPE"]),
            ("SVMETHOD", "merge.py"),
        ]),
        "samples": {sample: varFormat}  # length of samples
    }

    for j in range(len(vcfs) - 1):  # other samples
        # These variants have already been checked so we only need to add the empty genotypes
        sampleToCompare = sampleList[j]
        result["samples"][sampleToCompare] = emptyGenotype
    try:
        record = out.new_record(contig=var.chrom, start=var.start - 1, stop=var.stop,
                                alleles=(var.ref, var.alts[0]), filter="PASS",
                                id=str(n),
                                )
    except:
        print(traceback.format_exc())
        sys.exit(f"Error when writing record for last sample for var\n{var}")
    n += 1
    for field, value in result["info"].items():
        try:
            record.info[field] = value
        except:
            print(traceback.format_exc())
            sys.exit(
                    f"Error when setting INFO field {field} for sample {sample} at variant\n{var}"
                )
    for sampleName, sampleGeno in result["samples"].items():
        for field, value in sampleGeno.items():
            try:
                record.samples[sampleName][field] = value
            except Exception as e:
                print(traceback.format_exc())
                sys.exit(
                    f"Error when setting {field} with value {sampleGeno[field]} for sample {sampleName} (at last var) for var\n{var}"
                )

    out.write(record)

out.close()
