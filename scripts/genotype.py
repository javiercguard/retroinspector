#! /usr/env/python3

import math

import gzip
from pysam import VariantFile


def readCoverage(bedPath):
    result = {}
    with gzip.open(bedPath, "rt") as bed:
        for line in bed:
            chrom, pos, end, id, coverage = line.rstrip().split("\t")
            result[id] = float(coverage)
    return result


def probRatio(q1, q2):
    try:
        log = math.log10(q1 / q2)
        return log
    except ValueError:
        return 0


def binomSniffles(x, size, chance):
    """dbinom with a shortcut. It's possible to skip a calculation because the values will be compared to each other.
    Adapted from Sniffles2.

    Parameters
    ----------
    x : _int_
        _Number of reads supporting the INS._
    size : _int_
        _Coverage value_
    chance : _float_
        _Probability_
    """
    return chance ** x * (1 - chance) ** (size - x)


def genotype(supportValue, coverageValue,
             chances=[0.05, 0.5, 0.95],
             genotypes=[0, 1, 2],
             maxNorm=250):
    genotypeValues = [(0, 0), (0, 1), (1, 1)]
    if supportValue > coverageValue:
        coverageValue = supportValue
    maxValue = max(supportValue, coverageValue)
    if maxValue > maxNorm:
        factorNorm = 250 / maxValue
        supportValue *= factorNorm
        coverageValue *= factorNorm
    probs = [binomSniffles(supportValue, coverageValue, chance)
             for chance in chances]
    probs = [(name, value) for name, value in zip(genotypes, probs)]
    probs.sort(key=lambda x: x[1], reverse=True)
    sumProbs = sum([x[1] for x in probs])
    probs = [(x[0], x[1] / sumProbs) for x in probs]
    q1 = probs[0][1]
    q2 = probs[1][1]
    qual = min(60, math.floor(-10 * probRatio(q2, q1)))
    geno = genotypeValues[probs[0][0]]
    return geno, qual
    # return {key: val for key, val in zip(genotypes, )}

def isNumber(x):
    try:
        float(x)
        return True
    except:
        return False

def main(snakemake):
    coverages = readCoverage(snakemake.input.coverage)
    with VariantFile(snakemake.input.vcf) as vcf, \
            VariantFile(snakemake.output.vcf, 'w', header=vcf.header) as out:
        samples = list(vcf.header.samples)
        for record in vcf.fetch():
            try:
                support = [record.samples[x]["DR"][1]
                        for x in samples]  # second value, support for SV
                genotypes = [genotype(float(sup), coverages[record.id]) if isNumber(sup) else ((None, None), None)
                            for sup in support ]
                for i, (geno, qual) in enumerate(genotypes):
                    record.samples[i]["GT"] = geno
                    record.samples[i]["QV"] = qual
            except:
                print(record)
                import sys, traceback
                print(traceback.format_exc())
                sys.exit()
            out.write(record)


if __name__ == "__main__":
    main(snakemake)
