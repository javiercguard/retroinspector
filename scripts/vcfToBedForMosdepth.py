#! /usr/bin/env python3

from pysam import VariantFile

def main(snakemake):
    svType = ""
    if "svType" in snakemake.params.keys():
        svType = snakemake.params["svType"]
    vcfFile = snakemake.input.vcf
    bedFile = snakemake.output.bed
    with open(bedFile, "w") as bed, VariantFile(vcfFile) as vcf:
        for var in vcf.fetch():
            if svType and var.info["SVTYPE"] != svType:
                continue
            bed.write(f"{var.chrom}\t{max(0, var.pos - 61)}\t{var.stop + 60}\t{var.id}\n")

if __name__ == "__main__":
    main(snakemake)
