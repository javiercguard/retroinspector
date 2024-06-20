#! /usr/env/python2

from pysam import VariantFile

def main(snakemake):
    vcfFile = snakemake.input.vcf
    bedFile = snakemake.output.bed
    vcf = VariantFile(vcfFile)
    with open(bedFile, "w") as bed:
        for var in vcf.fetch():
            bed.write(f"{var.chrom}\t{var.pos - 61}\t{var.stop + 60}\t{var.id}\n")

if __name__ == "__main__":
    main(snakemake)
