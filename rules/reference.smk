rule get_reference_repeats:
  conda: "../env.yaml"
  output: 
    temp("tmp/hg38.fa.out.gz"),
    temp("data/repeatsReferenceTE.bed"),
    temp("data/repeatsReferenceTE.bed.gz"),
    temp("data/repeatsReferenceTE.bed.gz.csi"),
  shell:
    """
    wget -O {output[0]} \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz
    zcat {output[0]} | tail -n +4 | \
    grep -E "(S|L)INE|Retroposon|LTR|DNA" | \
    awk -v'OFS=\\t' '{{print $5,$6,$7,$9,$10,$11}}' > {output[1]}
    bgzip -k {output[1]}
    tabix -p bed --csi {output[2]}
    """

rule prepare_te_nooverlap:
  conda: "../env.yaml"
  input:
    "data/repeatsReferenceTE.bed.gz"
  output:
    temp("tmp/repeatsReferenceTENoOverlap.bed")
  shell:
    """
    bedtools merge -c 4,5,6 -o collapse,collapse,collapse \
      -i {input[0]} > \
      {output[0]}
    """