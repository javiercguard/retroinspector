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

rule get_reference_hg38:
  conda: "../env.yaml"
  output: 
    str(workflow.basedir) + "/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
    "data/ref.fna.gz"
  shell: 
    """
    wget -O {output[1]} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    gunzip {output[1]} > {output[0]}
    """