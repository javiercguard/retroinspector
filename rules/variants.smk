rule run_svim:
  conda: "../env.yaml"
  log:
    "logs/svim/{sample}.log"
  input:
    genome = config["referenceGenome"],
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai"
  output:
    directory("tmp/svim/{sample}"),
  shell:
    """
    svim alignment --sequence_alleles --read_names \
    {output} {input.bam} {input.genome} 2> {log}
    """

rule process_svim:
  conda: "../env.yaml"
  log:
    "logs/svim/{sample}.process.log"
  input:
    "tmp/svim/{sample}"
  output:
    "variants/svim/{sample}.svim.vcf.gz",
    "variants/svim/{sample}.svim.vcf.gz.csi"
  params:
    support = config["minimumReadSupport"],
    fixScript = str(workflow.basedir) + "/scripts/fixVCF.py",
  shell:
    """
    cat {input}/variants.vcf | \
    gawk -v 'OFS=\\t' \
    '{{if (substr($0, 1, 1) == "#") {{print}} \
    else {{ match($8, /SVTYPE=([^;]+)/, a); \
      if (a[1] == "INS" || a[1] == "BND") {{print}} \
      else {{\
        $5="<"a[1]">";$4="N";print}} \
      }} \
    }}' | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{print}} else {{ match($8, /SUPPORT=([0-9]+)/, a); if (a[1] >= {params.support}) {{print}} }} }}' | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{sub(/SUPPORT,/, "RE,", $0); print}} else {{ sub(/SUPPORT=/, "RE=", $8); print }} }}' | \
    python {params.fixScript} | \
    bcftools sort -O z -o {output[0]}; \
    bcftools index {output[0]}
    """

rule run_sniffles:
  conda: "../env.yaml"
  log:
    "logs/sniffles/{sample}.log"
  threads: config['threads']
  input:
    genome = config["referenceGenome"],
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai"
  output:
    temp("variants/{sample}.sniffles.temp.vcf")
    # "variants/{sample}.sniffles.temp.vcf"
  shell:
    """
    sniffles -t {threads} \
    -s 2 --report_BND -m {input.bam} -v {output} 2> {log}
    """

rule create_sample_file:
  output:
    temp("tmp/{sample}.sampleFileForHeader.txt")
  shell:
    """
    echo {wildcards.sample} > {output}
    """

rule process_sniffles:
  conda: "../env.yaml"
  log:
    "logs/sniffles/{sample}.process.log"
  input:
    "variants/{sample}.sniffles.temp.vcf",
    "tmp/{sample}.sampleFileForHeader.txt"
  output:
    "variants/sniffles/{sample}.sniffles.vcf.gz",
    "variants/sniffles/{sample}.sniffles.vcf.gz.csi"
  params:
    support = config["minimumReadSupport"]
  shell:
    """
    cat {input[0]} | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{print}} else {{ match($8, /RE=([0-9]+)/, a); if (a[1] >= {params.support}) {{print}} }} }}' 2> {log} | \
    gawk -v 'OFS=\\t' '{{\
    if ($0 ~ /^#CHROM/) {{\
      print "##FILTER=<ID=STRANDBIAS,Description=\\"Strand is biased if Strandbias_pval< 0.01.\\">\\n"$0\
    }} else {{print}} }}' | \
    bcftools reheader -s {input[1]} 2>> {log} | \
    bcftools sort -O z -o {output[0]} 2>> {log} ; \
    bcftools index {output[0]}
    """

rule run_sniffles2:
  conda: "../env_sniffles2.yaml" # Sniffles2 woud require its own environment or displace sniffles
  log:
    "logs/sniffles2/{sample}.log"
  threads: config['threads']
  input:
    genome = config["referenceGenome"],
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai"
  output:
    "variants/sniffles2/{sample}.sniffles2.temp.vcf"
  shell:
    """
    sniffles --threads {threads} --output-rnames \
    --input {input.bam} --vcf {output} &> {log}
    """

rule process_sniffles2:
  conda: "../env.yaml"
  log:
    "logs/sniffles2/{sample}.proc.log"
  input:
    rules.run_sniffles2.output[0]
  output:
    "variants/sniffles2/{sample}.sniffles2.vcf.gz",
    "variants/sniffles2/{sample}.sniffles2.vcf.gz.csi"
  params:
    support = config["minimumReadSupport"],
    fixScript = str(workflow.basedir) + "/scripts/fixVCF.py",
  shell:
    """
    cat {input[0]} | \
    gawk -v 'OFS=\\t' \
    '{{if (substr($0, 1, 1) == "#") {{print}} \
    else {{ match($8, /SVTYPE=([^;]+)/, a); \
      if (a[1] == "INS" || a[1] == "BND") {{print}} \
      else {{\
        $5="<"a[1]">";$4="N";print}} \
      }} \
    }}' | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{print}} else {{ match($8, /SUPPORT=([0-9]+)/, a); if (a[1]+0 >= {params.support}) {{print}} }} }}' | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{sub(/SUPPORT,/, "RE,", $0); print}} else {{ sub(/SUPPORT=/, "RE=", $8); print }} }}' | \
    python {params.fixScript} | \
    bcftools sort -O z -o {output[0]}; \
    bcftools index {output[0]}
    """


rule run_cutesv:
  conda: "../env.yaml"
  log:
    "logs/cutesv/{sample}.log"
  threads: config['threads']
  input:
    genome = config["referenceGenome"],
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai"
  output:
    # temp(directory("tmp/cutesv/{sample}")),
    directory("tmp/cutesv/{sample}"),
    temp("tmp/cutesv/{sample}/{sample}.cutesv.withseq.vcf"),
    "tmp/cutesv/{sample}/{sample}.cutesv.unproc.vcf",
  shell:
    """
    cuteSV -t {threads} \
    --max_cluster_bias_INS 100 -s2 -L -1 \
    --report_readid  --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
    {input.bam} {input.genome} {output[1]} {output[0]} 2> {log}
    cat {output[1]} 2> {log} | \
    gawk -v 'OFS=\\t' \
    '{{if (substr($0, 1, 1) == "#") {{print}} \
    else {{\
      match($8, /SVTYPE=([^;]+)/, a); \
      if (a[1] == "INS" || a[1] == "BND") {{print}} 
      else {{$5="<"a[1]">";$4="N";print}} }} }}' 2>> {log} > {output[2]}
    """

rule process_cutesv:
  conda: "../env.yaml"
  log:
    "logs/cutesv/{sample}.process.log"
  input:
    "tmp/cutesv/{sample}/{sample}.cutesv.unproc.vcf",
    "tmp/cutesv/{sample}",
  output:
    "variants/cutesv/{sample}.cutesv.vcf.gz",
    "variants/cutesv/{sample}.cutesv.vcf.gz.csi"
  params:
    support = config["minimumReadSupport"],
    fixScript = str(workflow.basedir) + "/scripts/fixVCF.py",
  shell:
    """
    cat {input[0]} 2> {log} | \
    gawk -v 'OFS=\\t' '{{if (substr($0, 1, 1) == "#") {{print}} else {{ match($8, /RE=([0-9]+)/, a); if (a[1] >= {params.support}) {{print}} }} }}' 2>> {log} | \
    python {params.fixScript} 2>> {log} | \
    bcftools sort -O z -o {output[0]} 2>> {log};
    bcftools index {output[0]} 2>> {log}
    """
