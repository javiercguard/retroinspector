def getFastqs (wildcards):
  return config["samples"][wildcards.sample]

rule align:
  conda: "../env.yaml"
  log:
    "logs/minimap/{sample}.log"
  threads: config['threads'] - 1
  input:
    config["referenceGenome"],
    getFastqs
  output:
    temp("tmp/sam/{sample}.sam")
  shell:
    """
    minimap2 -x map-ont -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
    --MD -a -t {threads} {input} > {output} 2> {log}
    """

rule sortAndCompress:
  conda: "../env.yaml"
  log:
    "logs/samtools/{sample}.sort.log"
  input:
    "tmp/sam/{sample}.sam"
  output:
    "alns/{sample}.bam"
  params:
    threads = 16
  shell:
    "samtools sort -@ {params.threads} -O BAM -o {output} {input} 2> {log}"

rule index:
  conda: "../env.yaml"
  log:
    "logs/samtools/{sample}.index.log"
  input:
    "alns/{sample}.bam"
  output:
    "alns/{sample}.bam.bai"
  shell:
    "samtools index {input}"
  