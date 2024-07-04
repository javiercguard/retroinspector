rule run_survivor_intrasample:
  conda: "../env.yaml"
  log:
    "logs/survivor/{sample}.log"
  input:
    vcfs = expand("variants/{infix}/{{sample}}.{infix}.vcf.gz", infix = config["callerInfixes"])[::-1],
    indexes = expand("variants/{infix}/{{sample}}.{infix}.vcf.gz.csi", infix = config["callerInfixes"]),
  output: 
    vcf = temp("tmp/{sample}.merged.survivor.ungenotyped.vcf"),
  params:
    distance = config["survivorInsertionDistanceLimitIntraPatient"],
  shell:
    """
    surpyvor merge --variants {input.vcfs} \
      -d {params.distance} \
      -c 1 \
      -l 1 \
      -o {output.vcf} 2> {log}; \
    """

rule create_mosdepth_bed_del:
  conda: "../env.yaml"
  log:
    "logs/geno/{sample}.createBed.del.log"
  input:
    vcf = "tmp/{sample}.merged.survivor.ungenotyped.vcf",
  output:
    bed = temp("tmp/{sample}.del.bed"),
  params:
    svType = "DEL",
  script: "../scripts/vcfToBedForMosdepth.py"

rule mosdepth_for_genotyping_del:
  conda: "../env.yaml"
  log:
    "logs/mosdepth/{sample}.del.log"
  threads: 4
  input:
    mosdepthBed = "tmp/{sample}.del.bed",
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai",
  output:
    "mosdepth/{sample}.del.regions.bed.gz",
  params:
    prefix = "mosdepth/{sample}.del",
  shell:
    """
    mosdepth -t {threads} -Q 20 -n -b {input.mosdepthBed} \
    {params.prefix} {input.bam} &> {log}
    """

rule genotype_del:
  conda: "../env.yaml"
  log:
    "logs/geno/{sample}.geno.del.log"
  input:
    coverage = "mosdepth/{sample}.del.regions.bed.gz",
    vcf = "tmp/{sample}.merged.survivor.ungenotyped.vcf",
  output: 
    vcf = temp("variants/survivor/{sample}.merged.survivor.vcf"),
  params:
    svType = "DEL"
  script: "../scripts/genotype.py"

rule sort_index_vcf_del:
  conda: "../env.yaml"
  input: 
    "variants/survivor/{sample}.merged.survivor.vcf",
  output:
    "variants/survivor/{sample}.merged.survivor.vcf.gz",
    "variants/survivor/{sample}.merged.survivor.vcf.gz.csi",
  shell:
    """
    bcftools sort -O z -o {output[0]} {input[0]}
    bcftools index {output[0]}
    """

rule run_survivor_intersample:
  conda: "../env.yaml"
  log:
    "logs/survivor/all.log"
  input:
    vcfs = [
      f"variants/survivor/{sample}.merged.survivor.vcf.gz" 
        for sample in config["samples"].keys()
      ],
    index = [
      f"variants/survivor/{sample}.merged.survivor.vcf.gz.csi" 
        for sample in config["samples"].keys()
      ],
  output: 
    vcf = f"variants/survivor/{config['allPrefix']}.merged.survivor.vcf.gz",
    index = f"variants/survivor/{config['allPrefix']}.merged.survivor.vcf.gz.csi",
  params:
    distance = config["survivorInsertionDistanceLimitInterPatient"],
  shell:
    """
    surpyvor merge --variants {input.vcfs} \
      -d {params.distance} \
      -c 1 \
      -l 1 \
      -o {output.vcf} 2> {log}; \
    bcftools index {output.vcf} 2>> {log}
    """

rule get_deletions:
  conda: "../env.yaml"
  log:
    "logs/getDeletions.log"
  input:
    script = str(workflow.basedir) + "/scripts/getMeDeletions.py",
    vcf = f"variants/survivor/{config['allPrefix']}.merged.survivor.vcf.gz",
    index = f"variants/survivor/{config['allPrefix']}.merged.survivor.vcf.gz.csi",
    bed = "data/repeatsReferenceTE.bed.gz",
    csi = "data/repeatsReferenceTE.bed.gz.csi",
  output: 
    vcfTemp = temp(f"variants/medeletions/{config['allPrefix']}.me.deletions.vcf.gz"),
    vcf = f"variants/{config['allPrefix']}.me.deletions.vcf.gz",
    index = f"variants/{config['allPrefix']}.me.deletions.vcf.gz.csi",
  shell:
    """
    python {input.script} \
      {input.vcf} \
      {input.bed} \
      {output.vcfTemp} 2> {log};\
    bcftools sort {output.vcfTemp} -O z -o {output.vcf} 2>> {log}
    bcftools index {output.vcf} 2>> {log}
    """
