rule run_survivor_intrasample:
  conda: "../env.yaml"
  log:
    "logs/survivor/{sample}.log"
  input:
    vcfs = [
      "variants/svim/{sample}.svim.vcf.gz",
      "variants/cutesv/{sample}.cutesv.vcf.gz",
      ],
    indexes = [
      "variants/svim/{sample}.svim.vcf.gz.csi",
      "variants/cutesv/{sample}.cutesv.vcf.gz.csi",
      ],
  output: 
    vcf = "variants/survivor/{sample}.merged.survivor.vcf.gz",
    index = "variants/survivor/{sample}.merged.survivor.vcf.gz.csi",
  params:
    distance = config["survivorInsertionDistanceLimitIntraPatient"],
  shell:
    """
    surpyvor merge --variants {input.vcfs} \
      -d {params.distance} \
      -c 1 \
      -l 1 \
      -o {output.vcf} 2> {log}; \
    bcftools index {output.vcf} 2>> {log}
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
