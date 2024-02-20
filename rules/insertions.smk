rule assemblyAlleles:
  conda: "../env.yaml"
  log:
    "logs/polish/{sample}.{aligner}.log"
  threads: config['threads']
  input:
    script = str(workflow.basedir) + "/scripts/getGoodAlts.py",
    vcf = "variants/{aligner}/{sample}.{aligner}.vcf.gz",
    vcfIndex = "variants/{aligner}/{sample}.{aligner}.vcf.gz.csi",
    bam = "alns/{sample}.bam",
    bai = "alns/{sample}.bam.bai"
  output: 
    temp("tmp/{sample}.{aligner}.polished.vcf"), 
    "variants/{aligner}/{sample}.{aligner}.polished.vcf.gz",
    "variants/{aligner}/{sample}.{aligner}.polished.vcf.gz.csi"
  shell:
    """
    python3 {input.script} \
      {input.vcf} \
      {input.bam} \
      {output[0]} \
      {threads};
    bcftools sort {output[0]} -O z -o {output[1]} 2>> {log};
    bcftools index {output[1]} 2>> {log}
    """

rule merge_insertions_intrapatient:
  conda: "../env.yaml"
  log:
    "logs/merge_intrapatient/{sample}.log"
  input:
    script = str(workflow.basedir) + "/scripts/merge.py",
    files = [
      "variants/cutesv/{sample}.cutesv.polished.vcf.gz",
      "variants/cutesv/{sample}.cutesv.polished.vcf.gz.csi",
      "variants/svim/{sample}.svim.polished.vcf.gz",
      "variants/svim/{sample}.svim.polished.vcf.gz.csi"
      ]
  output: 
    temp("tmp/{sample}.merged.both.vcf"),
    "variants/{sample}.merged.both.vcf.gz",
    "variants/{sample}.merged.both.vcf.gz.csi"
  params:
    distanceLimit = config["insertionDistanceLimitIntraPatient"]
  shell:
    """
    python3 {input.script} \
      -samples cuteSV SVIM \
      -vcf {input.files[0]} {input.files[2]} \
      -o {output[0]} \
      -d {params.distanceLimit} 2> {log};
    bcftools sort {output[0]} -O z -o {output[1]} 2>> {log};
    bcftools index {output[1]} 2>> {log}
    """

def getIntraPatientMerges (wilcards):
  return [f"variants/{sample}.merged.both.vcf.gz" for sample in config["samples"]]
def getIntraPatientMergesIndexes ():
  return [f"variants/{sample}.merged.both.vcf.gz.csi" for sample in config["samples"]]
rule merge_insertions_interpatient:
  conda: "../env.yaml"
  log:
    "logs/merge_interpatient.log"
  input:
    script = str(workflow.basedir) + "/scripts/merge.py",
    files = getIntraPatientMerges,
    indexes = getIntraPatientMergesIndexes()
  output: 
    temp(f"tmp/{config['allPrefix']}.merged.vcf"), 
    f"variants/{config['allPrefix']}.merged.vcf.gz",
    f"variants/{config['allPrefix']}.merged.vcf.gz.csi"
  params:
    distanceLimit = config["insertionDistanceLimitInterPatient"],
    samples = list(config["samples"].keys())
  shell:
    """
    python3 {input.script} \
      -samples {params.samples} \
      -vcf {input.files} \
      -o {output[0]} \
      -d {params.distanceLimit} 2> {log};
    bcftools sort {output[0]} -O z -o {output[1]} 2>> {log};
    bcftools index {output[1]} 2>> {log}
    """

rule repeatmasker:
  conda: "../env.yaml"
  log:
    "logs/repeatmasker.log"
  input:
    script = str(workflow.basedir) + "/scripts/checkInsertionsMultiSample.py",
    file = f"variants/{config['allPrefix']}.merged.vcf.gz",
    index = f"variants/{config['allPrefix']}.merged.vcf.gz.csi"
  params:
    direct = "repeatmasker",
    species = config["species"],
    threads = config['threads']
  output: 
    f"repeatmasker/{config['allPrefix']}.merged.rm.unproc.bed",
    # cannot reference params.direct, so it is harcoded
  shell:
    """
    mkdir -p {params.direct} ;
    python3 {input.script} \
      {input.file} \
      {params.direct} \
      {params.species} \
      {params.threads}
    """
rule svaf:
  conda: "../env.yaml"
  input: f"repeatmasker/{config['allPrefix']}.merged.rm.unproc.bed"
  output: 
    bed = f"repeatmasker/{config['allPrefix']}.merged.rm.bed",
    fasta = temp("tmp/svaf.fa"),
    db = temp(expand("tmp/svaf.db.{ext}", ext = ["ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"])),
    blast = temp("tmp/svaf.blast6"),
  params:
    dbPrefix = "tmp/svaf.db",
    svafTypes = str(workflow.basedir) + "/data/sva_types.fa",
  script: "../scripts/svaf.py"
