def conditionalTemp(what, infix = "rds"):
  return what if config["keepRds"] else temp(f"/tmp/{what}")

rule filter_rmbed:
  input:
    f"repeatmasker/{config['allPrefix']}.merged.rm.bed"
  output:
    temp(f"tmp/{config['allPrefix']}.merged.rm.bed")
  shell:
    """
    grep -E "^chr([[:digit:]]{{1,2}}|[XY])\t" {input[0]} > {output[0]} # Not useful for other species : (
    """

rule analysis_preparatory:
  conda: "../r.yaml"
  log:
    "logs/r_analysis.log"
  input:
    repeatMaskerBed = f"tmp/{config['allPrefix']}.merged.rm.bed",
    repeatMaskerVcfDel = f"variants/{config['allPrefix']}.me.deletions.vcf.gz",
  output:
    annotation = f"tmp/rds/annotation.{config['allPrefix']}.rds", #! make temp
    annotatedInsertionsMin3 = f"tmp/annotatedInsertionsMin3.unfiltered.{config['allPrefix']}.rds", #! make temp
    insertionsTable = f"rds/insertionsTable.unfiltered.{config['allPrefix']}.rds", #! make temp
    allIns = conditionalTemp(f"rds/allIns.rds"),
    # insertionsBed = temp(f"tmp/insertionsBed.{config['allPrefix']}.bed"),
    # insertionsBedRds = temp(f"insertionsBedRds.{config['allPrefix']}.rds"),
    meDeletionsMin3 = f"tmp/rds/meDeletionsMin3.unfiltered.{config['allPrefix']}.rds", #! make temp
    # deletionsBed = temp(f"tmp/deletionsBed.{config['allPrefix']}.bed"),
    # deletionsBedRds = temp(f"rds/deletionsBedRds.{config['allPrefix']}.rds"),
  params:
    samples = list(config["samples"].keys()),
  script:
    "../scripts/analysisPreparatory.R"

rule analysis_genotyping:
  conda: "../r.yaml"
  log:
    "logs/r_analysis_genotyping.log"
  input:
    # insertionsBedRds = rules.analysis_preparatory.output.insertionsBedRds,
    # deletionsBedRds = rules.analysis_preparatory.output.deletionsBedRds,
    # insertionBedFiles = [f"mosdepth/{sample}.ins.regions.bed.gz" for sample in config["samples"]],
    insertionsTable = rules.analysis_preparatory.output.insertionsTable,
    annotatedInsertionsMin3 = rules.analysis_preparatory.output.annotatedInsertionsMin3,
    # cppFile = str(workflow.basedir) + "/scripts/insSupport.cpp",
    # repeatsReferenceTENoOverlap = "tmp/repeatsReferenceTENoOverlap.bed",
    # mosdepthDelFiles = [f"mosdepth/{sample}.del.regions.bed.gz" for sample in config["samples"]],
    # repeatMaskerVcfDel = f"variants/{config['allPrefix']}.me.deletions.vcf.gz",
    # meDeletionsMin3 = rules.analysis_preparatory.output.meDeletionsMin3,
    # vcfs = [f"variants/survivor/{sample}.merged.survivor.vcf.gz" for sample in config["samples"]],
  output:
    # insertionsTable = conditionalTemp(f"rds/insertionsTable.{config['allPrefix']}.rds"),
    # annotatedInsertionsMin3 = conditionalTemp(f"rds/annotatedInsertionsMin3.{config['allPrefix']}.rds"),
    # meDeletionsMin3 = conditionalTemp(f"rds/meDeletionsMin3.{config['allPrefix']}.rds"),
    # observedMafTable = conditionalTemp(f"rds/observedMafTable.{config['allPrefix']}.rds"),
    genes = temp(f"rds/genes.{config['allPrefix']}.rds"),
    # genotypedDeletions = conditionalTemp(f"rds/genotypedDeletions.{config['allPrefix']}.rds"),

    vcfBody = f"tmp/{config['allPrefix']}.me.insertions.txt", #! make temp
    vcfBodyLax = f"tmp/{config['allPrefix']}.me.insertions.lax.txt", #! make temp
  threads: 16
  params:
    samples = list(config["samples"].keys()),
  script:
    "../scripts/analysisGenotyping.R"

rule analysis_enrichment:
  conda: "../r.yaml"
  input:
    dt = rules.analysis_genotyping.output.genes,
  output:
    egoMF = conditionalTemp("rds/egoMF.rds"),
    egoBP = conditionalTemp("rds/egoBP.rds"),
    egoCC = conditionalTemp("rds/egoCC.rds"),
    # kegg = conditionalTemp("rds/kegg.rds"),
    do = conditionalTemp("rds/do.rds"),
    ncg = conditionalTemp("rds/ncg.rds"),
  params:
    samples = list(config["samples"].keys()),
    plimit = config["enrichmentSignificanceThreshold"],
  script:
    "../scripts/enrichment.R"

def getPreparedDatasets():
  return [f"variants/{config['allPrefix']}.vs.{dataset}.vcf" for dataset in config["datasets"]]
rule generate_report:
  conda: "../r.yaml"
  log:
    "logs/r_report.log"
  input:
    insertionsTable = rules.analysis_preparatory.output.insertionsTable,
    allIns = rules.analysis_preparatory.output.allIns,
    annotatedInsertionsMin3 = rules.analysis_preparatory.output.annotatedInsertionsMin3,
    # observedMafTable = rules.analysis_genotyping.output.observedMafTable,
    meDeletionsMin3 = rules.analysis_preparatory.output.meDeletionsMin3,
    # genotypedDeletions = rules.analysis_genotyping.output.genotypedDeletions,
    # enrichment
    egoMF = rules.analysis_enrichment.output.egoMF,
    egoBP = rules.analysis_enrichment.output.egoBP,
    egoCC = rules.analysis_enrichment.output.egoCC,
    # kegg = rules.analysis_enrichment.output.kegg,
    do = rules.analysis_enrichment.output.do,
    ncg = rules.analysis_enrichment.output.ncg,
    # other studies
    otherSets = getPreparedDatasets(),
    # auxiliary scripts, just used to rerun the rule when they are updated
    auxScripts = str(workflow.basedir) + "/scripts/sets.Rmd",
  output:
    report = f"reports/report.{config['allPrefix']}.html"
  threads: 16
  params:
    samples = list(config["samples"].keys()),
    plimit = config["enrichmentSignificanceThreshold"],
    ownSet = config["allPrefix"],
    out = config["outputPath"],
  script:
    "../scripts/report.Rmd"

rule compare:
  conda: "../r.yaml"
  log:
    "logs/r_compare_{sample1}_vs_{sample2}.log"
  input:
    insertionsTable = rules.analysis_preparatory.output.insertionsTable,
    annotatedInsertionsMin3 = rules.analysis_preparatory.output.annotatedInsertionsMin3,
  output:
    "reports/{sample1}_vs_{sample2}.html"
  params:
    soi = ["{sample1}", "{sample2}"],
    samples = config["samples"].keys()
  script:
    "../scripts/comparison.Rmd"

rule generate_vcf:
  conda: "../env.yaml"
  input:
    vcfBody = rules.analysis_genotyping.output.vcfBody,
    vcfBodyLax = rules.analysis_genotyping.output.vcfBodyLax,
    header = str(workflow.basedir) + "/data/header.txt",
  output: 
    f"variants/{config['allPrefix']}.te.vcf.gz",
    f"variants/{config['allPrefix']}.te.lax.vcf.gz",
  shell:
    """
    cat {input.header} {input.vcfBody} | bcftools sort -O z -o {output[0]} ;
    bcftools index {output[0]}
    cat {input.header} {input.vcfBodyLax} | bcftools sort -O z -o {output[1]} ;
    bcftools index {output[1]}
    """

rule download_indigen:
  output: "data/nstd215.GRCh38.variant_call.vcf.gz"
  shell: 
    """
    wget -O {output[0]} \
    https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd215.GRCh38.variant_call.vcf.gz
    """
rule prepare_indigen:
  conda: "../env.yaml"
  input: rules.download_indigen.output
  output: "data/indigen.prep.vcf.gz", "data/indigen.prep.vcf.gz.csi"
  params:
    vcf = "data/indigen.prep.vcf", 
  shell:
    """
    zcat {input[0]} | awk -v'OFS=\\t' \
      '{{if ($0 ~ /^#/) {{print $0}} else {{if ($1 ~ /^([[:digit:]]{{1,2}}|[XYM])/) {{$1 = "chr"$1; print $0}} }} }}' | \
      awk  -v'OFS=\\t' \
      '{{if ($0 ~ /^#/) {{print $0}} else {{$5="<INS>"; print $0}}}}' > \
      {params.vcf}
    bgzip -f "{params.vcf}"
    tabix -f --csi "{output[0]}"
    """

rule mergeSet:
  conda: "../env.yaml"
  input:
    datafile= "data/{dataset}.prep.vcf.gz",
    vcf= f"variants/{config['allPrefix']}.te.vcf.gz",
    script = str(workflow.basedir) + "/scripts/merge.py",
  output:
    vcf= f"variants/{config['allPrefix']}.vs.{{dataset}}.vcf",
  params:
    distanceLimit = config["insertionDistanceLimitInterPatient"],
    ownSet = config["allPrefix"]
  shell:
    """
    python3 {input.script}  \
      -samples {params.ownSet} {wildcards.dataset} \
      -vcf {input.vcf} {input.datafile} \
      -o {output.vcf} \
      -d {params.distanceLimit};
    """

# rule testSets:
  # conda: "../r.yaml"
  # input:
  #   **dict(rules.generate_report.input),
  #   # indigen = "data/indigen.prep.vcf.gz",
  # params:
  #   **dict(rules.generate_report.params),
  # output:
  #   "reports/sets.html"
  # script: "../scripts/sets.Rmd"
