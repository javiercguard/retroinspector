configfile: "config.default.yaml"

import os, sys
from pprint import pformat

fastqDir = ''
if os.path.exists(config["outputPath"]+"/fastq"):
    fastqDir = config["outputPath"]+"/fastq" 
if "fastq_directory" in config.keys() and config["fastq_directory"] and not config["samples"]:
    fastqDir = config["fastq_directory"]
if fastqDir:
    samples = [x for x in os.listdir(fastqDir) if x.endswith(".fastq") or x.endswith(".fastq.gz")]
    for s in samples:
        config["samples"][s.split(".")[0]] = fastqDir + "/" + s

if config["referenceGenome"] == "":
    config["referenceGenome"] = str(workflow.basedir) + "/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
config["callerInfixes"] = [x.lower() for x in config["callers"]]

print(file = sys.stderr)
print(pformat(config["samples"], sort_dicts = False), file = sys.stderr)
print(file = sys.stderr)

include: "rules/alignment.smk"
include: "rules/variants.smk"
include: "rules/insertions.smk"
include: "rules/deletions.smk"
include: "rules/reference.smk"
include: "rules/r.smk"

workdir: config['outputPath']
config["minimumReadSupport"]: int(config["minimumReadSupport"])

def init():
    result = [
        f"variants/{config['allPrefix']}.te.vcf.gz",
        f"variants/{config['allPrefix']}.te.lax.vcf.gz",
    ]
    if config["mode"] == "full":
        result += [f"reports/{x}_vs_{y}.html" for x, y in config["comparisons"]] + \
        [
            f"reports/report.{config['allPrefix']}.html",
        ]
    return result

rule all:
    input: init()
