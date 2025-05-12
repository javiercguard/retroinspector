RetroInspector is a pipeline to detect and annotate TE insertions and deletions on human genomes sequenced by nanopore sequencing. It has been published on _[Advanced analysis of retrotransposon variation in the human genome with nanopore sequencing using RetroInspector](https://www.nature.com/articles/s41598-025-98847-7)_.

Please cite the article if you use it.

# Usage

## Installation

First, you will need to have `conda` or `mamba` on your computer. `mamba` is faster, but if you decide for `conda`, just replace it on the following commands.

Once you have `mamba`, you will need `snakemake`. You can install it on a conda/mamba environment like this:

```
mamba create -n snakemake snakemake
```

After that, clone this repository:

```
git clone https://github.com/javiercguard/retroinspector
```

That's it.

## Running

You need the hg38 genome, your FASTA files (one for each sample), and to choose where the output will go. 

We can check if everything is fine with the '-n' option, that will not run anything:

```
mamba activate snakemake
snakemake --use-conda -np --config \
    referenceGenome=<path to genome>
    outputPath=<outpath path> \
    fastq_directory=<input path>
```

__Note__: if you are using `conda`, you will need to add `--conda-frontend conda` (between `--use-conda` and `-np`, for example). This will be true for every command shown here, but will not be mentioned again.

Let's suppose you have saved the genome on `/home/yourname/hg38.fa`, you want to store the results at `/home/yourname/experiment_results`, and have the FASTA files at `/home/yourname/data/samples`. The command would be: 

```
snakemake --use-conda -np --config \
    referenceGenome="/home/yourname/hg38.fa"
    outputPath="/home/yourname/experiment_results" \
    fastq_directory="/home/yourname/data/samples"
```

After that, a list of jobs to run will be generated. If everything looks fine, remove `-n` to run them, although it's probably worth it to check the Configuration section down below in case you find any interesting options for your use case.

## Configuration

First, specifying a number of threads for `snakemake` to use is reccomended. This is done with `-c <number of threads>`.

If we are going to specify a few options, it will be easier to copy the template configuration file than writing all of them on the command line. 

```cp config.template.yaml config.<chosen name>.yaml```

To use the custom configuration file:

```
snakemake --use-conda -p \
    --config-file config.<chosen name>.yaml \
    -c 32
```

These are the options available, in __bold__ the ones more likely to need to be changed.

| Option | Meaning | Default value |
| --- | --- | --- |
| __`referenceGenome`__ | Full path to the reference genome | Example path |
| __`fastq_directory`__ | Full path to directory with FASTA files | Example path |
| __`outputPath`__ | Full path in which to store results | Example path |
| _`allPrefix`_ | An infix to put in the name of some files. For example, the report will be called `report.<allPrefix>.html` | retro |
| _`minimumReadSupport`_ | Minimum number of reads for variants to be called. | 3 |
| `insertionDistanceLimitIntraPatient` | When merging insertions from the same sample, maximum distance of coordinates for a group of insertions to be considered the same on.e | 20 |
| `insertionDistanceLimitInterPatient` | Same as above, but for merging from different samples. | 20 |
| `survivorInsertionDistanceLimitIntraPatient` | Same as above, but for SURVIVOR. | 10 |
| `survivorInsertionDistanceLimitInterPatient` | Same as above, but for SURVIVOR. | 500 |
| `keepRds` | Whether to keep R rds files, in case you want to load them latter. Use `True` or `False`, with the first letter capitalized. | `True` |
| `samples` | List of samples and their FASTQ, see below for more details. Leave as is for autodetection. See `config.example.yaml` for syntax. | `{}` |
| `comparisons` | List of samples to compare.  See `config.example.yaml` for syntax. | `[]` |

## Output

For each run, RetroInspector will generate a series of files. This is how it would look for a run with samples named S1 (S2, S3, etc.), with index files omitted (`.csi`, `.bai`, etc.):
```
├── alns
│   ├── S1.bam
│   └── S2.bam
├── logs
│   └── ... (logs for each step and sample)
├── rds
│   └── ... See tables.md for an exhaustive description.
├── repeatmasker
│   └── Standard RepeatMasker's output
├── reports
│   ├── S1_vs_S2.html
│   └── report.snakemake.html
└── variants
    ├── cutesv
    │   ├── S1.cutesv.polished.vcf.gz
    │   ├── S1.cutesv.vcf.gz
    │   └── ...
    ├── S1.merged.both.vcf.gz
    ├── ...
    ├── retro.me.deletions.vcf.gz
    ├── retro.merged.vcf.gz
    ├── public.te.lax.vcf.gz
    ├── public.te.vcf.gz
    ├── survivor
    │   ├── S1.merged.survivor.vcf.gz
    │   ├── ...
    │   └── retro.merged.survivor.vcf.gz
    └── svim
        ├── S1.svim.polished.vcf.gz
        ├── S1.svim.vcf.gz
        └── ...
```

- The alignments are in the `alns` directory.
- R files are in `rds`. See [tables.md](./tables.md) for a full description.
- Reports are inside `reports`.
- Inside `variants`:
    - Variants called by `cuteSV` are inside `cutesv`. The `polished` infix refers to the VCF that contains the reassembled inserted sequences.
    - Same for `svim`.
    - `S1.merged.both.vcf.gz` contains the union of the insertions called by both callers.
    - `retro.me.deletions.vcf.gz` contains all TE deletions (across samples) with allele frequency.
    - `public.te.vcf.gz` and `public.te.lax.vcf.gz` contain all TE insertions (across samples) for the strict and lax criterion, respectively.
    - `retro.merged.vcf.gz` contains all insertions merged across samples.
    - The directory `survivor` contains the result of merging both callers' results for each patient (`S1.merged.survivor.vcf.gz`, etc.), and an inter-sample merge (`retro.merged.survivor.vcf.gz`).

## Example

Please see [test_dataset/test.md](./test_dataset/test.md) for an example to run RetroInspector with publicly available data.
