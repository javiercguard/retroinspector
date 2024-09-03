options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(ggplot2)
library(magrittr)
library(qqman) %>% suppressMessages()

# Patients and folders
# ----
samples = snakemake@params[["samples"]]
# ----

if (F) {
  # For debugging
  saveRDS(snakemake, paste0(snakemake@config[["outputPath"]], "/snk.rds"))
  q()
  snakemake = readRDS("~/path/snk2.rds")
  setwd("~/path/")
}
if (F) {
  snakemake = readRDS("~/ontology/proj/public5/snk.rds")
}

# Annotation is generated and exported
print("Building annotation")
if (T) {
  library(annotatr)
  library(AnnotationHub)
  library("TxDb.Hsapiens.UCSC.hg38.knownGene") # This package is necessary, but not automatically installed
  annotationNames = builtin_annotations()[
    stri_detect_fixed(builtin_annotations(), "hg38") # all hg38 are built, although not used
  ]
  annotationNames = annotationNames[annotationNames != "hg38_cpg_inter" &
      annotationNames != "hg38_enhancers_fantom"] # We remove the "regions that are not CpG islands" annotation
  # and enhancers_fantom
  annotation = build_annotations(genome = "hg38", annotations = annotationNames)
  annotation = annotation[GenomicRanges::mcols(annotation)$type != "hg38_cpg_inter"] # The "regions that are not CpG islands" has to be removed again
  saveRDS(annotation, snakemake@output[["annotation"]])
}
print("Built annotation")

# Do the annotation
print("Annotating")
if (T) {
  print(paste0("read regions con ", snakemake@input[["repeatMaskerBed"]]))
  # repeatMaskerTable = annotatr::read_regions( # repeatMaskerTable: results from RepeatMasker, as a data.table
  #   con = snakemake@input[["repeatMaskerBed"]],
  #   # con = "~/datos-javier/retropip_out/tmp/all.merged.rm.bed",
  #   format = "bed",
  #   genome = "hg38",
  #   rename_score = "meaningless", # This column is necessary in BED format, but it has no meaning here
  #   extraCols = c(
  #     seqId = "character",
  #     sw_score = "numeric",
  #     repeat.class = "character",
  #     repeat.subclass = "character",
  #     repeat.start = "numeric",
  #     repeat.end = "numeric",
  #     repeat.left = "numeric",
  #     repeat.strains = "character",
  #     repeat.divergence_percentage = "numeric",
  #     repeat.deletion_percentage = "numeric",
  #     repeat.insertion_percentage = "numeric",
  #     vcf_alt = "character",
  #     vcf_info = "character",
  #     setNames(
  #       rep("character", times = length(samples)),
  #       paste0("id", samples)
  #       )
  #   )
  # ) %>% as.data.table()
  repeatMaskerTable = fread(
    snakemake@input[["repeatMaskerBed"]],
    col.names = c(
    "seqnames", "start", "end", "name", "score", "strand", "seqId", 
    "sw_score", 
    # "rename_score", 
    "repeat.class",
    "repeat.subclass",
    "repeat.start",
    "repeat.end",
    "repeat.left",
    "repeat.strains",
    "repeat.divergence_percentage",
    "repeat.deletion_percentage",
    "repeat.insertion_percentage",
    "vcf_alt",
    "vcf_info",
    paste0("id", samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))
    )
  )
  print("Annotated")
  # Generating the usable datasets: min3 and trusty, lax and strict critera, respectively
  genoFields = c("GT", "PSV", "LN", "DR", "ST", "QV", "TY", "ID", "RAL", "AAL", "CO", "SC")
  genoFieldsOfInt = c("GT", "PSV", "DR") # Used for counting occurrence and genotyping
  pos = which(genoFields %in% genoFieldsOfInt) # The position of "PSV" and "DR" in the genotype columns

  for (col in repeatMaskerTable %>% colnames() %>%
       .[. %in% paste0("id",
                       samples %>% stri_replace_all_fixed(., pattern = "-", replacement = "."))]
       
  ) {
    
    # create columns for GT, PSV and DR for each patient
    repeatMaskerTable[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
    # keep second value of DR, which is the one supporting the INS
    # Variant callers dont consistently report the support for the reference allele
    repeatMaskerTable[, paste0(col, "_DR") := tstrsplit(.SD[[paste0(col, "_DR")]], ",", fixed = T, keep = 2)]
    # Fix NAs
    repeatMaskerTable[,
        paste0(col, "_DR") :=
          fifelse(.SD[[paste0(col, "_DR")]] == ".", 0, as.numeric(.SD[[paste0(col, "_DR")]]) )
        ] # This will warn about NA by coercion, which are not included in the data.table
  }

  # Take the DR columns, create a bool vector with DR >=3, make it numeric (1 or 0), assign as column
  # this creates a column with the lax criterion
  # stored as numeric for addition
  repeatMaskerTable[, SUPP_VEC_min3 := apply(.SD, 2, `>=`, y = 3) %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_DR$", colnames(repeatMaskerTable))]
  repeatMaskerTable[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")] # Count SUPPORT with the min3 criteria

  # Initialize the column with only the two callers criteria
  # Since a data.table is a list of columns (like a data.frame)
  # apply with margin = 2 iterates rows instead of columns, and the opposite for margin = 1
  # because it is "transposed"
  repeatMaskerTable[
    , SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>%
      apply(., 2, as.integer) %>%
      apply(., 1, paste0, collapse = ""),
    .SDcols = patterns("_PSV$")]
  # The value for the column is provisional, and missing the RE >= 3 part

  # This should be in R, but alas
  # ("00", "11") => "00"
  # ("01", "11") => "01"
  # ("11", "11") => "11"
  bitwiseAndStr = function (a, b) {
    lapply (1:length(a), function (i) {
      a2 = unlist( strsplit(a[[i]], split = "", fixed = T) )
      b2 = unlist( strsplit(b[[i]], split = "", fixed = T) )
      paste0(ifelse(a2 == b2 & a2 == "1", "1", "0"), collapse = "")
    }
    )
  }

  # Now we do both callers AND min3 -> strict criterion
  correctTrusty = unlist( bitwiseAndStr(repeatMaskerTable$SUPP_VEC_min3, repeatMaskerTable$SUPP_VEC_trusty) )
  repeatMaskerTable[, SUPP_VEC_trusty := correctTrusty]
  rm(correctTrusty)
  repeatMaskerTable[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]

  # We will also need the INS' length
  repeatMaskerTable[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", vcf_info))]

  # We extract the genotype from the columns and store it as numeric
  # 0 is hom. ref., 1 is het., 2 is hom. alt.
  for (col in repeatMaskerTable %>% colnames() %>% 
       grep(pattern = "_GT$", x = ., value = T)) {
    repeatMaskerTable[, (col) := get(col) %>%
                                   lapply(function (x) {
                                     x %>%
                                       strsplit(split = "/", fixed = T) %>%
                                       unlist() %>%
                                       as.numeric() %>%
                                       sum()
                                   }) %>% unlist()]
  }
  repeatMaskerTable[
    , genotype := .SD %>% 
      simplify2array() %>% 
      t() %>% 
      split(., rep(1:ncol(.), each = nrow(.))) %>% 
      unname(), 
    .SDcols = grep(pattern = "_GT$", colnames(repeatMaskerTable))]
  repeatMaskerTable[, SUPP_mask_geno := genotype %>% lapply(`>=`, y = 1)]
  repeatMaskerTable[, SUPP_VEC_geno := SUPP_mask_geno %>% lapply(as.numeric) %>% lapply(paste0, collapse = "") %>% unlist()]
  repeatMaskerTable[, SUPP_geno := SUPP_mask_geno %>% lapply(sum) %>% unlist()]
  repeatMaskerTable[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2 ) ) %>% unlist() ] # humans are diploid, thus the length() * 2
  # A vector column of supports from callers, for genotyping
  # supports = repeatMaskerTable[, .SD, .SDcols = grep("DR$", colnames(repeatMaskerTable))]
  # supportColumn = supports %>% as.matrix() %>% t() %>% split(rep(1:ncol(.), each = nrow(.)))
  # repeatMaskerTable[, support := supportColumn]
  # rm(supports, supportColumn) ; gc()

  # Create vector bool columns
  repeatMaskerTable[, SUPP_mask_min3 := SUPP_VEC_min3 %>%
              lapply(strsplit, split = "") %>%
              lapply(function(x)as.logical(as.numeric(unlist(x))))]
  repeatMaskerTable[, SUPP_mask_trusty := SUPP_VEC_trusty %>%
              lapply(strsplit, split = "") %>%
              lapply(function(x)as.logical(as.numeric(unlist(x))))]

  repeatMaskerTable[, # split level2-level3 and replace NA for level 3 with "-"
    c("repeat.subclass", "repeat.family") := stri_match_all_regex(repeat.subclass, "([^-]+)(?:-([^-]+))?$") %>% rapply(., f = `[`, ... = 2:3, how = "r") %>%
        lapply(function (x) {x[is.na(x)] = "-"; x}) %>% data.table::transpose() ]

  # Remove unnecessary columns
  repeatMaskerTable[, grep("^id", colnames(repeatMaskerTable)) := NULL]
  repeatMaskerTable[, vcf_info := NULL]
  # annotatr::annotate_regions(regions = GenomicRanges::makeGRangesFromDataFrame(
  #       repeatMaskerTable,
  #       seqnames.field = "seqnames",
  #       keep.extra.columns = T),
  #                      annotations = annotation,
  #                      ignore.strand = F,
  #                      quiet = T
  #     )
  # Now lets annotate
  annotatedInsertions = data.table(data.frame(
      annotatr::annotate_regions(regions = GenomicRanges::makeGRangesFromDataFrame(
        repeatMaskerTable,
        seqnames.field = "seqnames",
        keep.extra.columns = T),
                       annotations = annotation,
                       ignore.strand = F,
                       quiet = T
      )
    ))
  rm(repeatMaskerTable)
  print("Annotated")
  annotatedInsertions[, annot.type := stri_replace_first_regex(annot.type, "hg38_(.*)", "$1")] # "hg38_genes_introns" => "genes_introns"
  # Lets separate "genes_introns" into 2 columns
  text = stri_match_all_regex(annotatedInsertions$annot.type, "([^-]*)_(.*)") %>%
    lapply(., function(x) {x[2:3]})
  annotatedInsertions[, `:=`(annot.class = text %>% lapply(`[`, 1) %>% unlist(),
                 annot.subclass = text %>% lapply(`[`, 2) %>% unlist()
                 )]
  annotatedInsertions[, annot.type := NULL] # Redundant
  rm(text)
  print("Added class and subclass")

  annotatedInsertions[, c("annot.strand", "annot.seqnames", "annot.width", "annot.id") := NULL]
  # No saveRDS yet, because
  # we only save the ones that pass at least the min3 criteria for at least one occurrence
}

rm(annotation);gc()
print("Finished annotation")

#
# Genotyping
#
# ----
chrs = paste0("chr", c(1:22, "X", "Y")) # Human chromosomes used to filter out Unknown and Random fragments
# chrM and Epstein-Barr are removed, too, but they had no presence in our study
class1 = c("Retroposon", "SINE", "LINE", "LTR") # Families of class I elements
class2 = c("DNA") # Class II elements

annotatedInsertionsMin3 = annotatedInsertions[SUPP_min3 > 0]
# if (F) {
#   repeatMaskerTable[SUPP_min3 > 0] %>% nrow()
#   annotatedInsertions[, .SD[1], by = seqId] %>% nrow()
# }
allIns = annotatedInsertions[ # An unfiltered version, with only one row per INS (so no annotation)
  , repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN
][, .SD[1], by = seqId]
print("Created <allIns> table")
rm(annotatedInsertions); gc()
annotatedInsertionsMin3[, repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN] # Percentage of INS covered by annotation
print("Added <repeat.percentage> column")
# Added percentage column
# We filter out the ones that are not >=85% repetitive elements (includes TEs and non TEs),
# the ones that are not TEs
# and the ones not in chr1-22, X or Y
annotatedInsertionsMin3 = annotatedInsertionsMin3[repeat.percentage >= 0.85 & repeat.class %in% c(class1, class2) & seqnames %in% chrs]
print("Filtered annotatedInsertionsMin3")
insertionsTable = annotatedInsertionsMin3[, .SD[1], by = seqId] # Only one annotation per INS, a "unique" version of the above
print("Created <insertionsTable> table")
insertionsTable[, grep("^annot\\.", colnames(insertionsTable)) := NULL] # Remove the annotation data (incomplete in this table)
print("Removed annotation columns from <insertionsTable>")
saveRDS(annotatedInsertionsMin3, snakemake@output[["annotatedInsertionsMin3"]])
saveRDS(insertionsTable, snakemake@output[["insertionsTable"]])
saveRDS(allIns, snakemake@output[["allIns"]])

# 
# Deletions
# 
meDeletionsMin3 = fread(snakemake@input[["repeatMaskerVcfDel"]], skip = "#CHROM")
# This removes "LTR?" elements and similar ones
meDeletionsMin3 = meDeletionsMin3[!INFO %like% "ME=[^;]*\\?[^;]*" & !INFO %like% "MEFAM=[^;]*\\?[^;]*"]
meDeletionsMin3 = meDeletionsMin3[`#CHROM` %in% chrs]

meDeletionsMin3[ # Ids from survivor output are not guaranteed to be unique, and indeed they are not
  , `:=`(
    survivorId = ID,
    ID = paste0("surv", 1:nrow(meDeletionsMin3)),
    svlen = INFO %>% sub(".*SVLEN=-([0-9]+).*", "\\1", x = .) %>% as.numeric()
  )
]

# Add columns
# ----
# The read supports gets lost in the intersample merge, but since the file is min3,
# all instances are min3, and if PSV != NaN, RE >= 3 -> min3 for that instance, since it existed
# then, if PSV == "11" -> trusty, since min3 is already checked
genoFields = meDeletionsMin3[1, 9] %>% unname() %>% stri_split_fixed(pattern = ":") %>% unlist()
genoFieldsOfInt = c("GT", "PSV")
pos = which(genoFields %in% genoFieldsOfInt)
for (col in meDeletionsMin3 %>% colnames() %>% 
     grep(pattern = "^Sample", x = ., ignore.case = T, value = T) ) {
  meDeletionsMin3[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
}
meDeletionsMin3[, SUPP_VEC_min3 := apply(.SD, 2, `!=`, y = "NaN") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_PSV$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")]
meDeletionsMin3[, SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = patterns("_PSV$")]
meDeletionsMin3[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
meDeletionsMin3[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", INFO)) %>% abs()]
# extract genotype for deletions
for (col in meDeletionsMin3 %>% colnames() %>% 
     grep(pattern = "_GT$", x = ., value = T)) {
  meDeletionsMin3[, (col) := get(col) %>% ifelse(. == "./.", "0/0", .) %>%
                    lapply(function (x) {
                      x %>%
                        strsplit(split = "/", fixed = T) %>%
                        unlist() %>%
                        as.numeric() %>%
                        sum()
                    }) %>% unlist()]
}
meDeletionsMin3[
  , genotype := .SD %>% 
    simplify2array() %>% 
    t() %>% 
    split(., rep(1:ncol(.), each = nrow(.))) %>% 
    unname(), 
  .SDcols = grep(pattern = "_GT$", colnames(meDeletionsMin3))]
meDeletionsMin3[, SUPP_VEC_geno := genotype %>% lapply(`>=`, y = 1)]
meDeletionsMin3[, SUPP_geno := SUPP_VEC_geno %>% lapply(sum) %>% unlist()]
meDeletionsMin3[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2) ) %>% unlist() ] # humans are diploid, thus the length() * 2

meDeletionsMin3[, grep("^(Sample|FORMAT)", colnames(meDeletionsMin3), ignore.case = T) := NULL]
meDeletionsMin3 = meDeletionsMin3[
  , .(
    seqnames = `#CHROM`,
    start = POS,
    end = sub(pattern = ".*;END=([^;]+).*", replacement = "\\1", x = INFO, perl = T) %>% as.numeric(),
    id = ID,
    SUPP_min3,
    SUPP_VEC_min3,
    SUPP_trusty,
    SUPP_VEC_trusty,
    genotype,
    maf,
    SUPP_geno,
    SUPP_VEC_geno,
    svlen,
    # info = INFO,
    name = sub(pattern = ".*ME=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.class = sub(pattern = ".*MEFAM=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.coords = sub(pattern = ".*MECOORDS=([^;]+).*", replacement = "\\1", x = INFO, perl = T)
  )
][, c("repeat.class", "repeat.subclass") := repeat.class %>% 
    stri_match_all_regex("([^/]+)(?:/(.+))?$") %>% 
    rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()][
      , repeat.subclass := fifelse(is.na(repeat.subclass), "-", repeat.subclass)
    ][
      , c("repeat.subclass", "repeat.family") := 
        stri_match_all_regex(repeat.subclass, "([^-]+|-)(?:-(.+))?$") %>% 
        rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()
    ][, repeat.family := fifelse(is.na(repeat.family), "-", repeat.family)]

# print("Annotating deletions") Not really interesting
# annotation = readRDS(snakemake@output[["annotation"]])
# # Annotation for deletions
# annotatedDeletionsMin3 = annotatr::annotate_regions(
#   regions = GenomicRanges::makeGRangesFromDataFrame(
#   meDeletionsMin3, 
#   seqnames.field = "seqnames", 
#   keep.extra.columns = T),
#   annotations = annotation,
#   ignore.strand = F,
#   quiet = T
# )
# annotatedDeletionsMin3 = annotatedDeletionsMin3 %>% as.data.frame(stringsAsFactors = F) %>% data.table()
# annotatedDeletionsMin3[, annot.type := stri_replace_first_regex(annot.type, "hg38_(.*)", "$1")]
# annotatedDeletionsMin3 = annotatedDeletionsMin3[stri_detect_fixed(annot.type, pattern = "genes_")]
# text = stri_match_all_regex(annotatedDeletionsMin3$annot.type, "([^-]*)_(.*)") %>% 
#   lapply(., function(x) {x[2:3]})
# annotatedDeletionsMin3[, `:=`(annot.class = text %>% lapply(`[`, 1) %>% unlist(),
#                annot.subclass = text %>% lapply(`[`, 2) %>% unlist()
# )]
# annotatedDeletionsMin3[, annot.type := NULL] # Redundant
# annotatedDeletionsMin3[, repeat.class := as.factor(repeat.class)]
# rm(text); gc()
# annotatedDeletionsMin3[, c("annot.strand", "annot.seqnames", "annot.width", "annot.id", "annot.tx_id") := NULL]

# deletionsBed = meDeletionsMin3[, c("seqnames", "start", "end", "id")]
# fwrite(deletionsBed, file = snakemake@output[["deletionsBed"]], sep = "\t", col.names = F, scipen=50)
# saveRDS(deletionsBed, snakemake@output[["deletionsBedRds"]])
saveRDS(meDeletionsMin3, snakemake@output[["meDeletionsMin3"]])
# saveRDS(annotatedDeletionsMin3, snakemake@output[["annotatedDeletionsMin3"]])
