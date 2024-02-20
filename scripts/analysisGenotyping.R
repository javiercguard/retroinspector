options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

print("Loading libraries")
library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(magrittr)
library(Rcpp)
library(parallel)
print("Loaded libraries")

if (F) {
  # For debugging
  saveRDS(snakemake, "~/path/snk2.rds")
  q()
  snakemake = readRDS("~/path/snk2.rds")
}

chrs = paste0("chr", c(1:22, "X", "Y")) # Human chromosomes used to filter out Unknown and Random fragments
class1 = c("Retroposon", "SINE", "LINE", "LTR") # Families of class I elements
class2 = c("DNA") # Class II elements

insertionsBed = readRDS(snakemake@input[["insertionsBedRds"]])
deletionsBed = readRDS(snakemake@input[["deletionsBedRds"]])
insertionCoverageFiles = snakemake@input[["insertionBedFiles"]]
insertionsTable = readRDS(snakemake@input[["insertionsTable"]])
annotatedInsertionsMin3 = readRDS(snakemake@input[["annotatedInsertionsMin3"]])
samples = snakemake@params[["samples"]]

seqids = insertionsBed[order(seqId)]$seqId
print(paste0("files: ", insertionCoverageFiles))
coveragesForGenotyping = lapply(insertionCoverageFiles, function (x) {
  fread(x,
        col.names = c("chr", "start", "end", "seqId", "coverage"))[order(seqId)]$coverage
})
coverageForGenotyping = data.table(
  seqId = seqids, 
  coverage = lapply(asplit(simplify2array(coveragesForGenotyping), 1), as.numeric)
)

colnamesoi = c("seqnames", "start", "end", "seqId", grep("^SUPP", colnames(insertionsTable), value = T), "support")

probRatio = function (q1, q2) {
  if (q1 / q2 > 0) {
    log = log10(q1 / q2)
    if (is.infinite(log) || is.nan(log)) return (0)
    else return(log)
  }
  else return(0)
}

sourceCpp(snakemake@input[["cppFile"]])

genotype = function (supportVal, coverageVal, 
                     chances = c(0.05, 0.5, 0.95), 
                     genotypes = c(0, 1, 2), 
                     genoNames = c("homRef", "het", "homAlt"),
                     maxNorm = 250) {
  if(supportVal > coverageVal) {
    coverageVal = supportVal 
  }
  maxValue = max(supportVal, coverageVal)
  if (maxValue > maxNorm) {
    factorNorm = 250 / maxValue
    supportVal = supportVal * factorNorm
    coverageVal = coverageVal * factorNorm
  }
  
  probs = lapply(chances, dbinomSnifflesCpp, x = supportVal, size = coverageVal) %>% unlist() %>% setNames(genoNames) %>% sort(decreasing = T)
  probs = probs / sum(probs)
  q1 = probs[1]
  q2 = probs[2]
  qual = min(60, floor(-10 * probRatio(q2, q1)))
  geno = q1 %>% attr("name")
  genoResult = (if (geno == "homRef") 0 else if (geno == "het") 1 else 2)
  c(genoResult, qual)
}
genotypeVec = Vectorize(genotype, vectorize.args = c("supportVal", "coverageVal"))

print("Finished initializing analysis_genotyping")
print("Genotyping")

genotypes = mclapply(1:length(samples), function (i) {
  x = samples[i]
  coverageFile = insertionCoverageFiles[i]
  support = unlist(lapply(insertionsTable[seqnames %in% chrs][order(seqId)]$support, `[`, ... = i))
  coverage = round(fread(coverageFile,
                         col.names = c("chr", "start", "end", "seqId", "coverage")
                         )[order(seqId)]$coverage)
  result = t(genotypeVec(support, coverage))
  result
}, mc.cores = snakemake@threads)

genotypesValues = mclapply(1:length(samples), function (i) {genotypes[[i]][, 1]}, mc.cores = snakemake@threads)
genotypesColumn = genotypesValues %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
if (seqids %>% unique() %>% length() != seqids %>% length()) {warning("Duplicated IDs in the genotyping!")}
genotyped = data.table(seqId = seqids, genotype = genotypesColumn)
genotyped = genotyped[insertionsTable, on = "seqId", nomatch = NULL]

genotyped = genotyped[repeat.percentage >= 0.85 &
                        repeat.class %in% c(class1, class2)]
genotyped[, SUPP_geno := genotype %>% lapply(`>=`, y = 1) %>% lapply(as.numeric) %>% lapply(sum) %>% unlist()]

print("Finsihed genotyping of insertions")

observedMafTable = copy(genotyped)
observedMafTable[, maf := genotype %>% lapply(function(x) sum(x) / (length(samples) * 2 ) ) %>% unlist() ] # humans are diploid, thus the length() * 2
# ----

# 
# Genes
# 
# ----
# Lets prepare a list of the genes affected by our INS
print("Processing gene info")
genes = annotatedInsertionsMin3[SUPP_trusty > 0 & annot.subclass %in% c("cds", "introns", "promoters") & !is.na(annot.symbol), 
                .SD[1], 
                by = c("seqId", "annot.symbol")][order(annot.symbol)]
genes = genes[, .(
  uniques = seqId %>% unique() %>% length(),
  effectLax = sum(SUPP_min3),
  effectStrict = sum(SUPP_trusty)
), by = annot.symbol][order(-effectStrict)]

repetitiveRef = fread(
  snakemake@input[["repeatsReferenceTENoOverlap"]], 
  select = c(1:3,5,6), 
  col.names = c("seqnames", "start", "end", "name", "fam")
  )
setkey(repetitiveRef, seqnames, start, end)

print("Finished gene processing")

# ----

# This is also used latter, liberally
roi = annotatedInsertionsMin3[repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") & SUPP_trusty > 0 & repeat.percentage >= 0.85] # repeats of interest for the article

# Insertions' targets
# ----
# A smaller version of ioi
ioi2 = roi[, .SD[1], by = seqId][,
           .(seqId, seqnames, start, end, 
           name, repeat.class, repeat.subclass,
           SVLEN,
           SUPP_min3, SUPP_trusty)]
setkey(ioi2, seqnames, start, end)

# This is a very inefficient step
# which takes a significant amount of time
# to do something practically useless
{
print("Starting splitting families and subfamilies")
fams = repetitiveRef$fam %>%
  stri_split_fixed(",") %>% 
  mclapply(function (x) {
    sub(pattern = "([^/]+)/.*", replacement = "\\1", perl = T, x = x) %>% 
      unique() %>% 
      paste0(collapse = ",")
  }, mc.cores = snakemake@threads)
subfams =
  repetitiveRef$fam %>%
  stri_split_fixed(",") %>% 
  mclapply(function (x) {
    sub(pattern = "[^/]+/(.*)", replacement = "\\1", perl = T, x = x) %>% 
      unique() %>% 
      paste0(collapse = ",")
  }, mc.cores = snakemake@threads)
repetitiveRef[, c("fam", "subfam") := .(fams, subfams)]
print("Finished splitting families and subfamilies")
}

o = foverlaps(repetitiveRef, ioi2, nomatch = NULL)
o[, c("start", "end", "seqnames", "i.start", "i.end") := NULL]

o2 = o[, .(
  inRepeats = .SD %>% nrow(),
  sameFam = .SD[stri_detect_regex(str = fam, pattern = paste0("(^|,)",repeat.class,"(,|$)"))] %>% nrow(),
  sameSubfam = .SD[stri_detect_regex(
    str = subfam, 
    pattern = paste0("(^|,)", repeat.subclass,"(\\?|(-[^,]*))?(,|$)"))] %>% nrow(),
  sameElement = .SD[stri_detect_regex(str = i.name, pattern = paste0("(^|,)",name,"(,|$)"))] %>% nrow() # Since name can be "Name1,Name2":
), by = c("repeat.class", "repeat.subclass")
][
  ioi2[, .(total = .N), by = repeat.subclass], on = "repeat.subclass", nomatch = NULL # join to get total count
] %T>% setcolorder(
  c("repeat.class", "repeat.subclass", "total", "inRepeats", "sameFam", "sameSubfam", "sameElement")
) %>% .[
  , .(
    class = repeat.class,
    subclass = repeat.subclass,
    total,
    inRepeats = inRepeats,
    sameFam = sameFam,
    sameSubfam = sameSubfam,
    sameElement = sameElement
  )
] %T>% setorder(class, subclass)
# ----

# ----

# 
# Deletions
# 

print("Processing deletions")
# We need the IDs for each DEL in each patient
deletedMeMin3Vcf = fread(snakemake@input[["repeatMaskerVcfDel"]], skip = "#CHROM")
deletedMeMin3Vcf = deletedMeMin3Vcf[`#CHROM` %in% chrs]
delMeMin3Genotypes = deletedMeMin3Vcf[
  , .SD,
  .SDcols = colnames(deletedMeMin3Vcf)[c(3, 10:(10 + length(samples) - 1))]
]
deletionsIndividiualIds = delMeMin3Genotypes[, .(ID)]
# will create a column of vectors
# position 1 is id for P1 in their file
deletionsIndividiualIds[, ids := delMeMin3Genotypes[, 2:(ncol(delMeMin3Genotypes))] %>%
               apply(2, function (x) {
                 x %>% 
                   stri_split_fixed(pattern = ":") %>% lapply(`[`, y = 8) %>% 
                   stri_split_fixed(pattern = ",") %>% lapply(`[`, y = 1)
               }) %>% simplify2array() %>% asplit(1) %>% lapply(unlist) %>% lapply(unname)]

vcfs = snakemake@input[["vcfs"]]
supports = mclapply(1:length(samples), function (i) {
  # i = 1
  patient = samples[i]
  file = vcfs[[i]]
  vcf = fread(file, skip = "#CHROM")
  setnames(
    vcf, 
    c("#CHROM", "ID", "POS"), 
    c("seqnames", "id", "start"))
  vcf[
    , `:=`(
      svtype = gsub(".*SVTYPE=([^;]+).*", "\\1", vcf$INFO, perl = T)
    )]
  idTable = deletionsIndividiualIds[ # A table with ids of merge and Pi's
    , .(
      ID = ID,
      ids = ids %>% lapply(`[`, y = i) %>% unlist()
    )
  ][order(ID)] # alphanumerical sort to be coherent with other tables
  mergedIdsPresent = idTable[ids != "NaN", ID] # The merged ids of deletions Pi has
  patientIds = idTable$ids %>% .[. != "NaN"] %>% unlist()  # The ids from Pi on their own file
  supportVector = vcf[id %in% patientIds, c(10, 11)] %>% # We take the Genotype columns from both callers in Pi file
    apply(2, function (x) {
      x %>% 
        stri_split_fixed(pattern = ":") %>% 
        lapply(`[`, y = 4) %>% # We select the field with read support
        stri_split_fixed(pattern = ",") %>% 
        lapply(`[`, y = 2) %>% unlist() %>% # and we get the support for the SV
        as.numeric()
    }) %>% 
    apply(1, max) # Then we get the support from the caller that returns the highest
  result = merge(
    idTable, 
    data.table(ID = mergedIdsPresent, support = supportVector),
    by = "ID",
    all.x = T)
  result[, `:=`(p = paste0("P", i), support = fifelse(is.na(support), 0, support))]
  result
}, mc.cores = snakemake@threads) # each element comes out sorted alphanumerically by ID (format: "survX")
supports = supports %>% rbindlist()
supports = supports %>% dcast(ID ~ p, value.var = "support");
setcolorder(supports, colnames(supports) %>% stri_sort(numeric = T))
supportsCol = supports[, -"ID"] %>% simplify2array() %>% asplit(1) %>% lapply(unname) # still sorted

seqids = deletionsBed[order(id)]$id # its sorted alphanumerically
mosdepthDelFiles = snakemake@input[["mosdepthDelFiles"]]
coveragesForGenotypingDel = mclapply(mosdepthDelFiles, function (x) {
  fread(
    x,
    col.names = c("chr", "start", "end", "seqId", "coverage")
  )[order(seqId)]$coverage # also sorted alphanumerically
}, mc.cores = snakemake@threads)
coverageForGenotypingDel = data.table(
  seqId = seqids, # has been sorted alphanumerically
  coverage = lapply(asplit(simplify2array(coveragesForGenotypingDel), 1), as.numeric),
  support = supportsCol
)
coverageForGenotypingDel[, totalReads := mapply(`+`, coverage, support, SIMPLIFY = F)]
genotypesDel = mclapply(1:length(samples), function (i) {
  support = coverageForGenotypingDel$support %>% lapply(`[`, y = i) %>% unlist()
  reads = coverageForGenotypingDel$totalReads %>% lapply(`[`, y = i) %>% unlist()
  
  result = t(genotypeVec(support, reads))
  result
}, mc.cores = snakemake@threads)

genotypesValuesDel = mclapply(1:length(samples), function (i) {genotypesDel[[i]][, 1]}, mc.cores = snakemake@threads)
genotypesDelColumn = genotypesValuesDel %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
meDeletionsMin3 = readRDS(snakemake@input[["meDeletionsMin3"]])
genotypedDel = meDeletionsMin3[data.table(id = seqids, genotype = genotypesDelColumn), on = "id"]
alleleN = (genotypedDel[1, genotype] %>% unlist() %>% length() ) * 2
genotypedDel[, maf := (genotype %>% lapply(sum) %>% unlist() ) / alleleN]

genotypesValuesDel = mclapply(1:length(samples), function (i) {genotypesDel[[i]][, 1]}, mc.cores = snakemake@threads)
genotypesDelColumn = genotypesValuesDel %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
genotypedDel = meDeletionsMin3[data.table(id = seqids, genotype = genotypesDelColumn), on = "id"]
alleleN = (genotypedDel[1, genotype] %>% unlist() %>% length() ) * 2
genotypedDel[, maf := (genotype %>% lapply(sum) %>% unlist() ) / alleleN]

# 
# Write VCF output
# 
print("Writing VCF")
vcf = ioi2[observedMafTable[, c("seqId", "SUPP_geno", "maf")], on = "seqId", nomatch = NULL]
vcf = vcf[insertionsTable[, c("seqId", "vcf_alt")], on = "seqId", nomatch = NULL]
vcf = vcf[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS = start,
  ID = seqId,
  REF = "N",
  ALT = vcf_alt,
  QUAL = ".",
  FILTER = ".",
  INFO = paste(
    "SVTYPE=INS",
    paste0("TE_FAMILY=", repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=", name),
    paste0("SVLEN=", SVLEN),
    paste0("END=", start),
    paste0("SUPP_strict=", SUPP_trusty),
    paste0("SUPP_lax=", SUPP_min3),
    paste0("SUPP_geno=", SUPP_geno),
    paste0("MAF=", maf),
    sep = ";")
)]
fwrite(vcf, snakemake@output[["vcfBody"]], sep = "\t") # Created

ioi3 = annotatedInsertionsMin3[repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") & repeat.percentage >= 0.85][, .SD[1], by = seqId][,
           .(seqId, seqnames, start, end, 
           name, repeat.class, repeat.subclass,
           SVLEN,
           SUPP_min3, SUPP_trusty)]
setkey(ioi3, seqnames, start, end)
vcf2 = ioi3[observedMafTable[, c("seqId", "SUPP_geno", "maf")], on = "seqId", nomatch = NULL]
vcf2 = vcf2[insertionsTable[, c("seqId", "vcf_alt")], on = "seqId", nomatch = NULL]
vcf2 = vcf2[SUPP_geno > 0, .(
  `#CHROM` = seqnames,
  POS = start,
  ID = seqId,
  REF = "N",
  ALT = vcf_alt,
  QUAL = ".",
  FILTER = ".",
  INFO = paste(
    "SVTYPE=INS",
    paste0("TE_FAMILY=", repeat.class),
    paste0("TE_SUBFAMILY=", repeat.subclass),
    paste0("TE_NAME=", name),
    paste0("SVLEN=", SVLEN),
    paste0("END=", start),
    paste0("SUPP_strict=", SUPP_trusty),
    paste0("SUPP_lax=", SUPP_min3),
    paste0("SUPP_geno=", SUPP_geno),
    paste0("MAF=", maf),
    sep = ";")
)]
fwrite(vcf2, snakemake@output[["vcfBodyLax"]], sep = "\t") # Created

saveRDS(insertionsTable, snakemake@output[["insertionsTable"]]) # Modified
saveRDS(annotatedInsertionsMin3,snakemake@output[["annotatedInsertionsMin3"]]) # Modified
saveRDS(meDeletionsMin3, snakemake@output[["meDeletionsMin3"]]) # Modified

saveRDS(genes, snakemake@output[["genes"]]) # Created
saveRDS(observedMafTable, snakemake@output[["observedMafTable"]]) # Created
saveRDS(genotypedDel, snakemake@output[["genotypedDeletions"]]) # Created
