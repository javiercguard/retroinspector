library(data.table)
library(magrittr)
library(ggplot2)
library(patchwork)
library(stringi)
# Enrichment:
library(GO.db) # 3.17 - 3.16
library(AnnotationHub) # 3.8 - 3.6
library(org.Hs.eg.db) # 
library(clusterProfiler) #
library(topGO) #
library(DOSE) #

if (F) {
  # For debugging
  saveRDS(snakemake, "~/path/snk4.rds")
  q()
  snakemake = readRDS("~/path/snk4.rds")
}

genes = readRDS(snakemake@input[["dt"]])

geneNames = genes$annot.symbol %>% unique()
geneDt = bitr(geneNames, fromType = "SYMBOL",
              toType = c("REFSEQ", "ENSEMBL", "ENTREZID"),
              OrgDb = org.Hs.eg.db) %>% data.table()

p = snakemake@params[["plimit"]] %>% as.numeric()

# GO enrichment

for (x in c("MF", "BP", "CC")) {
  name = paste0("ego", x)
  dtName = paste0(name, "Dt")
  assign(x = name,
         enrichGO(
           gene = geneDt$ENSEMBL,
           OrgDb = org.Hs.eg.db,
           keyType = 'ENSEMBL',
           ont = x,
           pAdjustMethod = "BH",
           pvalueCutoff = p)
  )
  object = get(name)
  if (is.null(object)) {
    saveRDS(get(name), snakemake@output[[name]])
    next()
  }
  object = object %>% setReadable(OrgDb = org.Hs.eg.db)
  saveRDS(get(name), snakemake@output[[name]])
}

# KEGG
# disabled because it wont work with the other bioconda packages

# kegg = enrichKEGG(
#   gene = geneDt$ENTREZID,
#   organism = 'hsa',
#   pvalueCutoff = p
# )
# saveRDS(kegg, snakemake@output[["kegg"]])

# DO

do <- enrichDO(gene          = geneDt$ENTREZID,
               ont           = "DO",
               pvalueCutoff  = p,
               pAdjustMethod = "BH",
               # universe      = names(geneList),
               minGSSize     = 5,
               maxGSSize     = 500,
               readable      = T)
saveRDS(do, snakemake@output[["do"]])

# NCG

ncg = enrichNCG(geneDt$ENTREZID %>% unique(),
                pvalueCutoff = p,
                readable = T)
saveRDS(ncg, snakemake@output[["ncg"]])
