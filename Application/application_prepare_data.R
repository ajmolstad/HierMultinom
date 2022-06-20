library(AnnotatedPBMC)
library(scater)
library(ggplot2)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data"
DATA_PATH = "data/"
dir.create(DATA_PATH, recursive = TRUE)

data = prepare_hao_2020(CACHE_PATH, sce = TRUE)

select_genes = function(sce_list) {

  genes = Reduce(intersect, lapply(sce_list, rownames))

  sce_list = lapply(sce_list, function(x) x[genes, ])

  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))]

  return(genes)

}

data_split = lapply(sort(unique(data$dataset)), function(x) data[, data$dataset == x])
genes = select_genes(data_split)
write.csv(genes, file.path(DATA_PATH, "genes.csv"))