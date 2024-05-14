# -----------------------------------------------
# Create Hao_2020_truncated dataset
# -----------------------------------------------
library(AnnotatedPBMC)
cache_path <- # set path 
data <- AnnotatedPBMC::get_hao_2020(cache_path)
SingleCellExperiment::altExp(data) = NULL
SingleCellExperiment::counts(data) = NULL
data$cell_type = ifelse(data$cell_type_2 == "Treg", data$cell_type_3, data$cell_type_2)
removed_labels = "*Proliferating*"
  data = data[, !grepl(removed_labels, data$cell_type)]
  
select_genes = function(sce_list) {
	genes = Reduce(intersect, lapply(sce_list, rownames))
	sce_list = lapply(sce_list, function(x) x[genes, ])
	ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))
	genes = genes[order(rowMeans(ranks))]
	return(genes)
}

data_split = lapply(sort(unique(data$dataset)), function(x) data[, data$dataset == x])
genes = select_genes(data_split)
# -----------------------------------------------------------------------
# truncated data includes 2000 most variable genes and 
# marker genes from Hao
# -----------------------------------------------------------------------
save_genes = genes[1:2000]
marker_genes = read.csv("marker_genes.csv", header=T)
