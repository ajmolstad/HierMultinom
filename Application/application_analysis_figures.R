RESULT_PATH = "results/application_Multires_Group_L1"
FIGURES_PATH = file.path(RESULT_PATH, "figures/")
dir.create(FIGURES_PATH, recursive = TRUE)

model = readRDS(file.path(RESULT_PATH, "real_data_application_n_genes_2000_Multires_1.rds"))

Beta = model$Beta_hat
rownames(Beta) = read.csv("data/genes.csv")[1:model$parameters$n_genes, 2]
colnames(Beta) = c("ASDC", "B intermediate",
                   "B memory", "B naive", "CD14 Mono", "CD16 Mono", "CD4 CTL", "CD4 Naive",
                   "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "cDC1",
                   "cDC2", "dnT", "Eryth", "gdT", "HSPC", "ILC", "MAIT", "NK", "NK_CD56bright",
                   "pDC", "Plasmablast", "Platelet", "Treg Memory", "Treg Naive")

groups = list(
  c("B intermediate", "B memory", "B naive", "Plasmablast"),
  c("CD14 Mono", "CD16 Mono"),
  c("CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory", "Treg Naive", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "gdT", "MAIT"),
  c("CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory", "Treg Naive"),
  c("CD4 Naive", "Treg Naive"),
  c("CD4 TCM", "CD4 TEM", "Treg Memory"),
  c("CD8 Naive", "CD8 TCM", "CD8 TEM"),
  c("CD8 TCM", "CD8 TEM"),
  c("NK", "NK_CD56bright"),
  c("ASDC", "cDC1", "cDC2", "pDC"),
  c("cDC1", "cDC2")
)
names(groups) = 1:length(groups)
groups = groups[order(sapply(groups, length), decreasing = TRUE)]

categories = c("B intermediate", "B memory", "B naive", "Plasmablast", "CD14 Mono",
               "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory",
               "Treg Naive", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "gdT",
               "MAIT", "NK", "NK_CD56bright", "ASDC", "cDC1", "cDC2", "pDC",
               "Eryth", "HSPC", "ILC", "Platelet")

Beta = Beta[rowSums(abs(Beta)) != 0, ]

for (i in 1:nrow(Beta)) {

  for (g in names(groups)) {

    if (isTRUE(all.equal(Beta[i, groups[[g]]], setNames(rep(Beta[i, groups[[g]][1]], length(Beta[i, groups[[g]]])), nm = names(Beta[i, groups[[g]]])))) & all(!(Beta[i, groups[[g]]] %in% names(groups)))) {

      Beta[i, groups[[g]]] = as.numeric(g)

    }

  }

}

Beta[!(Beta %in% names(groups))] = 0

unique = Beta[!duplicated(Beta), ]
summary = do.call(rbind, lapply(1:nrow(unique), function(i) Beta[apply(Beta, 1, function(x) all(x == unique[i, ])), ]))
summary[summary == 0] = NA

png(file.path(FIGURES_PATH, "Beta_patterns.png"), width = 6, height = 4, units = "in", res = 1200)
ComplexHeatmap::pheatmap(t(summary[, intersect(categories, unique(unlist(groups)))]), cluster_rows = F, cluster_cols = F, scale = "none", show_colnames = F, color = viridis::viridis(length(groups) + 1), legend = F, border_color = NA, fontsize_row = 8)
dev.off()

png(file.path(FIGURES_PATH, "Beta_patterns1.png"), height = 5, width = 4, units = "in", res = 1200)
ComplexHeatmap::pheatmap((summary[, intersect(categories, unique(unlist(groups)))]), cluster_rows = F, cluster_cols = F, scale = "none", show_rownames = F, color = viridis::viridis(length(groups) + 1), legend = F, border_color = NA, fontsize_col = 8)
dev.off()

marker_genes = read.csv("data/marker_genes.csv")
marker_genes = marker_genes[!grepl("Prolif", marker_genes$Label), ]
rownames(marker_genes) = marker_genes$Label
marker_genes = marker_genes[categories, ]
marker_genes_per_category = lapply(marker_genes$Markers, function(x) strsplit(x, ", ")[[1]][1])
names(marker_genes_per_category) = marker_genes$Label
marker_genes = sapply(unique(unlist(marker_genes_per_category)), function(y) paste(marker_genes$Label[grepl(y, marker_genes_per_category)], collapse = ", "))
marker_genes = marker_genes[names(marker_genes) %in% rownames(Beta)]

Beta_marker_genes = Beta[names(marker_genes), ]
rownames(Beta_marker_genes) = paste0(names(marker_genes), "\n(", marker_genes, ")")
Beta_marker_genes[Beta_marker_genes == 0] = NA

png(file.path(FIGURES_PATH, "Beta_marker_gene_patterns.png"), height = 5, width = 5.5, units = "in", res = 1200)
ComplexHeatmap::pheatmap(t(Beta_marker_genes[, intersect(categories, unique(unlist(groups)))]), cluster_rows = F, cluster_cols = F, scale = "none", color = viridis::viridis(length(groups) + 1), legend = F, fontsize_col = 6, fontsize_row = 8, border_color = "white")
dev.off()

png(file.path(FIGURES_PATH, "Beta_marker_gene_patterns1.png"), width = 5, height = 5.5, units = "in", res = 1200)
ComplexHeatmap::pheatmap((Beta_marker_genes[, intersect(categories, unique(unlist(groups)))]), cluster_rows = F, cluster_cols = F, scale = "none", color = viridis::viridis(length(groups) + 1), legend = F, fontsize_col = 8, fontsize_row = 6, border_color = "white")
dev.off()
