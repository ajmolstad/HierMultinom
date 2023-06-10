library(AnnotatedPBMC)
library(ggplot2)
library(ggridges)

# ---------------------------------------------------------
# get hao_2020_truncated from dropbox link in README
# ---------------------------------------------------------
data <- data <- readRDS("/Users/aaron/Dropbox/HierMultinomData/hao_2020_truncated.rds")
data$cell_type = ifelse(data$cell_type_2 == "Treg", data$cell_type_3, data$cell_type_2)
removed_labels = "*Proliferating*"
  data = data[, !grepl(removed_labels, data$cell_type)]

# get MS4A1 logcounts
MS4A1  <- SingleCellExperiment::logcounts(data)[which(SingleCellExperiment::logcounts(data)@Dimnames[[1]] == "MS4A1"),]
Y = as.character(data$cell_type)

# B cell boxplots
dat <- data.frame(MS4A1, Y)
Bs <- which(Y == "B memory" | Y == "B intermediate" | Y == "B naive")
dat <- data.frame("MS4A1" = MS4A1[Bs], "Type" = Y[Bs])

# Basic violin plot
p <- ggplot(dat, aes(x=Type, y= MS4A1, fill=Type)) + ylab("log(Counts)") + xlab("") + 
  geom_violin() + theme_minimal() + theme(legend.position="none") 
  
Bs <- which(Y == "B memory" | Y == "B intermediate" | Y == "B naive")
dat <- data.frame("MS4A1" = MS4A1[Bs], "Type" = Y[Bs])
p1 <- ggplot(dat, aes(x=Type, y= MS4A1, fill=Type)) + ylab("log(Counts)") + xlab("") + ylim(c(0,4.2)) + 
  geom_violin() + theme_minimal() + geom_boxplot(width=0.1, outlier.size = 0.1) + theme(legend.position="none") 
  
Tcells <- c("CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory", "Treg Naive", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "gdT", "MAIT")  
Ts <- which(Y%in%Tcells)  
dat2 <- data.frame("MS4A1" = MS4A1[Ts], "Type" = Y[Ts]) 
p2 <- ggplot(dat2, aes(x=Type, y= MS4A1, fill=Type)) + ylab("log(Counts)") + xlab("") + ylim(c(0,4.2)) +
  geom_violin() + theme_minimal() + theme(legend.position="none") 
  
library(cowplot)  
plot_grid(p1, p2, rel_widths=c(0.3, 0.7))


Monocytes <- c("CD14 Mono", "CD16 Mono")
NK <- c("NK", "NK_CD56bright")
Dendritic <- c("ASDC", "cDC1", "cDC2", "pDC")
Bcells <- c("B memory", "B intermediate", "B naive")
CoarseY <- Y
CoarseY[which(Y%in%Monocytes)] <- "Monocytes"
CoarseY[which(Y%in%NK)] <- "NK"
CoarseY[which(Y%in%Dendritic)] <- "Dendritic"
CoarseY[which(Y%in%Tcells)] <- "T cells"


dat3 <- data.frame("MS4A1" = MS4A1, "Type" = CoarseY) 

p3 <- ggplot(dat3, aes(y=Type, x= MS4A1)) + xlab("log(Counts)") + ylab("Type") + 
  geom_density_ridges() + theme_minimal() + theme(legend.position="none")  
  
library(cowplot)  
plot_grid(p, p3, rel_widths=c(0.3, 0.7))
 




XCL2  <- SingleCellExperiment::logcounts(data)[which(SingleCellExperiment::logcounts(data)@Dimnames[[1]] == "XCL2"),]
Y = as.character(data$cell_type)
CoarseY <- Y
CD4 <- c("CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory", "Treg Naive")
#CoarseY[which(Y%in%NK)] <- "NK"
CoarseY[which(Y%in%Dendritic)] <- "Dendritic"
CoarseY[which(Y%in%Bcells)] <- "B cells"
CoarseY[which(Y%in% CD4)] <- "CD4"  
CoarseY[which(Y%in%Monocytes)] <- "Monocytes"



dat4 <- data.frame("XCL2" = XCL2, "Type" = CoarseY) 
p4 <- ggplot(dat4, aes(y=Type, x= XCL2)) + xlab("log(Counts)") + ylab("Type") +  
  geom_density_ridges() + theme_minimal() + theme(legend.position="none")  
  
tmp <- c("NK", "NK_CD56bright")
dat5 <- data.frame("XCL2" = XCL2[which(Y%in%tmp)], "Type" = Y[which(Y%in%tmp)]) 
p5 <- ggplot(dat5, aes(x=Type, y= XCL2, fill=Type)) + ylab("log(Counts)") + xlab("") + ylim(c(0,4.2)) + 
  geom_violin() + theme_minimal() + theme(legend.position="none") 

plot_grid(p5, p4, rel_widths=c(0.3, 0.7))
