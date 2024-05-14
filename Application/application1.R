# --------------------------------------
# Setup analysis
# --------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(dplyr)
library(AnnotatedPBMC)
library(Matrix)
library(msda)
library(randomForest)

model.parameters1 <- expand.grid(n_train = 20000, n_genes = rep(c(500, 750, 1500, 2000), each=50))
model.parameters1$replicate <- rep(1:50, length=200)
model.parameters2 <- expand.grid(n_train = rep(1:5*10000, each=50), n_genes = 1000)
model.parameters2$replicate <- rep(1:50, length=250)
parameters <- bind_rows(model.parameters1, model.parameters2)


n_genes <- parameters[uu,2]
n_train <- parameters[uu,1]
replicate <- parameters[uu,3]
saveName <- paste("~/blue/HierMultinom/Application/results/M1_n", n_train/10000, "_p", n_genes, "_", replicate, ".RDS", sep="")
# data = readRDS("~/blue/HierMultinom/Application/data/hao_2020/hao_2020.rds")
# SingleCellExperiment::altExp(data) = NULL
# SingleCellExperiment::counts(data) = NULL
marker.genes <- read.csv("~/blue/HierMultinom/Application/data/hao_2020/marker_genes.csv")$Markers
marker.genes <- unique(unlist(lapply(marker.genes, function(x){unlist(strsplit(x, ", "))})))
genes = read.csv("~/blue/HierMultinom/Application/data/hao_2020/genes.csv")[, 2][1:n_genes]
# truncatedData =  data[unique(c(marker.genes, genes[1:5000])),]
# saveRDS(truncatedData, file="/Users/aaron/Desktop/HierMultinom/hao_2020/hao_2020_truncated.rds")

data = readRDS("~/blue/HierMultinom/Application/data/hao_2020/hao_2020_truncated.rds")
dataM = data[marker.genes, ]
data = data[genes,]
data$cell_type = ifelse(data$cell_type_2 == "Treg", data$cell_type_3, data$cell_type_2)
dataM$cell_type = ifelse(dataM$cell_type_2 == "Treg", dataM$cell_type_3, dataM$cell_type_2)
removed_labels = "*Proliferating*"
dataM = dataM[, !grepl(removed_labels, dataM$cell_type)]
data = data[, !grepl(removed_labels, data$cell_type)]

set.seed(replicate)
n_validation <- 20000
n_test <- 20000
sample.points <- sample(1:ncol(data), n_train + n_validation + n_test, replace = FALSE)
data = data[, sample.points]
dataM = dataM[, sample.points]

X = t(as.matrix(SingleCellExperiment::logcounts(data)))
XM = t(as.matrix(SingleCellExperiment::logcounts(dataM)))
Y = as.character(data$cell_type)
categories = sort(unique(Y))
Y = t(sapply(Y, function(y) categories == y) * 1)
colnames(Y) = categories

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
groups = lapply(groups, function(group) sapply(group, function(x) which(categories == x)))

train_indices = 1:n_train
validation_indices = (n_train + 1):(n_train + n_validation)
test_indices = (n_train + n_validation + 1):(n_train + n_validation + n_test)

X_train = X[train_indices, ]
X_validation = X[validation_indices, ]
X_test = X[test_indices, ]
XM_train = XM[train_indices, ]
XM_validation = XM[validation_indices, ]
XM_test = XM[test_indices, ]
Y_train = Y[train_indices, ]
Y_validation = Y[validation_indices, ]
Y_test = Y[test_indices, ]
groups = groups


if(any(apply(XM_train,2,sd) == 0)){
	rm.genes <- which(apply(XM_train, 2, sd) == 0)
	XM_train <- XM_train[,-rm.genes]
	XM_validation <- XM_validation[,-rm.genes]
	XM_test <- XM_test[,-rm.genes]
}

if(any(apply(X_train,2,sd) == 0)){
	rm.genes <- which(apply(X_train, 2, sd) == 0)
	X_train <- X_train[,-rm.genes]
	X_validation <- X_validation[,-rm.genes]
	X_test <- X_test[,-rm.genes]
}

Results <- NULL

# -------------------------------------
# oracle Ours 
# -------------------------------------
sourceCpp("~/blue/HierMultinom/Simulations_R1/Functions/matMult.cpp")
source("~/blue/HierMultinom/Simulations_R1/Functions/HierMultinomOverlapZero.R")
ptm <- proc.time()
t0 <- HierMultinomOverlapZero.path(XM_train, Y_train, groups, 
	lambda.vec = c(0, 10^seq(-5,-1, length = 25)), tol = 1e-7, 
	max.iter = 1e4, XM_validation, Y_validation)
oursOracle.compTime <- proc.time() - ptm

beta.est <- t0$beta
x.sd.temp <- apply(XM_train, 2, sd)
XtestInput <- cbind(1, (XM_test - rep(1, dim(XM_test)[1])%*%t(apply(XM_train, 2, mean)))/(rep(1, dim(XM_test)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
preds <- apply(prob.est, 1, which.max)
oursOracle.classification.error <- sum(apply(Y_test, 1, which.max) != preds)/dim(Y_test)[1]
oursOracle.deviance <- -2*sum(log(rowSums(Y_test*prob.est)))
oursOracle.hellinger <- (1/sqrt(2))*mean(sqrt(rowSums((Y_test - prob.est)^2)))
oursOracle.beta <- beta.est

Results <- append(Results, 
                  list("oursOracle.compTime" = oursOracle.compTime,
                    "oursOracle.classification.error" = oursOracle.classification.error,
                       "oursOracle.deviance" = oursOracle.deviance,
                       "oursOracle.hellinger" = oursOracle.hellinger, 
                       "oursOracle.beta" = oursOracle.beta))
saveRDS(Results, file = saveName)

# -------------------------------------
# Ours 
# -------------------------------------
library(HierMultinom)
ptm <- proc.time()
t0 <- HierMultinomOverlap.path(X_train, Y_train, groups, 
	ngamma = 100, delta = 1e-04, tol = 1e-7,
	lambda.vec = c(10^seq(-5,-3, length=5))[3:4], 
	max.iter = 1e4,
	X_validation, Y_validation)
ours.compTime <- proc.time() - ptm

ours.validation.error <- min(t0$val.err)
beta.est <- t0$beta
x.sd.temp <- apply(X_train, 2, sd)
XtestInput <- cbind(1, (X_test - rep(1, dim(X_test)[1])%*%t(apply(X_train, 2, mean)))/(rep(1, dim(X_test)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
preds <- apply(prob.est, 1, which.max)
ours.classification.error <- sum(apply(Y_test, 1, which.max) != preds)/dim(Y_test)[1]
ours.deviance <- -2*sum(log(rowSums(Y_test*prob.est)))
ours.hellinger <- (1/sqrt(2))*mean(sqrt(rowSums((Y_test - prob.est)^2)))
ours.beta <- beta.est


Results <- append(Results, 
                  list("ours.compTime" = ours.compTime,
                       "our.classification.error" = ours.classification.error,
                       "ours.deviance" = ours.deviance,
                       "ours.hellinger" = ours.hellinger, 
                       "ours.beta" = ours.beta,
                       "ours.validation.error" = ours.validation.error))
saveRDS(Results, file = saveName)

q("no")

