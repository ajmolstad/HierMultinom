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
saveName <- paste("~/blue/HierMultinom/Application/results/M0_n", n_train/10000, "_p", n_genes, "_", replicate, ".RDS", sep="")
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

# ------------------------------------------
# MSDA
# ------------------------------------------
Yfac <- as.factor(apply(Y_train, 1, function(x){which(x==1)}))
Yvalfac <- as.factor(apply(Y_validation, 1, function(x){which(x==1)}))
Ytestfac <- as.factor(apply(Y_test, 1, function(x){which(x==1)}))

t0.msda <- msda(x = X_train, y = Yfac, nlambda = 100, verbose = TRUE)
t0.msda.predictval <- predict(t0.msda, X_validation)
val.err <- rep(0, dim(t0.msda.predictval)[2])
for(kk in 1:dim(t0.msda.predictval)[2]){
  val.err[kk] <- sum(Yvalfac == t0.msda.predictval[,kk])/length(Yvalfac)
}
msda.ind <- min(which(val.err == max(val.err)))
msda.classification.error <- sum(Ytestfac != predict(t0.msda, X_test)[,msda.ind])/length(Ytestfac)
msda.theta <- t0.msda$theta[[msda.ind]]

Results <- append(Results, list("msda.theta" = msda.theta, "msda.classification.error" = msda.classification.error))
saveRDS(Results, file = saveName)

# -----------------------------------------
# Random forests with marker genes
# -----------------------------------------
library(randomForest)
Yfac <- as.factor(apply(Y_train, 1, function(x){which(x==1)}))
Yvalfac <- as.factor(apply(Y_validation, 1, function(x){which(x==1)}))
Ytestfac <- as.factor(apply(Y_test, 1, function(x){which(x==1)}))

dat.temp <- data.frame("Y" = Yfac, X_train)
names(dat.temp)[1] <- "Y"
rf <- randomForest(Y~., data=dat.temp, ntree = 1001)
dat.test <- data.frame("Y" = Ytestfac, X_test)
preds <- predict(rf, dat.test)
rf.classification.error <- sum(Ytestfac != preds)/length(Ytestfac)
preds.prob <- predict(rf, dat.test, type="prob")
rf.hellinger <- (1/sqrt(2))*mean(sqrt(rowSums((Y_test-preds.prob)^2)))

dat.temp <- data.frame("Y" = Yfac, XM_train)
names(dat.temp)[1] <- "Y"
rf <- randomForest(Y~., data=dat.temp, ntree = 1001)
dat.test <- data.frame("Y" = Ytestfac, XM_test)
preds <- predict(rf, dat.test)
rf.oracle.classification.error <- sum(Ytestfac != preds)/length(Ytestfac)
preds.prob <- predict(rf, dat.test, type="prob")
rf.oracle.hellinger <- (1/sqrt(2))*mean(sqrt(rowSums((Y_test-preds.prob)^2)))

Results <- append(Results, list(
    "rf.classification.error" = rf.classification.error, 
    "rf.oracle.classification.error" = rf.oracle.classification.error, 
    "rf.hellinger" = rf.hellinger,
		"rf.oracle.hellinger" = rf.oracle.hellinger))
saveRDS(Results, file = saveName)

# ----------------------------------------
# HieRFit
# ----------------------------------------
library(HieRFIT)
TreeTab <- read.delim("~/blue/HierMultinom/Application/data/hao_2020/HierFitTab.txt", header = F)
Yfac <- as.character(apply(Y_train, 1, function(x){which(x==1)}))
Yfac <- paste("X", Yfac, sep="")
Ytestfac <- as.character(apply(Y_test, 1, function(x){which(x==1)}))
Ytestfac <- paste("X", Ytestfac, sep="")

p <- dim(X_train)[2]
colnames(X_train) <- 1:p
rownames(X_train) <- names(Yfac) <- paste("n", 1:dim(X_train)[1], sep="")
ptm <- proc.time()
refmod <- CreateHieR(RefData = Matrix(t(X_train)),
                     ClassLabels = Yfac,
                     Tree = TreeTab,
                     species = "hsapiens")
hierfit.compTime <- proc.time() - ptm

colnames(X_test) <- 1:p
rownames(X_test) <- paste("n", (dim(X_train)[1]+1):(dim(X_train)[1]+ dim(X_test)[1]), sep="")
hierObj <- HieRFIT(Query = Matrix(t(X_test)), refMod = refmod)
evals <- which(hierObj@Evaluation$Projection %in% unique(Ytestfac))
HeirFIT.classification.error <- sum(hierObj@Evaluation$Projection[evals] != Ytestfac[evals])/length(Ytestfac[evals])
HeirFIT.classification.prop <- length(evals)/length(Ytestfac)

Results <- append(Results, list("hierfit.compTime" = hierfit.compTime,
  "HeirFIT.classification.prop" = HeirFIT.classification.prop, 
  "HeirFIT.classification.error" = HeirFIT.classification.error))
saveRDS(Results, file = saveName)

# ----------------------------------------
# Grouped
# ----------------------------------------
library(glmnet)
Ytestfac <- as.factor(apply(Y_test, 1, function(x){which(x==1)}))

ptm <- proc.time()
temp <- glmnet(y = Y_train, x = X_train, family="multinomial", 
               type.multinomial = "grouped", maxit = 1e7, 
  thresh = 1e-9, nlambda = 100, trace.it = T)
val.err <- rep(0, length(temp$lambda))
for(jj in 1:length(temp$lambda)){
  preds <- predict(temp, s=temp$lambda[jj], newx = X_validation, type="response")[,,1]
  val.err[jj] <- -2*sum(log(rowSums(Y_validation*preds)))
  cat(jj,":", val.err[jj], "\n")
}
glmnet.compTime <- proc.time() - ptm

prob.est.glmnet1 <- predict(temp, s = temp$lambda[which(val.err == min(val.err))], newx = X_test, type="response")[,,1]
preds <- apply(prob.est.glmnet1, 1, function(x){which.max(x)})
glmnet.classification.error <- sum(Ytestfac != preds)/length(Ytestfac)
glmnet.deviance <- -2*sum(log(rowSums(Y_test*prob.est.glmnet1)))
glmnet.hellinger <- (1/sqrt(2))*mean(sqrt(rowSums((Y_test - prob.est.glmnet1)^2)))
beta.glmnet <- coef(temp, s = temp$lambda[which(val.err == min(val.err))])
glmnet.validation.error <- val.err[which(val.err == min(val.err))]

Results <- append(Results, 
      list("glmnet.compTime" = glmnet.compTime,
        "glmnet.classification.error" = glmnet.classification.error,
           "glmnet.deviance" = glmnet.deviance,
           "glmnet.hellinger" = glmnet.hellinger, 
           "glmnet.beta" = beta.glmnet,
        "glmnet.validation.error" = glmnet.validation.error))

saveRDS(Results, file = saveName)


q("no")

