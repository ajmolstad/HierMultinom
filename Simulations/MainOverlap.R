# --------------------------------------------
# Get data
# --------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(uu)
library(glmnet)
library(Matrix)
library(HierMultinom) 
nreps <- 100
t0 <- expand.grid(p = rep(c(100, 200, 500, 1000), each=nreps),
  Model = c(1, 2, 3, 4, 5, 6))
p <- t0[uu,1]
Model <- t0[uu,2]
savename <- paste("~/blue/HierMultinom/Simulations_R1/ResultsOverlap/Complete_Model", Model, "_p", p, "_", uu %% nreps + 1, ".RDS", sep="")
n <- 500
K <- 12
SigmaX <- matrix(0, nrow=p, ncol=p)
for(jj in 1:p){
  for(kk in 1:p){
    SigmaX[jj,kk] <- 0.7^abs(jj-kk)
  }
}
eo <- eigen(SigmaX) 
SigmaXsqrt <- eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
X <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
Xval <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
Xtest <- matrix(rnorm(1e4*p), nrow=1e4)%*%SigmaXsqrt
groups <- list(
  7:9, 10:12, 1:6, 4:6, 1:3
)
source("~/blue/HierMultinom/Simulations_R1/GenerateDataOverlap.R")

getResults <- function(Ytest, Xtest, beta, prob.est, method){
  l0 <- exp(Xtest%*%beta)
  probs.true <- l0/rowSums(l0)
  hellinger <- sum(sqrt(rowSums((sqrt(prob.est) - sqrt(probs.true))^2)))/(sqrt(2)*dim(Xtest)[1])
  classification.error <- sum(apply(Ytest, 1, which.max) != apply(prob.est, 1, which.max))/dim(Ytest)[1]
  KL <- sum(log(prob.est/probs.true)*prob.est)/dim(Ytest)[1]
  Results <- list("hellinger" = hellinger, "classErr" = classification.error, "KL" = KL)
  names(Results) <- paste(method, names(Results), sep=".")
  return(Results)
}
eo <- NULL
SigmaX <- NULL
SigmaXsqrt <- NULL

# ------------------------------------------------------
# glmnet (grouped)
# ------------------------------------------------------
temp <- glmnet(y = Y, x = X, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(temp$lambda))
for(jj in 1:length(temp$lambda)){
  preds <- predict(temp, s=temp$lambda[jj], newx = Xval, type="response")[,,1]
  val.err[jj] <- -2*sum(log(rowSums(Yval*preds)))
}
prob.est.glmnet1 <- predict(temp, s = temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]
Results  <- getResults(Ytest, Xtest, beta, prob.est.glmnet1, "glmnetGroup")
beta.glmnet <- coef(temp, s = temp$lambda[which(val.err == min(val.err))])

# ------------------------------------------------------
# glmnet (ungrouped)
# ------------------------------------------------------
temp <- glmnet(y = Y, x = X, family="multinomial", type.multinomial = "ungrouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(temp$lambda))
for(jj in 1:length(temp$lambda)){
  preds <- predict(temp, s=temp$lambda[jj], newx = Xval, type="response")[,,1]
  val.err[jj] <- -2*sum(log(rowSums(Yval*preds)))
}
prob.est.glmnet1 <- predict(temp, s = temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est.glmnet1, "glmnetL1"))

# -------------------------------------------------------
# Two-step approach
# -------------------------------------------------------
groups.agg <- list(
  1:3, 4:6, 7:9, 10:12
)
coarseY <- matrix(0, nrow=n, ncol=length(groups.agg))
coarseYval <- matrix(0, nrow=n, ncol=length(groups.agg))
for(jj in 1:n){
  for(kk in 1:length(groups.agg)){
    if(any(Y[jj,groups.agg[[kk]]] == 1)){
      coarseY[jj,kk] <- 1
    }
    if(any(Yval[jj,groups.agg[[kk]]] == 1)){
      coarseYval[jj,kk] <- 1
    }
  }
}

coarse.temp <- glmnet(y = coarseY, x = X, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(coarse.temp$lambda))
for(jj in 1:length(coarse.temp$lambda)){
  preds <- predict(coarse.temp, s=coarse.temp$lambda[jj], newx = Xval, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(coarseYval*preds)))
}
coarse.preds <- predict(coarse.temp, s=coarse.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]
beta.coarse <- coef(coarse.temp, s=coarse.temp$lambda[which(val.err == min(val.err))])


# -----------------------
# fine-pred models
# ------------------------
Ytemp <- Y[which(rowSums(Y[,1:3])==1),1:3]
Xtemp <- X[which(rowSums(Y[,1:3])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,1:3])==1),1:3]
Xtemp.val <- Xval[which(rowSums(Yval[,1:3])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds1 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, thresh = 1e-9, type="response")[,,1]
beta.fine1 <- coef(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))])

Ytemp <- Y[which(rowSums(Y[,4:6])==1),4:6]
Xtemp <- X[which(rowSums(Y[,4:6])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,4:6])==1),4:6]
Xtemp.val <- Xval[which(rowSums(Yval[,4:6])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds2 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, thresh = 1e-9, type="response")[,,1]
beta.fine2 <- coef(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))])

Ytemp <- Y[which(rowSums(Y[,7:9])==1),7:9]
Xtemp <- X[which(rowSums(Y[,7:9])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,7:9])==1),7:9]
Xtemp.val <- Xval[which(rowSums(Yval[,7:9])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds3 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]
beta.fine3 <- coef(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))])


Ytemp <- Y[which(rowSums(Y[,10:12])==1),10:12]
Xtemp <- X[which(rowSums(Y[,10:12])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,10:12])==1),10:12]
Xtemp.val <- Xval[which(rowSums(Yval[,10:12])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds4 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]
beta.fine4 <- coef(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))])


preds.two.step <- cbind(fine.preds1*coarse.preds[,1], fine.preds2*coarse.preds[,2], fine.preds3*coarse.preds[,3], fine.preds4*coarse.preds[,4]) 
Results  <- append(Results, getResults(Ytest, Xtest, beta, preds.two.step, "twoStepGrouped"))







# -------------------------------------------------------
# Two-step approach
# -------------------------------------------------------
groups.agg <- list(
  1:6, 7:9, 10:12
)
coarseY <- matrix(0, nrow=n, ncol=length(groups.agg))
coarseYval <- matrix(0, nrow=n, ncol=length(groups.agg))
for(jj in 1:n){
  for(kk in 1:length(groups.agg)){
    if(any(Y[jj,groups.agg[[kk]]] == 1)){
      coarseY[jj,kk] <- 1
    }
    if(any(Yval[jj,groups.agg[[kk]]] == 1)){
      coarseYval[jj,kk] <- 1
    }
  }
}

coarse.temp <- glmnet(y = coarseY, x = X, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(coarse.temp$lambda))
for(jj in 1:length(coarse.temp$lambda)){
  preds <- predict(coarse.temp, s=coarse.temp$lambda[jj], newx = Xval, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(coarseYval*preds)))
}
coarse.preds <- predict(coarse.temp, s=coarse.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]



# ------------------------------------------
# Fine-prediction models
# ------------------------------------------
Ytemp <- Y[which(rowSums(Y[,1:6])==1),1:6]
Xtemp <- X[which(rowSums(Y[,1:6])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,1:6])==1),1:6]
Xtemp.val <- Xval[which(rowSums(Yval[,1:6])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds1 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]

Ytemp <- Y[which(rowSums(Y[,7:9])==1),7:9]
Xtemp <- X[which(rowSums(Y[,7:9])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,7:9])==1),7:9]
Xtemp.val <- Xval[which(rowSums(Yval[,7:9])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds2 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]


Ytemp <- Y[which(rowSums(Y[,10:12])==1),10:12]
Xtemp <- X[which(rowSums(Y[,10:12])==1),]
Ytemp.val <- Yval[which(rowSums(Yval[,10:12])==1),10:12]
Xtemp.val <- Xval[which(rowSums(Yval[,10:12])==1),]
fine.temp <- glmnet(y = Ytemp, x = Xtemp, family="multinomial", type.multinomial = "grouped", thresh = 1e-9, maxit = 1e6)
val.err <- rep(0, length(fine.temp$lambda))
for(jj in 1:length(fine.temp$lambda)){
  preds <- predict(fine.temp, s=fine.temp$lambda[jj], newx = Xtemp.val, type="response")[,,1]
  val.err[jj] <-  -2*sum(log(rowSums(Ytemp.val*preds)))
}
fine.preds3 <- predict(fine.temp, s=fine.temp$lambda[which(val.err == min(val.err))], newx = Xtest, type="response")[,,1]


preds.two.step <- cbind(fine.preds1*coarse.preds[,1], fine.preds2*coarse.preds[,2], fine.preds3*coarse.preds[,3]) 
Results  <- append(Results, getResults(Ytest, Xtest, beta, preds.two.step, "twoStepGroupedcoarse"))


# ------------------------------------------
# MSDA
# ------------------------------------------
library(msda)
Yfac <- as.factor(apply(Y, 1, function(x){which(x==1)}))
Yvalfac <- as.factor(apply(Yval, 1, function(x){which(x==1)}))
Ytestfac <- as.factor(apply(Ytest, 1, function(x){which(x==1)}))
if(any(table(Yfac) <= 2)){
  rm.cat <- which(colSums(Y) <= 2)
  rm.ind <- which(Yfac == rm.cat)
  Yfac <- apply(Y[-rm.ind,], 1, function(x){which(x==1)})
  adj <- which(apply(Y[-rm.ind,], 1, function(x){which(x==1)}) > rm.cat)
  Yfac[adj] <- Yfac[adj]- 1
  X.tmp <- X[-rm.ind,]
  makelast <- which(Yvalfac == rm.cat)
  Yvalfac <- apply(Yval, 1, function(x){which(x==1)})
  adj <- which(apply(Yval, 1, function(x){which(x==1)}) > rm.cat)
  Yvalfac[adj] <-  Yvalfac[adj] - 1
  Yvalfac[makelast] <- 12
  makelast <- which(Ytestfac == rm.cat)
  Ytestfac <- apply(Ytest, 1, function(x){which(x==1)})
  adj <- which(apply(Ytest, 1, function(x){which(x==1)}) > rm.cat)
  Ytestfac[adj] <-  Ytestfac[adj] - 1
  Ytestfac[makelast] <- 12
  t0.msda <- msda(x = X.tmp, y = Yfac, nlambda = 100)
} else {
  t0.msda <- msda(x = X, y = Yfac, nlambda = 100)
}
t0.msda.predictval <- predict(t0.msda, Xval)
val.err <- rep(0, dim(t0.msda.predictval)[2])
for(kk in 1:dim(t0.msda.predictval)[2]){
  val.err[kk] <- sum(Yvalfac == t0.msda.predictval[,kk])
}
msda.ind <- min(which(val.err == max(val.err)))
msda.classification.error <- sum(Ytestfac != predict(t0.msda, Xtest)[,msda.ind])/length(Ytestfac)
Results <- append(Results, list("msda.classErr" = msda.classification.error))
theta.msda <- t0.msda$theta[[msda.ind]]


# --------------------------------
# Random forests 
# --------------------------------
library(randomForest)
imp.preds <- which(rowSums(abs(beta))!=0)
Yfac <- as.factor(apply(Y, 1, function(x){which(x==1)}))
Yvalfac <- as.factor(apply(Yval, 1, function(x){which(x==1)}))
Ytestfac <- as.factor(apply(Ytest, 1, function(x){which(x==1)}))

dat.temp <- data.frame(as.factor(apply(Y, 1, function(x){which(x==1)})), X)
names(dat.temp)[1] <- "Y"
rf <- randomForest(Y~., data=dat.temp)
dat.test <- data.frame(as.factor(apply(Y, 1, function(x){which(x==1)})), Xtest)
preds <- predict(rf, dat.test)
rf.classification.error <- sum(Ytestfac != preds)/length(Ytestfac)
rf <- NULL

varX <- sort(apply(X, 2, sd), decreasing=T, index=T)$ix[1:50]
dat.temp <- data.frame(as.factor(apply(Y, 1, function(x){which(x==1)})), X[,unique(c(varX, imp.preds))])
names(dat.temp)[1] <- "Y"
rf.oracle <- randomForest(Y~., data=dat.temp)
dat.test <- data.frame(as.factor(apply(Y, 1, function(x){which(x==1)})), Xtest[,unique(c(varX, imp.preds))])
preds <- predict(rf.oracle, dat.test)
rf.oracle.classification.error <- sum(Ytestfac != preds)/length(Ytestfac)
Results <- append(Results, list("rf.classErr" = rf.classification.error, 
                                "rf.oracle.classErr" = rf.oracle.classification.error))
rf.oracle <- NULL


# --------------------------------
# Hierarchical random forests 
# --------------------------------
TreeTab <- read.delim("~/blue/HierMultinom/Simulations_R1/TreeTabOverlap.txt", header=F)
library(HieRFIT, lib.loc="~/blue/HierMultinom/Simulations_R1/packages/")

Yfac <- as.character(apply(Y, 1, function(x){which(x==1)}))
Yfac <- paste("X", Yfac, sep="")
Ytestfac <- as.character(apply(Ytest, 1, function(x){which(x==1)}))
Ytestfac <- paste("X", Ytestfac, sep="")

colnames(X) <- 1:p
rownames(X) <- names(Yfac) <- paste("n", 1:dim(X)[1], sep="")
ptm <- proc.time()
refmod <- CreateHieR(RefData = Matrix(t(X)),
                     ClassLabels = Yfac,
                     Tree = TreeTab,
                     species = "hsapiens")
comp.time.hrf <- proc.time() - ptm

colnames(Xtest) <- 1:p
rownames(Xtest) <- paste("n", (dim(X)[1]+1):(dim(X)[1]+ dim(Xtest)[1]), sep="")
hierObj <- HieRFIT(Query = Matrix(t(Xtest)), refMod = refmod)
evals <- which(hierObj@Evaluation$Projection %in% unique(Ytestfac))
HeirFIT.classification.error <- sum(hierObj@Evaluation$Projection[evals] != Ytestfac[evals])/length(Ytestfac[evals])
HeirFIT.classification.prop <- length(evals)/length(Ytestfac)

Results <- append(Results, list("HeirFIT.classErr" = HeirFIT.classification.error, 
                                "HeirFIT.classProp" = HeirFIT.classification.prop))
refmod <- NULL

# --------------------------------------------------------
# Our method with only multiresolution penalty
# ---------------------------------------------------------
sourceCpp("~/blue/HierMultinom/Simulations_R1/Functions/matMult.cpp")
source("~/blue/HierMultinom/Simulations_R1/Functions/HierMultinomOverlapZero.R")
groups.input <- list(1:6, 7:9, 10:12, 4:6, 1:3)
t0 <- HierMultinomOverlapZero.path(X[,imp.preds], Y, groups.input, lambda.vec = 10^seq(1,-3, length=25), tol = 1e-8, max.iter = 5000, Xval[,imp.preds], Yval)
beta.est.oracle <- t0$beta
XtestInput <- cbind(1, (Xtest[,imp.preds] - rep(1, dim(Xtest)[1])%*%t(apply(X[,imp.preds], 2, mean)))/(rep(1, dim(Xtest)[1])%*%t(apply(X[,imp.preds], 2, sd))))
l0 <- exp(XtestInput%*%beta.est.oracle)
prob.est <- l0/rowSums(l0)
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est, "OursOracle"))



# --------------------------------------------------------------------
# Our method
# --------------------------------------------------------------------
groups.input <- list(1:6, 7:9, 10:12, 4:6, 1:3)
ptm <- proc.time()
t0 <- HierMultinomOverlap.path(X, Y, groups.input, ngamma = 50, delta = 0.01, 
                                 lambda.vec = 10^seq(-5,-1, length=10), tol = 1e-9, 
                                 max.iter = 1000, Xval, Yval)
comp.time <- proc.time() - ptm

beta.est <- t0$beta.est
x.sd.temp <- apply(X, 2, sd)
XtestInput <- cbind(1, (Xtest - rep(1, dim(Xtest)[1])%*%t(apply(X, 2, mean)))/(rep(1, dim(Xtest)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est, "Ours"))

# -----------------------------------------------
# Our method with relaxed convergence tolerance
# -----------------------------------------------
Results <- append(Results, list("beta" = beta, "beta.est" = beta.est, "beta.glmnet" = beta.glmnet, "theta.msda" = theta.msda,
  "comp.time" = comp.time, "comp.time.hrf" = comp.time.hrf))


ptm <- proc.time()
t0 <- HierMultinomOverlap.path(X, Y, groups.input, ngamma = 50, delta = 0.01, 
                                 lambda.vec = 10^seq(-5,-1, length=10), tol = 1e-7, 
                                 max.iter = 1000, Xval, Yval)
comp.time.e7 <- proc.time() - ptm

beta.est <- t0$beta.est
x.sd.temp <- apply(X, 2, sd)
XtestInput <- cbind(1, (Xtest - rep(1, dim(Xtest)[1])%*%t(apply(X, 2, mean)))/(rep(1, dim(Xtest)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est, "Ours.e7"))
Results <- append(Results, list("comp.time.e7" = comp.time.e7, "beta.est.e7" = Matrix(beta.est, sparse=TRUE)))


# -----------------------------------------------
# Our method with more relaxed convergence tolerance
# -----------------------------------------------
Results <- append(Results, list("beta" = beta, "beta.est" = beta.est, "beta.glmnet" = beta.glmnet, "theta.msda" = theta.msda,
  "comp.time" = comp.time, "comp.time.hrf" = comp.time.hrf))


ptm <- proc.time()
t0 <- HierMultinomOverlap.path(X, Y, groups.input, ngamma = 50, delta = 0.01, 
                                 lambda.vec = 10^seq(-5,-1, length=10), tol = 1e-5, 
                                 max.iter = 1000, Xval, Yval)
comp.time.e5 <- proc.time() - ptm

beta.est <- t0$beta.est
x.sd.temp <- apply(X, 2, sd)
XtestInput <- cbind(1, (Xtest - rep(1, dim(Xtest)[1])%*%t(apply(X, 2, mean)))/(rep(1, dim(Xtest)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est, "Ours.e5"))
Results <- append(Results, list("comp.time.e5" = comp.time.e5, "beta.est.e5" = Matrix(beta.est, sparse=TRUE)))



saveRDS(Results, file=savename)