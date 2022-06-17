# --------------------------------------------
# Get data
# --------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(uu)
library(glmnet)
nreps <- 100
t0 <- expand.grid(p = rep(c(100, 200, 500, 1000), each=nreps),
  Model = c(1, 2, 3, 4, 5, 6))
p <- t0[uu,1]
Model <- t0[uu,2]
savename <- paste("~/blue/HierMultinom/Simulations/ResultsOverlap/Model", Model, "_p", p, "_", uu %% nreps + 1, ".RDS", sep="")
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

source("~/blue/HierMultinom/Simulations/GenerateDataOverlap.R")
source("~/blue/HierMultinom/Simulations/HierMultinomNestedOverlap.R")

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

# -----------------------
# fine-pred models
# ------------------------
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









# --------------------------------------------------------------------
# For fixed lambda, need to workout lambda yielding complete sparsity
# --------------------------------------------------------------------
groups.input <- list(1:6, 7:9, 10:12, 4:6, 1:3)
t0 <- HierMultNestedOverlap.Path(X, Y, groups.input, ngamma = 50, delta = 0.01, 
                                 lambda.vec = 10^seq(-5,-1, length=10), tol = 1e-9, 
                                 max.iter = 1000, Xval, Yval)
beta.est <- t0$beta
x.sd.temp <- apply(X, 2, sd)
XtestInput <- cbind(1, (Xtest - rep(1, dim(Xtest)[1])%*%t(apply(X, 2, mean)))/(rep(1, dim(Xtest)[1])%*%t(x.sd.temp)))
if(any(x.sd.temp==0)){
  XtestInput[,which(x.sd.temp==0)] <- 0
}
l0 <- XtestInput%*%beta.est
prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
Results  <- append(Results, getResults(Ytest, Xtest, beta, prob.est, "Ours"))


saveRDS(Results, file=savename)