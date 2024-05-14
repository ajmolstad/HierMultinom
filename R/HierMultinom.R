library(MASS)
library(Rcpp)

# -------------------------------------------------------
# Compute gradient and log-likelihood
# -------------------------------------------------------
getGradLik <- function(X, Y, beta){
  t0 <- eigenMatMult(X, beta)
  if(any(t0 > 700)){ #710 is exact threshold for Inf
    t1 <- exp(t0 - apply(t0, 1, max))
    P <- t1/eigenMatMult(armaRowSums(t1), t(rep(1, dim(beta)[2])))
  } else {
    P <- exp(t0)/eigenMatMult(armaRowSums(exp(t0)), t(rep(1, dim(beta)[2])))
  }
  grad <- crossprod(X, P - Y)/dim(Y)[1]
  return(list("grad" = grad, "lik" = -sum(log(armaRowSums(Y*P)))/dim(Y)[1]))
}

# -------------------------------------------------------
# Compute objective function value
# -------------------------------------------------------
getObj <- function(X, Y, beta, groups, lambda, gamma){
  t0 <- eigenMatMult(X, beta)
  if(any(t0 > 700)){ #710 is exact threshold for Inf
    t1 <- exp(t0 - apply(t0, 1, max))
    P <- t1/eigenMatMult(armaRowSums(t1), t(rep(1, dim(beta)[2])))
    t1 <- NULL
  } else {
    P <- exp(t0)/eigenMatMult(armaRowSums(exp(t0)), t(rep(1, dim(beta)[2])))
  }
  t0 <- -sum(log(armaRowSums(Y*P)))/dim(Y)[1]
  if(lambda == 0 & gamma == 0){
    return(t0)
  } else {
    p1 <- gamma*sum(apply(beta[-1,], 1, function(x){sqrt(sum(x^2))}))
    p2 <- 0
    for(kk in 1:length(groups)){
      p2 <- p2 + sum(apply(beta[-1,groups[[kk]]] - rowMeans(beta[-1,groups[[kk]], drop=FALSE])%*%t(rep(1, length(groups[[kk]]))), 1, function(x){sqrt(sum(x^2))}))
    }
    return(t0 + p1 + lambda*p2)
  }
}

# -------------------------------------------------------
# Evaluate penalty
# -------------------------------------------------------
evalPen <- function(beta, lambda, gamma, groups){
  p1 <- gamma*sum(apply(beta[-1,], 1, function(x){sqrt(sum(x^2))}))
  p2 <- 0
  for(kk in 1:length(groups)){
    p2 <- p2 + sum(apply(beta[-1,groups[[kk]]] - rowMeans(beta[-1,groups[[kk]], drop=FALSE])%*%t(rep(1, length(groups[[kk]]))), 1, function(x){sqrt(sum(x^2))}))
  }
  return(p1 + lambda*p2)
}

# ------------------------------------------
# Compute first part of prox operator
# ------------------------------------------
proxOp1 <- function(y, lambda) {
  s <- length(y)
  tildeW <- diag(1, s) - matrix(1, s, s)/s
  v <- crossprod(svd(tildeW)$u, y)[-s]
  t1 <- crossprod(tildeW, y)
  t0 <- sum(t1^2)
  if (t0 <= lambda^2) {
    result <- mean(y)*rep(1, s)
  } else {
    result <- y - crossprod(tildeW, y)/(sqrt(sum(v^2))/lambda)
  }
  return(result)
}


# ------------------------------------------
# Compute second part of prox operator
# ------------------------------------------
outerProx <- function(y, lambda, gamma, groups){
  b.out <- y
  if (lambda!=0) {
    for (jj in 1:length(groups)) {
      b.out[groups[[jj]]] <- proxOp1(y[groups[[jj]]], lambda)
    }
  } else {
    b.out <- y
  }
  t0 <- sqrt(sum(b.out^2))
  if (t0 > gamma) {
    b.out <- (1 - gamma/t0)*b.out
  } else {
    b.out <- rep(0, length(y))
  }
  return(b.out)
  
}

# --------------------------------------------
# Compute complete proximal operator
# --------------------------------------------
fullProx <- function(input, lambda, gamma, groups){
  out <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])
  out[1,] <- input[1,]
  for (j in 2:dim(input)[1]) {
    if (sqrt(sum(input[j,]^2)) < gamma) {
      out[j,] <- 0
    } else {
      out[j,] <- outerProx(input[j,], lambda, gamma, groups)
    }
  }
  return(out)
}

# ---------------------------------------------------------
# Evaluate the negative log-likelihood on validation set
# --------------------------------------------------------
validationLik <- function(X, Xval, Yval, beta){
  Xval.temp <- cbind(1, (Xval - tcrossprod(rep(1, dim(Xval)[1]), colMeans(X)))/tcrossprod(rep(1, dim(Xval)[1]), apply(X, 2, sd)))
  lleval <- crossprod(t(Xval.temp),beta)
  probs <- exp(lleval)/tcrossprod(rowSums(exp(lleval)), rep(1, dim(beta)[2]))
  return(-2*sum(log(rowSums(Yval*probs))))
}


# ---------------------------------------------------------------
# Minimize objective function for one pair of tuning parameters
# --------------------------------------------------------------
AccPGD.HierMult <- function(X, Y, groups, beta.init, 
      lambda, gamma, max.iter = 1e4, tol = 1e-9){
  
  obj.vals <- rep(Inf, max.iter)
  obj.vals[1] <- getObj(X, Y, beta.init, groups, lambda, gamma)
  beta.iter <- beta.init
  beta.prev <- beta.iter
  alpha <- 1
  alpha.prev <- 1
  L0 <- (dim(Y)[1])/(sqrt(dim(Y)[2])*sum(X^2))
  step.size <- 1e4*L0
  lik.current <- getObj(X, Y, beta.iter, groups, lambda = 0, gamma = 0)
  for (kk in 1:max.iter) {
    
    search.point <- beta.iter + ((alpha.prev - 1)/(alpha))*(beta.iter - beta.prev)
    temp <- getGradLik(X, Y, search.point)
    temp.grad <- temp$grad
    lik.current <- temp$lik
    linesearch <- TRUE
    
    while (linesearch) {
      
      temp.step <- search.point - step.size*temp.grad
      beta.try <- fullProx(temp.step, step.size*lambda, step.size*gamma, groups)
      lik.try <- getObj(X, Y, beta.try, groups, lambda = 0, gamma = 0)
      tmp <- sum(temp.grad*(beta.try - search.point))
      check <- (lik.try <= lik.current + tmp + (1/(2*step.size))*sum((beta.try - search.point)^2))
      if (check | step.size == L0) {
        
        linesearch <- FALSE
        beta.prev <- beta.iter
        beta.iter <- beta.try
        alpha.prev <- alpha
        alpha <- (1 + sqrt(1 + 4*alpha.prev^2))/2
        
      } else {
        
        step.size <- max(step.size*.5, L0)
        
      }
    }
    
    obj.vals[kk+1] <- lik.try + evalPen(beta.iter, lambda, gamma, groups)
    if (kk > 10) {
      if (all(abs(obj.vals[kk:(kk-3)] - obj.vals[(kk+1):(kk-2)]) < tol * abs(obj.vals[1]))) {
        break
      }
      if (kk > max.iter) {
        break
      }
    }
  }
  return(beta.iter)
}


# ---------------------------------------------------------------
# Compute solution path, return validation errors 
# ----------------------------------------------------------------
HierMultinom.path <- function(X, Y, groups, ngamma = 100, delta = 0.01, 
  lambda.vec = 10^seq(-3, 0, length=10), tol = 1e-8, 
  max.iter = 1e4, Xval, Yval, quiet=FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2] + 1
  K <- dim(Y)[2]
  X.mean <- colMeans(X)
  X.sd <- apply(X, 2, sd)
  Xtrain <- cbind(1, (X - tcrossprod(rep(1, n), X.mean))/(tcrossprod(rep(1, n), X.sd)))
  grad.temp <- crossprod(Xtrain, Y - rep(1, n)%*%t(colMeans(Y))/n)
  
  # ----------------------------------------------------
  # get candidate tuning parameters for gamma 
  # ----------------------------------------------------
  gamma.max <- max(apply(abs(grad.temp[-1,]), 1, function(x){sqrt(sum(x^2))}))/10
  gamma.min <- 1e-4*gamma.max
  gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=round(ngamma/2))
  margPi <- colSums(Y)/n
  beta.init0 <- matrix(0, p, K)
  beta.init0[1,] <- log(margPi) 
  gamma.matrix <- matrix(0, nrow=ngamma, ncol=length(lambda.vec))
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), length(gamma.vec)))
  
  for (kk in 1:length(lambda.vec)) {
    beta.init <- beta.init0
    for (jj in 1:length(gamma.vec)) {
      temp <- AccPGD.HierMult(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], gamma = gamma.vec[jj], max.iter = max.iter, tol = tol)
      if (sum(temp[-1,] == 0) < length(c(temp[-1,]))) {
        gamma.matrix[1,kk] <- gamma.vec[jj-1]
        break
      } else {
        beta.init <- temp
      }
    }
  }

  for (kk in 1:length(lambda.vec)) {
    gamma.matrix[,kk] <- 10^seq(log10(gamma.matrix[1,kk]), log10(delta*gamma.matrix[1,kk]), length=ngamma)
  }

  # ---------------------------------------------------- 
  # Compute complete solution path
  # ----------------------------------------------------
  margPi <- colSums(Y)/n
  beta.init0 <- matrix(0, p, K)
  beta.init0[1,] <- log(margPi) 
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), ngamma))
  val.errs <- matrix(Inf, nrow=length(lambda.vec), ncol=ngamma)
  
  for (kk in 1:length(lambda.vec)) {
    beta.init <- beta.init0
    for (jj in 1:ngamma) {
      temp <- AccPGD.HierMult(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], 
        gamma = gamma.matrix[jj,kk], max.iter = max.iter, tol = tol)
      beta.array[,,kk,jj] <- temp
      val.errs[kk,jj] <- validationLik(X, Xval, Yval, temp)
      if (jj > 10 & val.errs[kk,jj] < val.errs[kk,1]) {
        if (all(val.errs[kk,(jj):(jj-4)] > val.errs[kk,(jj-1):(jj-5)])) {
          break
        }
      }
      if (!quiet) {
        cat(val.errs[kk,jj], "\n")
      }
      beta.init <- temp
    }
    if (kk > 3) {
      if (all(val.errs[kk,] > min(val.errs[kk-1,])) & 
          all(val.errs[kk-1,] > min(val.errs[kk-2,])) & 
          all(val.errs[kk-2,] > min(val.errs[kk-3,]))) {
        break
      }
    }
  }
  
  ind1 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,1]
  ind2 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,2]

  return(list("beta.est" = beta.array[,,ind1,ind2], 
    "beta.array" = beta.array, 
    "val.errs" = val.errs,
    "X.train.mean" = X.mean, 
    "X.train.sd" = X.sd))
}



HierMultinom.predict <- function(mod.fit, Xtest, ind1 = NULL, ind2 = NULL){

  if (is.null(ind1)) {
    beta.est <- mod.fit$beta.est
  } else {
    beta.est <- mod.fit$beta.array[,,ind1, ind2]
  }
  XtestInput <- cbind(1, (Xtest - rep(1, dim(Xtest)[1])%*%t(mod.fit$X.train.mean))/(rep(1, dim(Xtest)[1])%*%t(mod.fit$X.train.sd)))
  if (any(mod.fit$X.train.sd==0)) {
    XtestInput[,which(mod.fit$X.train.sd==0)] <- 0
  }
  l0 <- XtestInput%*%beta.est
  prob.est <- exp(l0 - apply(l0, 1, max))/rowSums(exp(l0 - apply(l0, 1, max)))
  preds <- apply(prob.est, 1, which.max)
  return(list("probs" = prob.est, "preds" = preds))

}

HierMultinom.coef <- function(mod.fit, ind1 = NULL, ind2 = NULL){

  if (is.null(ind1)) {
    beta.est <- mod.fit$beta.est
  } else {
    beta.est <- mod.fit$beta.array[,,ind1, ind2]
  }
  beta.hat <- beta.est[-1,]/(mod.fit$X.train.sd%*%t(rep(1, dim(beta.est)[2])))
  intercept.hat <- beta.est[1,] - mod.fit$X.train.mean%*%beta.hat
  return(list("beta.hat" = beta.hat, "intercept.hat" = intercept.hat))

}

















proxFunction.overlap <- function(eta, lambda, groups, vmat = NULL){

  if (is.null(vmat)) {
    vmat <- matrix(0, length(eta), length(groups))
  }
  D.mats <- list()
  eo <- list()
  ginvs <- list()
  for(kk in 1:length(groups)){
    s <- length(groups[[kk]])
    D.mats[[kk]] <- diag(1, s) - matrix(1, s, s)/s
    eo[[kk]] <- eigen(D.mats[[kk]])
    ginvs[[kk]] <- ginv(D.mats[[kk]])
  }

  y.resid <- eta
  for (kk in length(groups):1) {
    y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] - crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
  }


  for (jj in 1:1e4) {
    vmatold <- vmat
    for (kk in length(groups):1) {
      y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] + crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
      tmp <- crossprod(D.mats[[kk]], y.resid[groups[[kk]]])
      vmat[groups[[kk]],kk] <- tmp * min(1, lambda/sqrt(sum(tmp^2)))
      y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] - crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
    }
    if (sum((vmatold - vmat)^2) < 1e-9){
      break
    }
  }

  t0 <- rep(0, length(eta))
  for (kk in length(groups):1) {
    t0[groups[[kk]]] <- t0[groups[[kk]]] + crossprod(D.mats[[kk]], vmat[groups[[kk]],kk])
  }
  return(list("result" = eta - t0, "vmat" = vmat))
}


# -------------------------------------------------
# Going to write a quick sub-algorithm to compute 
# -------------------------------------------------


outerProx.overlap <- function(y, lambda, gamma, groups, vmat){
  
  b.out <- rep(0, length(y))
  if (lambda!=0) {
    t0 <- proxFunction.overlap(y, lambda, groups, vmat)
    b.out <- t0$result
    vmat <- t0$vmat
  } else {
    b.out <- y
  }
  t0 <- sqrt(sum(b.out^2))
  if (t0 > gamma) {
    b.out <- (1 - gamma/t0)*b.out
  } else {
    b.out <- rep(0, length(y))
  }
  return(list("b.out" = b.out, "vmat" = vmat))
}

fullProx.overlap <- function(input, lambda, gamma, groups, vmat.full = NULL){
  
  out <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])
  if (is.null(vmat.full)) {
   vmat.full <- array(0, dim=c(dim(input)[1], dim(input)[2], length(groups))) 
  }
  out[1,] <- input[1,]
  for (j in 2:(dim(input)[1])) {
    if (sqrt(sum(input[j,]^2)) < gamma) {
      out[j,] <- 0
    } else {
      t0 <- outerProx.overlap(input[j,], lambda, gamma, groups, vmat.full[j,,])
      out[j,] <- t0$b.out
      vmat.full[j,,] <- t0$vmat
    }
  }
  return(list("out" = out, "vmat.full" = vmat.full))
}

















AccPGD.HierMultinomOverlap <- function(X, Y, groups, beta.init, lambda, gamma, max.iter = 1e4, tol = 1e-9){
  
  obj.vals <- rep(Inf, max.iter)
  obj.vals[1] <- getObj(X, Y, beta.init, groups, lambda, gamma)
  beta.iter <- beta.init
  beta.prev <- beta.iter
  alpha <- 1
  alpha.prev <- 1
  step.size <- 1e4*(dim(Y)[1]/(sqrt(dim(Y)[2])*sum(X^2)))
  lik.current <- getObj(X, Y, beta.iter, groups, lambda = 0, gamma = 0)
  L0 <- dim(Y)[1]/(sqrt(dim(Y)[2])*sum(X^2))
  vmat.full <- NULL
  
  for (kk in 1:max.iter) {
    
    search.point <- beta.iter + ((alpha.prev - 1)/(alpha))*(beta.iter - beta.prev)
    temp <- getGradLik(X, Y, search.point)
    temp.grad <- temp$grad
    lik.current <- temp$lik
    linesearch <- TRUE
    
    while (linesearch) {
      
      temp.step <- search.point - step.size*temp.grad
      t0 <- fullProx.overlap(temp.step, step.size*lambda, step.size*gamma, groups, vmat.full = vmat.full)
      beta.try <- t0$out
      vmat.full <- t0$vmat.full
      lik.try <- getObj(X, Y, beta.try, groups, lambda = 0, gamma = 0)
      check <- (lik.try <= lik.current + sum(temp.grad*(beta.try - search.point)) + (1/(2*step.size))*sum((beta.try - search.point)^2))

      if (check | step.size == L0) {
        linesearch <- FALSE
        beta.prev <- beta.iter
        beta.iter <- beta.try
        alpha.prev <- alpha
        alpha <- (1 + sqrt(1 + 4*alpha.prev^2))/2
      } else {
        #cat(check, "\n")
        step.size <- max(step.size/2, L0)
      }
    }
    
    obj.vals[kk+1] <- lik.try + evalPen(beta.iter, lambda, gamma, groups)
    if (kk > 1) {
      if (obj.vals[kk+1] > obj.vals[kk]) {
        beta.iter <- beta.prev
        obj.vals[kk+1] <- obj.vals[kk]
      }
    }
    if (kk > 10) {
      if (all(abs((obj.vals[kk:(kk-5)] - obj.vals[(kk+1):(kk-4)])) < tol * abs(obj.vals[1]))) {
        break
      }
    }
  }
  return(beta.iter)
}

HierMultinomOverlap.path <- function(X, Y, groups, ngamma = 100, delta = 0.01, 
  lambda.vec = 10^seq(-3, 0, length=10), 
  tol = 1e-8, max.iter = 1e4, Xval, Yval, quiet=FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2] + 1
  K <- dim(Y)[2]
  x.sd.temp <- apply(X, 2, sd)
  Xtrain <- cbind(1, (X - tcrossprod(rep(1, n), colMeans(X)))/(tcrossprod(rep(1, n), x.sd.temp)))
  if (any(x.sd.temp==0)) {
    Xtrain[,which(x.sd.temp==0)] <- 0
  }
  grad.temp <- crossprod(Xtrain, Y - rep(1, n)%*%t(colMeans(Y))/n)
  gamma.max <- max(apply(abs(grad.temp[-1,]), 1, function(x){sqrt(sum(x^2))}))/10
  gamma.min <- 1e-4*gamma.max
  gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=ngamma)
  margPi <- colSums(Y)/n
  beta.init0 <- matrix(0, p, K)
  beta.init0[1,] <- log(margPi) 
  gamma.matrix <- matrix(0, nrow=ngamma, ncol=length(lambda.vec))
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), length(gamma.vec)))
  for (kk in 1:length(lambda.vec)) {
    beta.init <- beta.init0
    for (jj in 1:length(gamma.vec)) {
      temp <- AccPGD.HierMultinomOverlap(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], gamma = gamma.vec[jj], max.iter = max.iter, tol = tol)
      if (sum(temp[-1,] == 0) < length(c(temp[-1,]))) {
        gamma.matrix[1,kk] <- gamma.vec[jj-1]
        gamma.vec <- 10^seq(log10(gamma.matrix[1,kk]), log10(delta*gamma.matrix[1,kk]), length=ngamma)
        break
      } else {
        beta.init <- temp
      }
    }
  }

  for (kk in 1:length(lambda.vec)) {
    gamma.matrix[,kk] <- 10^seq(log10(gamma.matrix[1,kk]), log10(delta*gamma.matrix[1,kk]), length=ngamma)
  }

  # ----------------------------------------
  # Run one pass to get appropriate TP grid
  # ----------------------------------------
  margPi <- colSums(Y)/n
  beta.init0 <- matrix(0, p, K)
  beta.init0[1,] <- log(margPi) 
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), ngamma))
  val.errs <- matrix(Inf, nrow=length(lambda.vec), ncol=ngamma)
  for (kk in 1:length(lambda.vec)) {
    beta.init <- beta.init0
    for (jj in 1:ngamma) {
      temp <- AccPGD.HierMultinomOverlap(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], gamma = gamma.matrix[jj,kk], max.iter = max.iter, tol = tol)
      beta.array[,,kk,jj] <- temp
      val.errs[kk,jj] <- validationLik(X, Xval, Yval, temp)
      if (jj > 10 & val.errs[kk,jj] < val.errs[kk,1]) {
        if (all(val.errs[kk,(jj):(jj-3)] > val.errs[kk,(jj-1):(jj-4)])) {
          break
        }
      }
      if (!quiet) {
        cat(val.errs[kk,jj], "\n")
      }
      beta.init <- temp
    }
    if (kk > 3) {
      if (all(val.errs[kk,] > min(val.errs[kk-1,])) & 
          all(val.errs[kk-1,] > min(val.errs[kk-2,])) & 
          all(val.errs[kk-2,] > min(val.errs[kk-3,]))) {
        break
      }
    }
  }
  ind1 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,1]
  ind2 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,2]

  return(list("beta.est" = beta.array[,,ind1,ind2], 
    "val.errs" = val.errs, 
    "beta.array" = beta.array, 
    "X.train.mean" = apply(X, 2, mean), 
    "X.train.sd" = apply(X, 2, sd)))
}
