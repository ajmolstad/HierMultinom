library(MASS)
library(Rcpp)

getGrad_Lik <- function(X, Y, beta){
  t0 <- eigenMatMult(X, beta)
  if(any(t0 > 700)){ #710 is exact threshold for Inf
    t1 <- exp(t0 - apply(t0, 1, max))
    P <- t1/eigenMatMult(armaRowSums(t1), t(rep(1, dim(beta)[2])))
  } else {
    P <- exp(t0)/eigenMatMult(armaRowSums(exp(t0)), t(rep(1, dim(beta)[2])))
  }
  grad <- eigenMatMult(t(X), P - Y)/dim(Y)[1]
  return(list("grad" = grad, "lik" = -sum(log(armaRowSums(Y*P)))/dim(Y)[1]))
}


evalPen <- function(beta, lambda, gamma, groups){
  p1 <- gamma*sum(apply(beta[-1,], 1, function(x){sqrt(sum(x^2))}))
  p2 <- 0
  for(kk in 1:length(groups)){
    p2 <- p2 + sum(apply(beta[-1,groups[[kk]]] - rowMeans(beta[-1,groups[[kk]], drop=FALSE])%*%t(rep(1, length(groups[[kk]]))), 1, function(x){sqrt(sum(x^2))}))
  }
  return(p1 + lambda*p2)
}


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


# Psi  <- function(eta, w, group) {
  
#   y <- eta[group]
#   s <- length(y)
#   tildey <- y - mean(y)
#   tildeW <- diag(1, s) - matrix(1, s, s)/s
#   v <- crossprod(svd(tildeW)$u, tildey)
#   t1 <- crossprod(tildeW, tildey)
#   t0 <- sum(t1^2)
#   if (t0 <= w^2) {
#     result <- mean(y)*rep(1, s)
#   } else {
#     theta <- solve(crossprod(tildeW) + diag(w/(sqrt(sum(v[-s]^2)) - w), s), t1)
#     result <- theta + mean(y - theta)
#   }
#   return(result)
# }


prox.function <- function(eta, lambda, groups, D.input, vmat = NULL){

  if(is.null(vmat)){
    vmat <- matrix(0, length(eta), length(groups))
  }
  D.mats <- list()
  for(kk in 1:length(groups)){
    s <- length(groups[[kk]])
    D.mats[[kk]] <- diag(1, s) - matrix(1, s, s)/s
  }

  y.resid <- eta
  for(kk in length(groups):1){
    y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] - crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
  }

  for(jj in 1:1e4){
    
    vmatold <- vmat
    for(kk in length(groups):1){
      y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] + crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
      vmat[groups[[kk]],kk] <- crossprod(ginv(D.mats[[kk]]), crossprod(D.mats[[kk]], y.resid[groups[[kk]]]))
      if(sum(vmat[groups[[kk]],kk]^2) < lambda^2){
        y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] - crossprod(D.mats[[kk]],vmat[groups[[kk]], kk])
        #cat("passed", "\n")
      } else {
        eo <- eigen(D.mats[[kk]])
        w <- crossprod(y.resid[groups[[kk]]], eo$vec[,which(abs(eo$values) > 1e-8)])
        t0 <- sum(w^2)
        tau <- sqrt(t0/lambda^2) - 1 
        #cat("# ------------------", "\n")
        #tryFunc <- function(x){c(crossprod(solve(crossprod(D.mats[[kk]]) + diag(x, dim(D.mats[[kk]])[1]), crossprod(D.mats[[kk]], y.resid[groups[[kk]]])))) - lambda^2}
        #out <- uniroot(tryFunc, c(1e-12,1e12), tol=1e-8)
        #cat(tau, out$root, "\n")
        # Replace with exact solution 
        vmat[groups[[kk]],kk] <- crossprod(t(ginv(crossprod(D.mats[[kk]]) + diag(tau, dim(D.mats[[kk]])[1]))), crossprod(D.mats[[kk]], y.resid[groups[[kk]]]))
        y.resid[groups[[kk]]] <- y.resid[groups[[kk]]] - crossprod(D.mats[[kk]], vmat[groups[[kk]], kk])
      }
    }
    #if(jj > 2){
      if(sum((vmatold - vmat)^2) < 1e-9){
        #cat("Broken at iteration ", jj,":" "\n")
        break
      }
    #}
  }

  t0 <- rep(0, length(eta))
  for(kk in length(groups):1){
    t0[groups[[kk]]] <- t0[groups[[kk]]] + crossprod(D.mats[[kk]], vmat[groups[[kk]],kk])
  }
  return(list("result" = eta - t0, "vmat" = vmat))
}


# -------------------------------------------------
# Going to write a quick sub-algorithm to compute 
# -------------------------------------------------


outer_prox <- function(y, lambda, gamma, groups, vmat){
  
  b.out <- rep(0, length(y))
  if(lambda!=0){
    t0 <- prox.function(y, lambda, groups, vmat)
    b.out <- t0$result
    vmat <- t0$vmat
  } else {
    b.out <- y
  }
  t0 <- sqrt(sum(b.out^2))
  if(t0 > gamma){
    b.out <- (1 - gamma/t0)*b.out
  } else {
    b.out <- rep(0, length(y))
  }
  return(list("b.out" = b.out, "vmat" = vmat))
}

full_prox <- function(input, lambda, gamma, groups, vmat.full = NULL){
  
  out <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])
  if(is.null(vmat.full)){
   vmat.full <- array(0, dim=c(dim(input)[1], dim(input)[2], length(groups))) 
  }
  out[1,] <- input[1,]
  for(j in 2:(dim(input)[1])){
    if(sqrt(sum(input[j,]^2)) < gamma){
      out[j,] <- 0
    } else {
      t0 <- outer_prox(input[j,], lambda, gamma, groups, vmat.full[j,,])
      out[j,] <- t0$b.out
      vmat.full[j,,] <- t0$vmat
    }
  }
  return(list("out" = out, "vmat.full" = vmat.full))
}

validationLik <- function(X, Xval, Yval, beta){
  x.sd.temp <- apply(X, 2, sd)
  Xval.temp <- cbind(1, (Xval - tcrossprod(rep(1, dim(Xval)[1]), colMeans(X)))/tcrossprod(rep(1, dim(Xval)[1]), x.sd.temp))
  if(any(x.sd.temp==0)){
    Xval.temp[,which(x.sd.temp==0)] <- 0
  }
  lleval <- crossprod(t(Xval.temp),beta)
  probs <- exp(lleval - apply(lleval, 1, max))/tcrossprod(rowSums(exp(lleval - apply(lleval, 1, max))), rep(1, dim(beta)[2]))
  return(-2*sum(log(rowSums(Yval*probs))))
}



AccPGD_HierMultNestedOverlap <- function(X, Y, groups, beta.init, lambda, gamma, max.iter = 1e4, tol = 1e-9){
  
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
  
  for(kk in 1:max.iter){
    
    search.point <- beta.iter + ((alpha.prev - 1)/(alpha))*(beta.iter - beta.prev)
    temp <- getGrad_Lik(X, Y, search.point)
    temp.grad <- temp$grad
    lik.current <- temp$lik
    linesearch <- TRUE
    
    while(linesearch){
      
      temp.step <- search.point - step.size*temp.grad
      t0 <- full_prox(temp.step, step.size*lambda, step.size*gamma, groups, vmat.full = vmat.full)
      beta.try <- t0$out
      vmat.full <- t0$vmat.full
      lik.try <- getObj(X, Y, beta.try, groups, lambda = 0, gamma = 0)
      check <- (lik.try <= lik.current + sum(temp.grad*(beta.try - search.point)) + (1/(2*step.size))*sum((beta.try - search.point)^2))

      if(check | step.size == L0){
        linesearch <- FALSE
        beta.prev <- beta.iter
        beta.iter <- beta.try
        alpha.prev <- alpha
        alpha <- (1 + sqrt(1 + 4*alpha.prev^2))/2
        #lik.current <- lik.try
      } else {
        #cat(check, "\n")
        step.size <- max(step.size/2, L0)
      }
    }
    
    # cat(step.size, "\n")

    obj.vals[kk+1] <- lik.try + evalPen(beta.iter, lambda, gamma, groups)
    if(kk > 1){
      if(obj.vals[kk+1] > obj.vals[kk]){
        beta.iter <- beta.prev
        obj.vals[kk+1] <- obj.vals[kk]
      }
    }
    # cat(obj.vals[kk+1],"\n")
    if(kk > 10){
      if(all(abs((obj.vals[kk:(kk-5)] - obj.vals[(kk+1):(kk-4)])) < tol * abs(obj.vals[1]))){
        break
      }
    }
  }
  return(beta.iter)
}

HierMultinomOverlapZero.path <- function(X, Y, groups, lambda.vec = 10^seq(-3, 0, length=10), tol = 1e-8, max.iter = 1e4, Xval, Yval){
  
  n <- dim(X)[1]
  p <- dim(X)[2] + 1
  K <- dim(Y)[2]
  x.sd.temp <- apply(X, 2, sd)
  Xtrain <- cbind(1, (X - tcrossprod(rep(1, n), colMeans(X)))/(tcrossprod(rep(1, n), x.sd.temp)))
  if(any(x.sd.temp==0)){
    Xtrain[,which(x.sd.temp==0)] <- 0
  }
  gamma <- 0 
  margPi <- colSums(Y)/n
  gamma.matrix <- matrix(0, nrow=1, ncol=length(lambda.vec))
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), 1))

  # ----------------------------------------
  # Run one pass to get appropriate TP grid
  # ----------------------------------------
  margPi <- colSums(Y)/n
  beta.init <- matrix(0, p, K)
  beta.init[1,] <- log(margPi) 
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), 1))
  val.errs <- matrix(Inf, nrow=length(lambda.vec), ncol=1)
  for(jj in 1:length(lambda.vec)){
    temp <- AccPGD_HierMultNestedOverlap(Xtrain, Y, groups, beta.init, lambda = lambda.vec[jj], gamma = 0, max.iter = max.iter, tol = tol)
    beta.array[,,jj,1] <- temp
    val.errs[jj,1] <- validationLik(X, Xval, Yval, temp)
    cat(val.errs[jj,1], "\n")
    if(jj > 10 & val.errs[jj,1] < val.errs[1,1]){
      if(all(val.errs[(jj):(jj-4),1] > val.errs[(jj-1):(jj-5),1])){
        break
      }
    }
    beta.init <- temp
  }
  ind1 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,1]
  ind2 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,2]
  return(list("beta" = beta.array[,,ind1,ind2], "val.errs" = val.errs))
}