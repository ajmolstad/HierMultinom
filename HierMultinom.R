
getGrad_Lik <- function(X, Y, beta){
  t0 <- exp(crossprod(t(X), beta))
  P <- t0/tcrossprod(rowSums(t0), rep(1, dim(beta)[2]))
  t0 <- -sum(log(rowSums(Y*P)))/dim(Y)[1]
  grad <- crossprod(X, P - Y)/dim(Y)[1]
  return(list("grad" = grad, "lik" = t0))
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
  n <- dim(X)[1]
  lleval <- crossprod(t(X), beta)
  t0 <- -sum(log(rowSums(Y*exp(lleval))) - log(rowSums(exp(lleval))))/n
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

prox_op_1 <- function(y, lambda) {
  
  s <- length(y)
  tildeW <- diag(1, s) - matrix(1, s, s)/s
  v <- crossprod(svd(tildeW)$u, y)[-s]
  t1 <- crossprod(tildeW, y)
  t0 <- sum(t1^2)
  if(t0 <= lambda^2){
    result <- mean(y)*rep(1, s)
  } else {
    result <- y - crossprod(tildeW, y)/(sqrt(sum(v^2))/lambda)
  }
  return(result)

  # -----------------------------------------------
  # compare to CVXR
  # -------------------------------------------------
  # out <- result
  # 0.5*sum((out - y)^2) + lambda*sqrt(sum((out - mean(out))^2))
  # s <- length(y)
  # b <- Variable(s)
  # obj <- .5*sum((b - y)^2) + lambda*norm2(b - mean(b))
  # prob <- Problem(Minimize(obj))
  # result <- solve(prob)
  # 0.5*sum((result[[1]] - y)^2) + lambda*sqrt(sum((result[[1]] - mean(result[[1]]))^2))

}

outer_prox <- function(y, lambda, gamma, groups){
  
  b.out <- rep(0, length(y))
  if(lambda!=0){
    for(jj in 1:length(groups)){
      b.out[groups[[jj]]] <- prox_op_1(y[groups[[jj]]], lambda)
    }
  } else {
    b.out <- y
  }
  t0 <- sqrt(sum(b.out^2))
  if(t0 > gamma){
    b.out <- (1 - gamma/t0)*b.out
  } else {
    b.out <- rep(0, length(y))
  }
  return(b.out)
  
}


full_prox <- function(input, lambda, gamma, groups){
  out <- matrix(0, nrow=dim(input)[1], ncol=dim(input)[2])
  out[1,] <- input[1,]
  # need to screen first, then use this function. 
  for(j in 2:dim(input)[1]){
    if(sqrt(sum(input[j,]^2)) < gamma){
      out[j,] <- 0
    } else {
      out[j,] <- outer_prox(input[j,], lambda, gamma, groups)
    }
  }
  return(out)
}

validationLik <- function(X, Xval, Yval, beta){
  Xval.temp <- cbind(1, (Xval - tcrossprod(rep(1, dim(Xval)[1]), colMeans(X)))/tcrossprod(rep(1, dim(Xval)[1]), apply(X, 2, sd)))
  lleval <- crossprod(t(Xval.temp),beta)
  probs <- exp(lleval)/tcrossprod(rowSums(exp(lleval)), rep(1, dim(beta)[2]))
  return(-2*sum(log(rowSums(Yval*probs))))
}

# # ======================================
# #
# # ======================================
# max.iter <- 2000
# lambda <- 0.00
# gamma <- 0.01
# tol <- 1e-8
# beta.init <- matrix(0, nrow = p, ncol = K)



AccPGD_HeirMult <- function(X, Y, groups, beta.init, lambda, gamma, max.iter = 1e4, tol = 1e-9){
  
  obj.vals <- rep(Inf, max.iter)
  obj.vals[1] <- getObj(X, Y, beta.init, groups, lambda, gamma)
  beta.iter <- beta.init
  beta.prev <- beta.iter
  alpha <- 1
  alpha.prev <- 1
  L0 <- (dim(Y)[1])/(sqrt(dim(Y)[2])*sum(X^2))
  step.size <- 1e4*L0
  lik.current <- getObj(X, Y, beta.iter, groups, lambda = 0, gamma = 0)
  for(kk in 1:max.iter){
    
    search.point <- beta.iter + ((alpha.prev - 1)/(alpha))*(beta.iter - beta.prev)
    temp <- getGrad_Lik(X, Y, search.point)
    temp.grad <- temp$grad
    lik.current <- temp$lik
    linesearch <- TRUE
    
    while(linesearch){
      
      temp.step <- search.point - step.size*temp.grad
      beta.try <- full_prox(temp.step, step.size*lambda, step.size*gamma, groups)
      lik.try <- getObj(X, Y, beta.try, groups, lambda = 0, gamma = 0)
      check <- (lik.try <= lik.current + sum(diag(crossprod(temp.grad, beta.try - search.point))) + (1/(2*step.size))*sum((beta.try - search.point)^2))
      if(check | step.size == L0){
        linesearch <- FALSE
        beta.prev <- beta.iter
        beta.iter <- beta.try
        alpha.prev <- alpha
        alpha <- (1 + sqrt(1 + 4*alpha.prev^2))/2
        #lik.current <- lik.try
      } else {
        step.size <- max(step.size*.5, L0)
      }
    }
    
    # cat(step.size, "\n")
    obj.vals[kk+1] <- lik.try + evalPen(beta.iter, lambda, gamma, groups)
    cat(obj.vals[kk+1],"\n")
    if(kk > 10){
      # cat(all(abs(obj.vals[kk:(kk-3)] - obj.vals[(kk+1):(kk-2)]) < tol * abs(obj.vals[1])), "\n")
      if(all(abs(obj.vals[kk:(kk-3)] - obj.vals[(kk+1):(kk-2)]) < tol * abs(obj.vals[1]))){
        break
      }
      if(kk > max.iter){
        break
      }
    }
  }
  return(beta.iter)
}

HeirMult.Path <- function(X, Y, groups, ngamma = 100, delta = 0.01, lambda.vec = 10^seq(-3, 0, length=10), tol = 1e-8, max.iter = 1e4, Xval, Yval){
  
  n <- dim(X)[1]
  p <- dim(X)[2] + 1
  K <- dim(Y)[2]
  Xtrain <- cbind(1, (X - tcrossprod(rep(1, n), colMeans(X)))/(tcrossprod(rep(1, n), apply(X, 2, sd))))
  grad.temp <- crossprod(Xtrain, Y - rep(1, n)%*%t(colMeans(Y))/n)
  gamma.max <- max(apply(abs(grad.temp[-1,]), 1, function(x){sqrt(sum(x^2))}))/10
  gamma.min <- 1e-4*gamma.max
  gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=ngamma)
  margPi <- colSums(Y)/n
  beta.init0 <- matrix(0, p, K)
  beta.init0[1,] <- log(margPi) 
  gamma.matrix <- matrix(0, nrow=ngamma, ncol=length(lambda.vec))
  beta.array <- array(0, dim=c(p, K, length(lambda.vec), length(gamma.vec)))
  for(kk in 1:length(lambda.vec)){
    beta.init <- beta.init0
    for(jj in 1:length(gamma.vec)){
      temp <- AccPGD_HeirMult(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], gamma = gamma.vec[jj], max.iter = max.iter, tol = tol)
      #cat(sum(temp[-1,] == 0), "\n")
      if(sum(temp[-1,] == 0) < length(c(temp[-1,]))){
        gamma.matrix[1,kk] <- gamma.vec[jj-1]
        break
      } else {
        beta.init <- temp
      }
    }
    #cat("# --------", kk,"\n")
  }



  for(kk in 1:length(lambda.vec)){
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
  for(kk in 1:length(lambda.vec)){
    beta.init <- beta.init0
    for(jj in 1:ngamma){
      temp <- AccPGD_HeirMult(Xtrain, Y, groups, beta.init, lambda = lambda.vec[kk], gamma = gamma.matrix[jj,kk], max.iter = max.iter, tol = tol)
      beta.array[,,kk,jj] <- temp
      val.errs[kk,jj] <- validationLik(X, Xval, Yval, temp)
      if(jj > 10 & val.errs[kk,jj] < val.errs[kk,1]){
        if(all(val.errs[kk,(jj):(jj-8)] > val.errs[kk,(jj-1):(jj-9)])){
          break
        }
      }
      cat(val.errs[kk,jj], "\n")
      beta.init <- temp
    }
    if(kk > 3){
      if(all(val.errs[kk,] > min(val.errs[kk-1,])) & all(val.errs[kk-1,] > min(val.errs[kk-2,])) & all(val.errs[kk-2,] > min(val.errs[kk-3,]))){
        break
      }
    }
  }
  ind1 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,1]
  ind2 <- which(val.errs == min(val.errs), arr.ind=TRUE)[1,2]
  return(beta.array[,,ind1,ind2])
}