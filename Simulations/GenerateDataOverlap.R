# -------------------------------------------
# Generate data function
# -------------------------------------------
if(Model == 1){
  beta <- matrix(0, nrow=p, K)
  fine.inds <- sample(1:p, 18)
  for(kk in c(1:2, 4:5)){
    for(ll in 1:length(fine.inds)){
      beta[fine.inds[ll], groups[[kk]]] <- rnorm(length(groups[[kk]]), sd = sqrt(5))
    }
  }
} 


if(Model > 1){
  s <- 3 + (Model - 2)*3
  beta <- matrix(0, nrow=p, K)
  coarse.inds <- sample(1:p, s)
  for(kk in 1:3){
    t0 <- rnorm(length(coarse.inds), sd = sqrt(5))
    beta[coarse.inds, groups[[kk]]] <- t0%*%t(rep(1, length(groups[[kk]])))
  }
  fine.inds <- sample(c(1:p)[-coarse.inds], 18-s)
  for(kk in c(1:2, 4:5)){
    for(ll in 1:length(fine.inds)){
      beta[fine.inds[ll], groups[[kk]]] <- rnorm(length(groups[[kk]]), sd = sqrt(5))
    }
  }
} 


t0 <- exp(X%*%beta)
probs <- t0/rowSums(t0)%*%t(rep(1, K))
Y <- matrix(0, nrow = n, ncol = K)
for(jj in 1:n){
  Y[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}

t0 <- exp(Xval%*%beta)
probs <- t0/rowSums(t0)%*%t(rep(1, K))
Yval <- matrix(0, nrow = n, ncol = K)
for(jj in 1:n){
  Yval[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}


t0 <- exp(Xtest%*%beta)
probs <- t0/rowSums(t0)%*%t(rep(1, K))
Ytest <- matrix(0, nrow = 1e4, ncol = K)
for(jj in 1:1e4){
  Ytest[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}