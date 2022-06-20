expand_parameters = function(run_name,
                             considered_values,
                             defaults,
                             n_replicates,
                             methods) {

  parameter_list = list()

  for (name in names(considered_values)) {
    values = considered_values[[name]]
    for (replicate in 1:n_replicates) {
      for (value in values) {
        for (method in methods) {
          params = c(
            defaults,
            experiment = name,
            replicate = replicate,
            method = method,
            run = run_name
          )
          params[name] = value
          parameter_list = c(parameter_list, list(params))
        }
      }
    }
  }

  return(parameter_list)
}

prepare_real_data_application = function(n_genes,
                                         n_train,
                                         n_validation,
                                         n_test,
                                         cache_path,
                                         replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  prepare_hao_2020(cache_path, n_genes, n_train, n_validation, n_test)

}

prepare_hao_2020 = function(cache_path, n_genes = NA, n_train = NA, n_validation = NA, n_test = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_hao_2020(cache_path)

  SingleCellExperiment::altExp(data) = NULL
  SingleCellExperiment::counts(data) = NULL

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = ifelse(data$cell_type_2 == "Treg", data$cell_type_3, data$cell_type_2)

  removed_labels = "*Proliferating*"
  data = data[, !grepl(removed_labels, data$cell_type)]

  if (sce) return(data)

  data = data[, sample(1:ncol(data), n_train + n_validation + n_test, replace = FALSE)]

  X = t(as.matrix(SingleCellExperiment::logcounts(data)))
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

  list(
    X_train = X[train_indices, ],
    X_validation = X[validation_indices, ],
    X_test = X[test_indices, ],
    Y_train = Y[train_indices, ],
    Y_validation = Y[validation_indices, ],
    Y_test = Y[test_indices, ],
    groups = groups
  )

}

fit_Multires = function(data) {

  t0 <- HierMultNestedOverlap.Path(data$X_train, data$Y_train, data$groups, ngamma = 100, delta = 1e-04,
                                   lambda.vec = c(0, 10^seq(-5,-3, length=5)), tol = 1e-7,
                                   max.iter = 1e4, data$X_validation, data$Y_validation)

  return(list(fit = t0$beta, Beta_hat = t0$beta[-1, ], alpha_hat = t0$beta[1, ], val_error = t0$val.errs))

}

fit_grouped = function(data) {

  t0 <- HierMultNestedOverlap.Path(data$X_train, data$Y_train, data$groups, ngamma = 100, delta = 1e-04,
                                   lambda.vec = 0, tol = 1e-7,
                                   max.iter = 1e4, data$X_validation, data$Y_validation)

  return(list(fit = t0$beta, Beta_hat = t0$beta[-1, ], alpha_hat = t0$beta[1, ], val_error = t0$val.errs))

}

fit_Group = function(data) {

  temp <- glmnet(y = data$Y_train, x = data$X_train, family="multinomial", type.multinomial = "grouped", trace.it = 1, maxit = 1e7, thresh = 1e-9, nlambda = 100)
  val.err <- rep(0, length(temp$lambda))
  for(jj in 1:length(temp$lambda)){
    preds <- predict(temp, s=temp$lambda[jj], newx = data$X_validation, type="response")[,,1]
    val.err[jj] <- -2*sum(log(rowSums(data$Y_validation*preds)))
    print(val.err[jj])
  }
  lambda = which(val.err == min(val.err, na.rm = TRUE))
  attributes(temp)$selected_lambda = lambda

  print(temp$lambda)

  return(list(fit = temp, Beta_hat = do.call(cbind, lapply(temp$beta[colnames(data$Y_train)], function(x) x[, lambda])), alpha_hat = temp$a0[colnames(data$Y_train), lambda]))

}

fit_L1 = function(data) {

  temp <- glmnet(y = data$Y_train, x = data$X_train, family="multinomial", type.multinomial = "ungrouped", trace.it = 1, maxit = 1e7, thresh = 1e-9, nlambda = 100)
  val.err <- rep(0, length(temp$lambda))
  for(jj in 1:length(temp$lambda)){
    preds <- predict(temp, s=temp$lambda[jj], newx = data$X_validation, type="response")[,,1]
    val.err[jj] <- -2*sum(log(rowSums(data$Y_validation*preds)))
    print(val.err[jj])
  }
  lambda = which(val.err == min(val.err, na.rm = TRUE))
  attributes(temp)$selected_lambda = lambda

  print(temp$lambda)

  return(list(fit = temp, Beta_hat = do.call(cbind, lapply(temp$beta[colnames(data$Y_train)], function(x) x[, lambda])), alpha_hat = temp$a0[colnames(data$Y_train), lambda]))

}

compute_performance = function(fit, Y_test, X_test, X_train) {

  if (any(class(fit) == "glmnet")) {

    prob.est <- predict(fit, s = fit$lambda[attributes(fit)$selected_lambda], newx = X_test, type="response")[,,1]

  } else {

    x.sd.temp <- apply(X_train, 2, sd)
    X_test <- cbind(1, (X_test - rep(1, dim(X_test)[1])%*%t(apply(X_train, 2, mean)))/(rep(1, dim(X_test)[1])%*%t(apply(X_train, 2, sd))))
    if(any(x.sd.temp==0)){
      X_test[,(which(x.sd.temp==0) + 1)] <- 0
    }

    l0 <- exp(X_test%*%fit)
    prob.est <- l0/rowSums(l0)

  }

  classification.error <- sum(apply(Y_test, 1, which.max) != apply(prob.est, 1, which.max))/dim(Y_test)[1]
  deviance = -2*sum(log(rowSums(Y_test*prob.est)))

  return(list(error = classification.error, deviance = deviance))

}

evaluate_parameters = function(parameters) {

  method = parameters$method

  method_function = get(paste0("fit_", method))

  print("Preparing data")

  data = do.call(prepare_real_data_application, parameters[intersect(names(parameters), formalArgs(prepare_real_data_application))])
  data$parameters = parameters

  print("Fitting model")

  fit = method_function(data)
  colnames(fit$Beta_hat) = colnames(data$Y_train)
  names(fit$alpha_hat) = colnames(data$Y_train)

  print("Evaluating performance")

  performance = compute_performance(fit$fit,
                                    data$Y_test,
                                    data$X_test,
                                    data$X_train)

  result = list(parameters = parameters, performance = performance, Beta_hat = fit$Beta_hat, alpha_hat = fit$alpha_hat, val_error = fit$val_error)

  return(result)

}
