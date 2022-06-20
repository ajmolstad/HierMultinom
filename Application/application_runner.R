library(glmnet)
source("scripts/HierMultinomNestedOverlap.R")
source("scripts/application_setup.R")

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])
RESULT_PATH = "results/application_Multires_Group_L1"
dir.create(RESULT_PATH, recursive = TRUE)

parameters = expand_parameters(
  run_name = "real_data_application",
  considered_values = list(
    n_train = 1:5 * 10000,
    n_genes = c(500, 750, 1000, 1500, 2000)
  ),
  defaults = list(
    cache_path = "../AnnotatedPBMC/data",
    n_train = 20000,
    n_validation = 20000,
    n_test = 20000,
    n_genes = 1000
  ),
  50,
  c("Multires", "Group", "L1")
)

current_parameters = parameters[[PARAMETER_ID]]

print(PARAMETER_ID)
print(current_parameters)

system.time({result = evaluate_parameters(current_parameters)})
saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))
