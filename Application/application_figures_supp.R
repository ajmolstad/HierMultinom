library(tidyverse)
library(patchwork)

almost.unique =
  function(x,
           tolerance = sqrt(.Machine$double.eps),
           ...)
  {
    y <- round(x/tolerance, 0)
    d <- duplicated(y, ...)
    x[!d]
  }

RESULT_PATH = "results/application_Multires_Group_L1"
FIGURES_PATH = file.path(RESULT_PATH, "figures")
sapply(FIGURES_PATH, function(path) dir.create(path, recursive = TRUE))

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, function(x) {
  tryCatch({
    print(x)
    result = readRDS(x)
    data.frame(
      run = result$parameters$run,
      experiment = result$parameters$experiment,
      n_train = result$parameters$n_train,
      n_genes = result$parameters$n_genes,
      value = result$parameters[[result$parameters$experiment]],
      method = result$parameters$method,
      replicate = result$parameters$replicate,
      error = result$performance$error,
      deviance = result$performance$deviance,
      nonzero = sum(apply(result$Beta_hat, 1, function(x) all.equal(setNames(x, nm = NULL), rep(0, ncol(result$Beta_hat))) != TRUE)),
      df = sum(apply(result$Beta_hat, 1, function(x) length(setdiff(almost.unique(x), 0))))
    )
  },
  error = function(x) NULL)
})

result = do.call(rbind, results)

result = result %>% group_by(run, experiment, value, method, n_train, n_genes) %>% summarize(
  mean_error = mean(error),
  mean_deviance = mean(deviance),
  se_error = sd(error) / sqrt(n()),
  se_deviance = sd(deviance) / sqrt(n()),
  mean_nonzero = mean(nonzero),
  mean_df = mean(df),
  se_nonzero = sd(nonzero) / sqrt(n()),
  se_df = sd(df) / sqrt(n())
)

methods = c("Group", "L1", "Multires")
plasma_pal = c(2, 3, 1)
names(plasma_pal) = methods

plasma_pal = plasma_pal[methods]
result = result %>% filter(method %in% methods)
result$method = factor(result$method, levels = methods)

p1 = ggplot(
  result[result$experiment == "n_genes", ],
  aes(
    x = n_genes,
    y = mean_error,
    group = method,
    linetype = method,
    ymin = mean_error - 2 * se_error,
    ymax = mean_error + 2 * se_error
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("p") +
  ylab("Error rate") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p2 = ggplot(
  result[result$experiment == "n_genes", ],
  aes(
    x = n_genes,
    y = mean_deviance,
    group = method,
    linetype = method,
    ymin = mean_deviance - 2 * se_deviance,
    ymax = mean_deviance + 2 * se_deviance
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("p") +
  ylab("Deviance") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p3 = ggplot(
  result[result$experiment == "n_genes", ],
  aes(
    x = n_genes,
    y = mean_df,
    group = method,
    linetype = method,
    ymin = mean_df - 2 * se_df,
    ymax = mean_df + 2 * se_df
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("p") +
  ylab("Degrees of freedom") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p7 = ggplot(
  result[result$experiment == "n_genes", ],
  aes(
    x = n_genes,
    y = mean_nonzero,
    group = method,
    linetype = method,
    ymin = mean_nonzero - 2 * se_nonzero,
    ymax = mean_nonzero + 2 * se_nonzero
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("p") +
  ylab("Nonzero genes") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p4 = ggplot(
  result[result$experiment == "n_train", ],
  aes(
    x = n_train / 1000,
    y = mean_error,
    group = method,
    linetype = method,
    ymin = mean_error - 2 * se_error,
    ymax = mean_error + 2 * se_error
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("n/1000") +
  ylab("Error rate") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p5 = ggplot(
  result[result$experiment == "n_train", ],
  aes(
    x = n_train / 1000,
    y = mean_deviance,
    group = method,
    linetype = method,
    ymin = mean_deviance - 2 * se_deviance,
    ymax = mean_deviance + 2 * se_deviance
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("n/1000") +
  ylab("Deviance") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p6 = ggplot(
  result[result$experiment == "n_train", ],
  aes(
    x = n_train / 1000,
    y = mean_df,
    group = method,
    linetype = method,
    ymin = mean_df - 2 * se_df,
    ymax = mean_df + 2 * se_df
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("n/1000") +
  ylab("Degrees of freedom") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

p8 = ggplot(
  result[result$experiment == "n_train", ],
  aes(
    x = n_train / 1000,
    y = mean_nonzero,
    group = method,
    linetype = method,
    ymin = mean_nonzero - 2 * se_nonzero,
    ymax = mean_nonzero + 2 * se_nonzero
  )
) +
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  xlab("n/1000") +
  ylab("Nonzero genes") +
  labs(linetype = "Method") +
  theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = plasma_pal) +
  theme(plot.title = element_text(hjust = 0.5, size=10))

library(cowplot)
pdf(file.path(FIGURES_PATH, "application_figure_supp.pdf"), width=9.6, height=3.8)
wrap_plots(p1, p2, p3, p7, p4, p5, p6, p8, nrow = 2) & theme(legend.position = "")
# plot_grid(p1 + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .05), "cm")), 
#   p2 + theme(legend.position="") + theme(plot.margin = unit(c(.05, .0, .0, .0), "cm")), 
#     p3 + theme(legend.position="")+ theme(plot.margin = unit(c(.05, .05, .0, .0), "cm")), 
#     p4 + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .05), "cm")), 
#     p5 + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .0, .05, .0), "cm")),
#     p6 + theme(legend.position="")+ theme(plot.margin = unit(c(.0, .05, .05, .0), "cm")), nrow=2, align='vh')
dev.off()

pdf(file.path(FIGURES_PATH, "application_figure_supp_legend.pdf"), width=9.6, height=0.6)
ggdraw(get_legend(p1 + theme(legend.position = "bottom")))
dev.off()
