## Worked examples (Table 3 of the paper): closed-form estimates of
## prevalence (p), sensitivity (Se) and specificity (Sp), with 95% CIs.
## Run from the repository root:  Rscript R/worked_examples.R

source("R/load_functions.R")
load_estimator_functions("R/simulation_closed_form_vs_bfgs.R")

datasets <- read.csv("data/example_datasets.csv")

fmt <- function(est, ci) sprintf("%.2f (%.2f-%.2f)", est, ci[1], ci[2])

table3 <- do.call(rbind, lapply(seq_len(nrow(datasets)), function(i) {
  y <- unlist(datasets[i, c("y0", "y1", "y2", "y3")])
  fit <- stephen_mle_3rep(y, ci_type = "theta")
  data.frame(
    Dataset = datasets$dataset[i],
    p  = fmt(fit$p_hat,  fit$p_CI_95),
    Se = fmt(fit$Se_hat, fit$Se_CI_95),
    Sp = fmt(fit$Sp_hat, fit$Sp_CI_95)
  )
}))

print(table3, row.names = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(table3, "results/tables/table3_worked_examples.csv", row.names = FALSE)
