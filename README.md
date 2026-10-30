# Replicated diagnostic tests without a gold standard: closed-form solution

Code and results for

> Walter SD, Kostoulas PP. *Evaluation of replicated diagnostic and screening tests without a gold standard: an exact solution and its numerical evaluation.* (submitted to *Biometrical Journal*)

The repository reproduces the worked examples (Table 3) and the Monte Carlo simulation that compares the closed-form maximum likelihood solution with direct numerical maximization by multistart quasi-Newton BFGS (Table 1, Figures 1–3).

## Contents

| Path | What it is |
|---|---|
| `R/simulation_closed_form_vs_bfgs.R` | Full simulation: 320 scenarios × 1,000 datasets × 2 methods. Contains the closed-form estimator (`stephen_mle_3rep`), its Fisher information (`info_SW`) and the multistart BFGS fit (`mle_3rep_grad_logit_multistart`). |
| `R/worked_examples.R` | Closed-form estimates and 95% CIs for the two datasets (Table 3). |
| `R/load_functions.R` | Loads the estimator functions without starting the simulation. |
| `data/example_datasets.csv` | Counts of individuals with 0, 1, 2 or 3 positive results for the two examples (from the published sources cited in the paper). |
| `results/tables/` | Scenario-level and overall summaries from the final run used in the paper. |
| `results/figures/` | Bias, MSE, coverage and inadmissibility figures from the final run. |
| `results/FINAL_QC_SUMMARY.txt` | Quality-control summary of the final run. |

The replicate-level results (640,000 rows) are not stored here because of their size (31 MB as `.rds`, 238 MB as `.csv`). They are regenerated exactly by `R/simulation_closed_form_vs_bfgs.R` (see Reproducibility) and are available from the authors on request.

## Requirements

R (≥ 4.2) with `dplyr`, `tidyr`, `ggplot2`, `openxlsx` and `scales`; `parallel` ships with R. The worked examples need base R only.

## How to run

From the repository root:

```bash
# Table 3 (a few seconds)
Rscript R/worked_examples.R

# Full simulation (all 320 scenarios; uses all but two cores; restartable)
BFGS_RERUN_OUTPUT_DIR=results/full_run Rscript R/simulation_closed_form_vs_bfgs.R

# Quick check on the first 4 scenarios
BFGS_RERUN_SCENARIO_LIMIT=4 BFGS_RERUN_OUTPUT_DIR=results/test_run Rscript R/simulation_closed_form_vs_bfgs.R
```

Other settings (`BFGS_RERUN_NSIM`, `BFGS_RERUN_BASE_SEED`, `BFGS_RERUN_CORES`) are documented at the top of the simulation script.

## Simulation design

Counts of individuals with 0–3 positive results were generated from the multinomial distribution for all combinations of p ∈ {0.05, 0.10, 0.15, 0.20}, Se = 1 − α ∈ {0.60, 0.70, 0.80, 0.90}, Sp = γ ∈ {0.90, 0.93, 0.96, 0.99} and N ∈ {100, 250, 500, 750, 1000}, with 1,000 datasets per scenario. Scenario *i* uses the seed `20260723 + i`, so every scenario can be regenerated on its own. Both methods are applied to the same simulated counts.

## Reproducibility

The simulated counts are reproduced exactly from the seeds, and so are the closed-form estimates and the admissibility status of both methods. In a few sparse-data datasets, where the likelihood is nearly flat, numerical BFGS can stop at slightly different points depending on the platform and R version.

## Licence

MIT (see `LICENSE`).
