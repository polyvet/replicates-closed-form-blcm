## ============================================================
## CORRECTED FULL SIMULATION: CLOSED FORM VS MULTISTART BFGS
## ============================================================
## This script is a corrected, self-contained replacement for
## all_tries_Walter_final_REVISED_raw_results.R.
##
## Corrections:
##   1. Drops inherited R names during the mirror transformation.
##   2. Uses deterministic multistart BFGS, including the closed-form
##      estimate as a start whenever it is admissible.
##   3. Selects the admissible fit with the largest log-likelihood.
##   4. Separates optimization, point-estimation, Hessian/variance,
##      and confidence-interval failure stages.
##   5. Saves restartable scenario checkpoints.
##   6. Presents all four prevalence values together in the figures,
##      matching the clearer layout of the original all-tries analysis.
##   7. Rejects nonpositive moment variance a2-a1^2 rather than
##      replacing a negative value by its absolute value.
##   8. Uses a floating-point tolerance to identify estimates that
##      are mathematically on the boundary of the unit interval.
##   9. Can rapidly regenerate final outputs from a completed RDS
##      without rerunning the simulated datasets or BFGS fits.
##
## Default run: 320 scenarios x 1,000 replicates x 2 methods.
## All outputs are written beside this script in:
##   simulation_final_outputs_corrected/
##
## Optional environment variables:
##   BFGS_RERUN_NSIM             (default 1000)
##   BFGS_RERUN_BASE_SEED        (default 20260723)
##   BFGS_RERUN_CORES            (default detected cores minus 2)
##   BFGS_RERUN_SCENARIO_LIMIT   (default 0 = all scenarios)
##   BFGS_RERUN_OUTPUT_DIR       (default folder beside this script)
##   BFGS_REPROCESS_INPUT_RDS    (optional completed replicate RDS)
## ============================================================

script_file_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_dir <- if (length(script_file_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_file_argument)))
} else {
  normalizePath(getwd())
}

simulate_y_3rep <- function(N, p, Se, Sp) {
  # Diseased: Y|D=1 ~ Binomial(3, Se)
  a <- c((1 - Se)^3,
         3 * Se * (1 - Se)^2,
         3 * Se^2 * (1 - Se),
         Se^3)
  # Non-diseased: Y|D=0 ~ Binomial(3, 1-Sp)
  b <- c(Sp^3,
         3 * (1 - Sp) * Sp^2,
         3 * (1 - Sp)^2 * Sp,
         (1 - Sp)^3)
  
  pi_k <- p * a + (1 - p) * b      # mixture probabilities for k=0..3
  
  # Multinomial counts y0..y3
  y <- as.vector(rmultinom(1, size = N, prob = pi_k))
  names(y) <- paste0("y", 0:3)
  y
}

invlogit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))

stephen_ci_logit_from_Itheta <- function(theta_hat, I_theta) {
  # theta_hat: named c(p=, Se=, Sp=) strictly inside (0,1)
  # I_theta : 3x3 Fisher information for theta=(p,Se,Sp)
  
  stopifnot(all(c("p","Se","Sp") %in% names(theta_hat)))
  if (any(theta_hat <= 0) || any(theta_hat >= 1)) {
    return(matrix(NA_real_, 3, 2, dimnames=list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))))
  }
  if (!all(is.finite(I_theta))) {
    return(matrix(NA_real_, 3, 2, dimnames=list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))))
  }
  
  # Jacobian dtheta/deta
  d <- theta_hat * (1 - theta_hat)   # elementwise
  J <- diag(as.numeric(d))
  
  # Transform Fisher info: I_eta = J^T I_theta J = J I_theta J
  I_eta <- J %*% I_theta %*% J
  
  # Var(eta) = I_eta^{-1}
  V_eta <- tryCatch(solve(I_eta), error = function(e) matrix(NA_real_, 3, 3))
  if (!all(is.finite(V_eta))) {
    return(matrix(NA_real_, 3, 2, dimnames=list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))))
  }
  
  se_eta <- sqrt(diag(V_eta))
  eta_hat <- logit(theta_hat)
  
  z <- qnorm(0.975)
  ci_eta <- cbind(eta_hat - z * se_eta, eta_hat + z * se_eta)
  
  # Back-transform endpoints
  ci_theta <- apply(ci_eta, 2, invlogit)
  rownames(ci_theta) <- c("p","Se","Sp")
  colnames(ci_theta) <- c("CI_lower_95","CI_upper_95")
  ci_theta
}

compute_a <- function(Se) {
  c((1 - Se)^3,
    3 * Se * (1 - Se)^2,
    3 * Se^2 * (1 - Se),
    Se^3)
}
compute_b <- function(Sp) {
  c(Sp^3,
    3 * (1 - Sp) * Sp^2,
    3 * (1 - Sp)^2 * Sp,
    (1 - Sp)^3)
}

da_dSe <- function(Se) {
  c(-3 * (1 - Se)^2,
    3 * (1 - Se) * (1 - 3 * Se),
    3 * Se * (2 - 3 * Se),
    3 * Se^2)
}
d2a_dSe2 <- function(Se) {
  c(6 * (1 - Se),
    18 * Se - 12,
    6 - 18 * Se,
    6 * Se)
}

db_dSp <- function(Sp) {
  c(3 * Sp^2,
    3 * Sp * (2 - 3 * Sp),
    3 * (1 - Sp) * (1 - 3 * Sp),
    -3 * (1 - Sp)^2)
}
d2b_dSp2 <- function(Sp) {
  c(6 * Sp,
    6 - 18 * Sp,
    18 * Sp - 12,
    6 * (1 - Sp))
}

# Returns loglik, gradient, and Hessian w.r.t eta=(logit p, logit Se, logit Sp)
loglik_grad_hess_logit <- function(eta, y) {
  eta_p  <- eta[1]; eta_Se <- eta[2]; eta_Sp <- eta[3]
  
  p  <- invlogit(eta_p)
  Se <- invlogit(eta_Se)
  Sp <- invlogit(eta_Sp)
  
  a <- compute_a(Se)
  b <- compute_b(Sp)
  pi <- p * a + (1 - p) * b
  if (any(!is.finite(pi)) || any(pi <= 0)) {
    return(list(loglik = -Inf,
                grad   = rep(NA_real_, 3),
                hess   = matrix(NA_real_, 3, 3)))
  }
  
  # loglik
  loglik <- sum(y * log(pi))
  
  # first derivatives in theta-space
  dpi_dp  <- a - b
  dpi_dSe <- p * da_dSe(Se)
  dpi_dSp <- (1 - p) * db_dSp(Sp)
  
  # second derivatives in theta-space
  d2pi_dp2   <- rep(0, 4)
  d2pi_dpSe  <- da_dSe(Se)
  d2pi_dpSp  <- -db_dSp(Sp)
  d2pi_dSe2  <- p * d2a_dSe2(Se)
  d2pi_dSp2  <- (1 - p) * d2b_dSp2(Sp)
  d2pi_dSeSp <- rep(0, 4)
  
  # theta-gradient of loglik: g_theta = sum_k y_k/pi_k * dpi/dtheta
  w <- y / pi
  g_theta <- c(
    sum(w * dpi_dp),
    sum(w * dpi_dSe),
    sum(w * dpi_dSp)
  )
  
  # theta-Hessian of loglik:
  # H_theta = sum_k [ y_k/pi_k * d2pi_k  -  y_k/pi_k^2 * (dpi_k dpi_k^T) ]
  H_theta <- matrix(0, 3, 3)
  for (k in 1:4) {
    gk <- c(dpi_dp[k], dpi_dSe[k], dpi_dSp[k])
    Hk <- matrix(
      c(d2pi_dp2[k],  d2pi_dpSe[k],  d2pi_dpSp[k],
        d2pi_dpSe[k], d2pi_dSe2[k],  d2pi_dSeSp[k],
        d2pi_dpSp[k], d2pi_dSeSp[k], d2pi_dSp2[k]),
      3, 3, byrow = TRUE
    )
    H_theta <- H_theta + (y[k] / pi[k]) * Hk - (y[k] / (pi[k]^2)) * (gk %o% gk)
  }
  
  # chain rule: eta -> theta
  # theta_i = invlogit(eta_i), so dtheta_i/deta_i = theta_i(1-theta_i)
  d <- c(p * (1 - p), Se * (1 - Se), Sp * (1 - Sp))
  D <- diag(d)
  
  # second derivatives of theta_i wrt eta_i: d2theta/deta2 = d*(1-2*theta)
  d2 <- c(d[1] * (1 - 2 * p),
          d[2] * (1 - 2 * Se),
          d[3] * (1 - 2 * Sp))
  
  # gradient in eta-space
  g_eta <- as.vector(D %*% g_theta)
  
  # Hessian in eta-space: H_eta = D H_theta D + diag( g_theta * d2 )
  H_eta <- D %*% H_theta %*% D
  diag(H_eta) <- diag(H_eta) + g_theta * d2
  
  list(loglik = loglik, grad = g_eta, hess = H_eta)
}



mle_3rep_grad_logit <- function(y,
                                init = c(p = 0.2, Se = 0.7, Sp = 0.9)) {
  eta_init <- logit(init)
  
  negloglik_eta <- function(eta) {
    out <- loglik_grad_hess_logit(eta, y)
    -out$loglik
  }
  
  neggrad_eta <- function(eta) {
    out <- loglik_grad_hess_logit(eta, y)
    -out$grad
  }
  
  fit <- optim(
    par    = eta_init,
    fn     = negloglik_eta,
    gr     = neggrad_eta,
    method = "BFGS",
    control = list(fnscale = 1, maxit = 1000)
  )
  
  eta_hat <- fit$par
  
  theta_core <- c(
    p  = unname(invlogit(eta_hat[1])),
    Se = unname(invlogit(eta_hat[2])),
    Sp = unname(invlogit(eta_hat[3]))
  )
  
  if (fit$convergence != 0) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = matrix(NA_real_, 3, 3, dimnames = list(c("p","Se","Sp"), c("p","Se","Sp"))),
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "optimizer_nonconvergence",
      switched = FALSE
    ))
  }
  
  out_hat <- loglik_grad_hess_logit(eta_hat, y)
  H_eta <- -out_hat$hess
  
  if (any(!is.finite(H_eta))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = matrix(NA_real_, 3, 3, dimnames = list(c("p","Se","Sp"), c("p","Se","Sp"))),
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "hessian_nonfinite",
      switched = FALSE
    ))
  }
  
  eig <- tryCatch(
    eigen(H_eta, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 3)
  )
  
  if (any(!is.finite(eig))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = matrix(NA_real_, 3, 3, dimnames = list(c("p","Se","Sp"), c("p","Se","Sp"))),
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "hessian_eigenvalues_nonfinite",
      switched = FALSE
    ))
  }
  
  if (min(eig) <= 1e-10) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = matrix(NA_real_, 3, 3, dimnames = list(c("p","Se","Sp"), c("p","Se","Sp"))),
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "hessian_not_positive_definite",
      switched = FALSE
    ))
  }
  
  vcov_eta <- tryCatch(solve(H_eta), error = function(e) NULL)
  if (is.null(vcov_eta) || any(!is.finite(vcov_eta))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = matrix(NA_real_, 3, 3, dimnames = list(c("p","Se","Sp"), c("p","Se","Sp"))),
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "vcov_eta_inversion_failed",
      switched = FALSE
    ))
  }
  
  p  <- theta_core["p"]
  Se <- theta_core["Se"]
  Sp <- theta_core["Sp"]
  
  J <- diag(c(p * (1 - p), Se * (1 - Se), Sp * (1 - Sp)))
  vcov_theta <- J %*% vcov_eta %*% t(J)
  dimnames(vcov_theta) <- list(c("p","Se","Sp"), c("p","Se","Sp"))
  
  if (any(!is.finite(vcov_theta))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = setNames(c(NA_real_, NA_real_, NA_real_), c("p","Se","Sp")),
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = vcov_theta,
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "vcov_theta_nonfinite",
      switched = FALSE
    ))
  }
  
  se <- sqrt(diag(vcov_theta))
  names(se) <- c("p","Se","Sp")
  
  switched <- FALSE
  if (is.finite(Se) && is.finite(Sp) && (Se + Sp < 1)) {
    switched <- TRUE
    theta_core <- c(
      p = unname(1 - theta_core["p"]),
      Se = unname(1 - theta_core["Sp"]),
      Sp = unname(1 - theta_core["Se"])
    )
    
    A <- matrix(
      c(-1,  0,  0,
        0,   0, -1,
        0,  -1,  0),
      nrow = 3, byrow = TRUE
    )
    
    vcov_theta <- A %*% vcov_theta %*% t(A)
    dimnames(vcov_theta) <- list(c("p","Se","Sp"), c("p","Se","Sp"))
    se <- sqrt(diag(vcov_theta))
    names(se) <- c("p","Se","Sp")
  }
  
  if (any(!is.finite(theta_core))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = se,
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = vcov_theta,
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "point_estimate_nonfinite_after_orientation",
      switched = switched
    ))
  }
  
  if (any(theta_core <= 0 | theta_core >= 1)) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = se,
      ci_95 = matrix(NA_real_, 3, 2, dimnames = list(c("p","Se","Sp"), c("CI_lower_95","CI_upper_95"))),
      loglik = -fit$value,
      vcov = vcov_theta,
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "point_estimate_outside_open_unit_interval",
      switched = switched
    ))
  }
  
  z <- qnorm(0.975)
  ci_lower <- pmax(0, theta_core - z * se)
  ci_upper <- pmin(1, theta_core + z * se)
  
  ci <- cbind(ci_lower, ci_upper)
  rownames(ci) <- c("p","Se","Sp")
  colnames(ci) <- c("CI_lower_95", "CI_upper_95")
  
  if (any(is.na(ci))) {
    return(list(
      par_hat = c("p.p" = unname(theta_core["p"]), "Se.Se" = unname(theta_core["Se"]), "Sp.Sp" = unname(theta_core["Sp"])),
      se = se,
      ci_95 = ci,
      loglik = -fit$value,
      vcov = vcov_theta,
      convergence = fit$convergence,
      inadmissible = TRUE,
      failure_reason = "ci_na",
      switched = switched
    ))
  }
  
  par_hat_out <- c("p.p" = unname(theta_core["p"]),
                   "Se.Se" = unname(theta_core["Se"]),
                   "Sp.Sp" = unname(theta_core["Sp"]))
  
  list(
    par_hat = par_hat_out,
    se      = se,
    ci_95   = ci,
    loglik  = -fit$value,
    vcov    = vcov_theta,
    convergence = fit$convergence,
    inadmissible = FALSE,
    failure_reason = "ok",
    switched = switched
  )
}

bfgs_failure_stage <- function(failure_reason) {
  if (is.null(failure_reason) || length(failure_reason) == 0L) {
    return("unknown")
  }
  reason <- as.character(failure_reason[1])
  if (identical(reason, "ok")) return("none")
  if (grepl("^error:", reason)) return("error")
  if (grepl("optim|convergence", reason)) return("optimization")
  if (grepl("^point_estimate", reason)) return("point_estimate")
  if (grepl("^hessian", reason)) return("hessian")
  if (grepl("^vcov", reason)) return("variance")
  if (grepl("^ci_|^confidence_interval", reason)) {
    return("confidence_interval")
  }
  "other"
}

bfgs_fit_is_admissible <- function(fit) {
  if (is.null(fit) || isTRUE(fit$inadmissible)) return(FALSE)
  par_hat <- fit$par_hat
  ci <- fit$ci_95
  all(
    is.finite(c(
      par_hat["p.p"], par_hat["Se.Se"], par_hat["Sp.Sp"],
      ci["p", ], ci["Se", ], ci["Sp", ]
    ))
  ) &&
    all(
      c(par_hat["p.p"], par_hat["Se.Se"], par_hat["Sp.Sp"]) > 0 &
        c(par_hat["p.p"], par_hat["Se.Se"], par_hat["Sp.Sp"]) < 1
    )
}

clip_bfgs_start <- function(x, epsilon = 1e-6) {
  x <- pmin(pmax(as.numeric(x), epsilon), 1 - epsilon)
  setNames(x, c("p", "Se", "Sp"))
}

mle_3rep_grad_logit_multistart <- function(
    y,
    closed_form_init = NULL) {

  deterministic_starts <- list(
    original_fixed = c(p = 0.50, Se = 0.70, Sp = 0.90),
    low_prevalence = c(p = 0.05, Se = 0.70, Sp = 0.95),
    moderate = c(p = 0.20, Se = 0.80, Sp = 0.95),
    high_accuracy = c(p = 0.10, Se = 0.90, Sp = 0.99),
    lower_accuracy = c(p = 0.20, Se = 0.60, Sp = 0.90)
  )

  closed_form_available <- !is.null(closed_form_init) &&
    length(closed_form_init) == 3L &&
    all(is.finite(closed_form_init)) &&
    all(closed_form_init > 0 & closed_form_init < 1)

  starts <- if (closed_form_available) {
    list(
      closed_form = clip_bfgs_start(closed_form_init),
      original_fixed = deterministic_starts$original_fixed
    )
  } else {
    deterministic_starts
  }

  run_start <- function(start_name, init) {
    error_message <- NA_character_
    fit <- tryCatch(
      suppressWarnings(mle_3rep_grad_logit(y, init = clip_bfgs_start(init))),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(fit)) {
      fit <- list(
        par_hat = setNames(rep(NA_real_, 3), c("p.p", "Se.Se", "Sp.Sp")),
        se = setNames(rep(NA_real_, 3), c("p", "Se", "Sp")),
        ci_95 = matrix(
          NA_real_, 3, 2,
          dimnames = list(
            c("p", "Se", "Sp"),
            c("CI_lower_95", "CI_upper_95")
          )
        ),
        loglik = NA_real_,
        vcov = matrix(NA_real_, 3, 3),
        convergence = NA_integer_,
        inadmissible = TRUE,
        failure_reason = paste0("error: ", error_message),
        switched = NA
      )
    }
    fit$start_name <- start_name
    fit
  }

  fits <- Map(run_start, names(starts), starts)

  # If the two economical starts available with an admissible closed-form
  # estimate do not yield valid inference, expand to the complete
  # deterministic set. This keeps the full simulation practical while
  # retaining a robust rescue path.
  if (
    closed_form_available &&
      !any(vapply(fits, bfgs_fit_is_admissible, logical(1)))
  ) {
    additional_names <- setdiff(
      names(deterministic_starts),
      names(starts)
    )
    additional_fits <- Map(
      run_start,
      additional_names,
      deterministic_starts[additional_names]
    )
    fits <- c(fits, additional_fits)
  }

  admissible <- vapply(fits, bfgs_fit_is_admissible, logical(1))
  loglik_values <- vapply(
    fits,
    function(fit) {
      value <- safe_num(fit$loglik)
      if (is.finite(value)) value else -Inf
    },
    numeric(1)
  )

  candidate_indexes <- if (any(admissible)) {
    which(admissible)
  } else {
    which(is.finite(loglik_values))
  }
  if (length(candidate_indexes) == 0L) {
    candidate_indexes <- seq_along(fits)
  }
  best_index <- candidate_indexes[
    which.max(loglik_values[candidate_indexes])
  ]
  best <- fits[[best_index]]

  point_values <- c(
    best$par_hat["p.p"],
    best$par_hat["Se.Se"],
    best$par_hat["Sp.Sp"]
  )
  best$selected_start <- best$start_name
  best$n_starts_attempted <- length(fits)
  best$closed_form_start_available <- closed_form_available
  best$point_estimation_success <- isTRUE(best$convergence == 0L) &&
    all(is.finite(point_values)) &&
    all(point_values > 0 & point_values < 1)
  best$inference_success <- bfgs_fit_is_admissible(best)
  best$failure_stage <- bfgs_failure_stage(best$failure_reason)
  best
}

info_SW <- function(p, Se, Sp, N) {
  # Expected Fisher information for theta = (p, Se, Sp)
  # under the 3-replicate independence model, using
  # the full expression:
  # I(θ) = N * sum_k [ (1 / π_k) * (∂π_k/∂θ)(∂π_k/∂θ)^T  -  ∂²π_k/∂θ∂θ^T ]
  
  # Basic validity checks
  if (p <= 0 || p >= 1 || Se <= 0 || Se >= 1 || Sp <= 0 || Sp >= 1) {
    return(matrix(NA_real_, nrow = 3, ncol = 3))
  }
  
  ## --- a_k and b_k ---
  a <- c(
    (1 - Se)^3,
    3 * Se * (1 - Se)^2,
    3 * Se^2 * (1 - Se),
    Se^3
  )
  
  b <- c(
    Sp^3,
    3 * (1 - Sp) * Sp^2,
    3 * (1 - Sp)^2 * Sp,
    (1 - Sp)^3
  )
  
  ## --- mixture probabilities π_k = p a_k + (1-p) b_k ---
  pi <- p * a + (1 - p) * b
  if (any(pi <= 0)) {
    return(matrix(NA_real_, nrow = 3, ncol = 3))
  }
  
  ## --- first derivatives of a_k wrt Se ---
  da_dSe <- c(
    -3 * (1 - Se)^2,
    3 * (1 - Se) * (1 - 3 * Se),
    3 * Se * (2 - 3 * Se),
    3 * Se^2
  )
  
  ## --- second derivatives of a_k wrt Se ---
  # a0'' = 6(1 - Se)
  # a1'' = 18 Se - 12
  # a2'' = 6 - 18 Se
  # a3'' = 6 Se
  d2a_dSe2 <- c(
    6 * (1 - Se),
    18 * Se - 12,
    6 - 18 * Se,
    6 * Se
  )
  
  ## --- first derivatives of b_k wrt Sp ---
  db_dSp <- c(
    3 * Sp^2,
    3 * Sp * (2 - 3 * Sp),
    3 * (1 - Sp) * (1 - 3 * Sp),
    -3 * (1 - Sp)^2
  )
  
  ## --- second derivatives of b_k wrt Sp ---
  # b0'' = 6 Sp
  # b1'' = 6 - 18 Sp
  # b2'' = 18 Sp - 12
  # b3'' = 6(1 - Sp)
  d2b_dSp2 <- c(
    6 * Sp,
    6 - 18 * Sp,
    18 * Sp - 12,
    6 * (1 - Sp)
  )
  
  ## --- first derivatives of π_k wrt (p, Se, Sp) ---
  dpi_dp  <- a - b              # ∂π_k/∂p
  dpi_dSe <- p * da_dSe         # ∂π_k/∂Se
  dpi_dSp <- (1 - p) * db_dSp   # ∂π_k/∂Sp
  
  ## --- second derivatives of π_k wrt (p, Se, Sp) ---
  # d²π_k/dp² = 0 (π_k linear in p)
  d2pi_dp2   <- rep(0, 4)
  # d²π_k/dp dSe = d/dSe(a_k - b_k) = da_k/dSe
  d2pi_dpSe  <- da_dSe
  # d²π_k/dp dSp = d/dSp( (1-p) db_k/dSp ) = -(db_k/dSp)
  d2pi_dpSp  <- -db_dSp
  # d²π_k/dSe² = p * d²a_k/dSe²
  d2pi_dSe2  <- p * d2a_dSe2
  # d²π_k/dSp² = (1-p) * d²b_k/dSp²
  d2pi_dSp2  <- (1 - p) * d2b_dSp2
  # d²π_k/dSe dSp = 0 (a depends only on Se, b only on Sp)
  d2pi_dSeSp <- rep(0, 4)
  
  ## --- assemble Fisher information ---
  I <- matrix(0, nrow = 3, ncol = 3)
  # parameter order: θ = (p, Se, Sp)
  
  for (k in 1:4) {
    # gradient vector g_k = (∂π_k/∂p, ∂π_k/∂Se, ∂π_k/∂Sp)
    gk <- c(dpi_dp[k], dpi_dSe[k], dpi_dSp[k])
    
    # Hessian of π_k wrt (p, Se, Sp)
    Hk <- matrix(
      c(
        d2pi_dp2[k],   d2pi_dpSe[k],   d2pi_dpSp[k],
        d2pi_dpSe[k],  d2pi_dSe2[k],   d2pi_dSeSp[k],
        d2pi_dpSp[k],  d2pi_dSeSp[k],  d2pi_dSp2[k]
      ),
      nrow = 3, byrow = TRUE
    )
    
    I <- I + (1 / pi[k]) * (gk %o% gk) - Hk
  }
  
  I * N
}

stephen_mle_3rep <- function(y, ci_type = c("theta","logit")) {
  
  logit  <- function(x) qlogis(x)
  expit  <- function(x) plogis(x)
  clamp01 <- function(x, eps = 1e-8) pmin(pmax(x, eps), 1 - eps)
  
  ci_type <- match.arg(ci_type)
  
  if (length(y) != 4L) stop("y must be a length-4 vector: c(y0,y1,y2,y3)")
  N <- sum(y)
  if (N <= 0) stop("Total N must be > 0")
  
  ## 1. Lambdas
  lambda <- y / N
  lam0 <- lambda[1]
  lam1 <- lambda[2]
  lam2 <- lambda[3]
  lam3 <- lambda[4]
  
  ## 2. a0, a1, a2
  a0 <- lam0
  a1 <- (3 * lam0 + 2 * lam1 + lam2) / 3
  a2 <- (3 * lam0 + lam1) / 3
  
  ## 3. k
  k_num    <- a0 + 2 * a1^3 - 3 * a1 * a2
  k_den0   <- a2 - a1^2
  numerical_tolerance <- sqrt(.Machine$double.eps)
  
  if (!is.finite(k_den0) || k_den0 <= numerical_tolerance) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "k_den_nonpositive_or_nonfinite",
      p_hat = NA, Se_hat = NA, Sp_hat = NA,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      orientation = NA
    ))
  }

  k_den <- k_den0^(1.5)
  
  k <- k_num / k_den
  
  ## 4. X and prevalence p
  X <- (sqrt(k^2 + 4) - k) / 2
  
  if (!is.finite(X) || X <= 0) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "X_nonfinite_or_nonpositive",
      p_hat = NA, Se_hat = NA, Sp_hat = NA,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      orientation = NA
    ))
  }
  
  p_hat <- (X^2) / (1 + X^2)
  
  ## 5. Gamma and alpha roots
  s <- sqrt(k_den0)
  
  gamma_pos <- a1 + X * s
  gamma_neg <- a1 - X * s
  alpha_pos <- a1 + s / X
  alpha_neg <- a1 - s / X
  
  ## 6. Choose root by fit to a0,a1,a2
  a_vec <- c(a0 = a0, a1 = a1, a2 = a2)
  
  fitted_a <- function(alpha, gamma, p) {
    a0_hat <- p * alpha^3 + (1 - p) * gamma^3
    a1_hat <- p * alpha    + (1 - p) * gamma
    a2_hat <- p * alpha^2  + (1 - p) * gamma^2
    c(a0_hat = a0_hat, a1_hat = a1_hat, a2_hat = a2_hat)
  }
  
  combos <- list(
    pos_gamma_pos_alpha = list(gamma = gamma_pos, alpha = alpha_pos),
    pos_gamma_neg_alpha = list(gamma = gamma_pos, alpha = alpha_neg),
    neg_gamma_pos_alpha = list(gamma = gamma_neg, alpha = alpha_pos),
    neg_gamma_neg_alpha = list(gamma = gamma_neg, alpha = alpha_neg)
  )
  
  errs <- sapply(names(combos), function(nm) {
    g <- combos[[nm]]$gamma
    a <- combos[[nm]]$alpha
    if (!all(is.finite(c(g, a, p_hat)))) return(Inf)
    a_hat <- fitted_a(alpha = a, gamma = g, p = p_hat)
    sum((a_hat - a_vec)^2)
  })
  
  if (all(!is.finite(errs))) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "all_root_combinations_invalid",
      p_hat = NA, Se_hat = NA, Sp_hat = NA,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      orientation = NA
    ))
  }
  
  best_name  <- names(combos)[which.min(errs)]
  best_gamma <- combos[[best_name]]$gamma
  best_alpha <- combos[[best_name]]$alpha
  
  ## 7. Map to Se and Sp, apply identifiability orientation
  Se_raw <- 1 - best_alpha
  Sp_raw <- best_gamma
  p_raw  <- p_hat
  
  p_alt  <- 1 - p_raw
  Se_alt <- 1 - Sp_raw
  Sp_alt <- 1 - Se_raw
  
  if (is.finite(Se_raw) && is.finite(Sp_raw) && (Se_raw + Sp_raw >= 1)) {
    p_hat_final  <- p_raw
    Se_hat_final <- Se_raw
    Sp_hat_final <- Sp_raw
    orientation  <- "raw"
  } else {
    p_hat_final  <- p_alt
    Se_hat_final <- Se_alt
    Sp_hat_final <- Sp_alt
    orientation  <- "switched"
  }
  
  p_hat  <- p_hat_final
  Se_hat <- Se_hat_final
  Sp_hat <- Sp_hat_final
  
  if (any(!is.finite(c(p_hat, Se_hat, Sp_hat)))) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "point_estimate_nonfinite_after_orientation",
      p_hat = NA, Se_hat = NA, Sp_hat = NA,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      orientation = orientation
    ))
  }
  
  if (p_hat <= numerical_tolerance ||
      p_hat >= 1 - numerical_tolerance ||
      Se_hat <= numerical_tolerance ||
      Se_hat >= 1 - numerical_tolerance ||
      Sp_hat <= numerical_tolerance ||
      Sp_hat >= 1 - numerical_tolerance) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "point_estimate_outside_open_unit_interval",
      p_hat = NA, Se_hat = NA, Sp_hat = NA,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      orientation = orientation
    ))
  }
  
  p_CI <- Se_CI <- Sp_CI <- rep(NA_real_, 2)
  se   <- rep(NA_real_, 3)
  
  I_mat <- tryCatch(
    info_SW(p = p_hat, Se = Se_hat, Sp = Sp_hat, N = N),
    error = function(e) matrix(NA_real_, 3, 3)
  )
  
  if (!all(is.finite(I_mat))) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "info_matrix_nonfinite",
      p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      info_mat = I_mat,
      orientation = orientation
    ))
  }
  
  Sigma_hat <- tryCatch(solve(I_mat), error = function(e) matrix(NA_real_, 3, 3))
  
  if (!all(is.finite(Sigma_hat))) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "info_matrix_inversion_failed",
      p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
      p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
      info_mat = I_mat,
      orientation = orientation
    ))
  }
  
  theta_hat <- c(p = p_hat, Se = Se_hat, Sp = Sp_hat)
  se <- sqrt(diag(Sigma_hat))
  
  if (ci_type == "theta") {
    z <- qnorm(0.975)
    p_CI  <- pmin(1, pmax(0, c(p_hat  - z * se[1], p_hat  + z * se[1])))
    Se_CI <- pmin(1, pmax(0, c(Se_hat - z * se[2], Se_hat + z * se[2])))
    Sp_CI <- pmin(1, pmax(0, c(Sp_hat - z * se[3], Sp_hat + z * se[3])))
  } else {
    z <- qnorm(0.975)
    th <- clamp01(theta_hat)
    eta_hat <- qlogis(th)
    J <- diag(1 / (th * (1 - th)))
    Sigma_eta <- J %*% Sigma_hat %*% t(J)
    
    if (!all(is.finite(Sigma_eta))) {
      return(list(
        inadmissible = TRUE,
        failure_reason = "logit_variance_nonfinite",
        p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
        p_CI_95 = c(NA, NA), Se_CI_95 = c(NA, NA), Sp_CI_95 = c(NA, NA),
        info_mat = I_mat,
        orientation = orientation
      ))
    }
    
    v_eta <- pmax(diag(Sigma_eta), 0)
    se_eta <- sqrt(v_eta)
    
    eta_lo <- eta_hat - z * se_eta
    eta_hi <- eta_hat + z * se_eta
    
    CI <- rbind(eta_lo, eta_hi)
    CI <- apply(CI, 2, plogis)
    
    p_CI  <- CI[, 1]
    Se_CI <- CI[, 2]
    Sp_CI <- CI[, 3]
  }
  
  if (any(is.na(c(p_CI, Se_CI, Sp_CI)))) {
    return(list(
      inadmissible = TRUE,
      failure_reason = "ci_na",
      p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
      p_CI_95 = p_CI, Se_CI_95 = Se_CI, Sp_CI_95 = Sp_CI,
      info_mat = I_mat,
      orientation = orientation
    ))
  }
  
  list(
    inadmissible = FALSE,
    failure_reason = "ok",
    ci_type  = ci_type,
    N        = N,
    y        = y,
    lambda   = lambda,
    a        = a_vec,
    k        = k,
    X        = X,
    p_hat    = p_hat,
    Se_hat   = Se_hat,
    Sp_hat   = Sp_hat,
    p_CI_95  = p_CI,
    Se_CI_95 = Se_CI,
    Sp_CI_95 = Sp_CI,
    se_p     = se[1],
    se_Se    = se[2],
    se_Sp    = se[3],
    gamma_pos = gamma_pos,
    gamma_neg = gamma_neg,
    alpha_pos = alpha_pos,
    alpha_neg = alpha_neg,
    best_combo = best_name,
    best_gamma = best_gamma,
    best_alpha = best_alpha,
    fit_errors = errs,
    info_mat   = I_mat,
    orientation = orientation
  )
}


# ============================================================
# FINAL RE-RUN: SAVE REPLICATE-LEVEL RESULTS AND ALL SUMMARIES
# ============================================================
# This section replaces the previous run_scenario_both() and all
# code below it. It compares:
#   1. Closed-form analytic estimator with theta-scale Wald CIs
#   2. Numerical BFGS estimator
#
# It saves each simulated dataset and both method results, enabling:
#   - bias, empirical SD, MSE and coverage
#   - 2x2 admissibility cross-tabs
#   - failure-reason summaries
#   - matched head-to-head comparisons when both methods work
# ============================================================

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(scales)
})

safe_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(default)
  as.numeric(x[1])
}

safe_ci <- function(x) {
  if (is.null(x) || length(x) != 2L || any(!is.finite(x))) {
    return(c(NA_real_, NA_real_))
  }
  as.numeric(x)
}

extract_closed_form <- function(fit, error_message = NA_character_) {
  if (is.null(fit)) {
    return(list(
      admissible = FALSE,
      failure_reason = ifelse(
        is.na(error_message), "closed_form_error", paste0("error: ", error_message)
      ),
      p_hat = NA_real_, Se_hat = NA_real_, Sp_hat = NA_real_,
      p_lower = NA_real_, p_upper = NA_real_,
      Se_lower = NA_real_, Se_upper = NA_real_,
      Sp_lower = NA_real_, Sp_upper = NA_real_,
      orientation = NA_character_,
      optimizer_convergence = NA_integer_,
      point_estimation_success = FALSE,
      inference_success = FALSE,
      failure_stage = "error",
      selected_start = NA_character_,
      n_starts_attempted = NA_integer_,
      closed_form_start_available = FALSE,
      loglik = NA_real_
    ))
  }

  p_ci  <- safe_ci(fit$p_CI_95)
  Se_ci <- safe_ci(fit$Se_CI_95)
  Sp_ci <- safe_ci(fit$Sp_CI_95)

  p_hat  <- safe_num(fit$p_hat)
  Se_hat <- safe_num(fit$Se_hat)
  Sp_hat <- safe_num(fit$Sp_hat)

  admissible <- !isTRUE(fit$inadmissible) &&
    all(is.finite(c(p_hat, Se_hat, Sp_hat, p_ci, Se_ci, Sp_ci))) &&
    all(c(p_hat, Se_hat, Sp_hat) > 0 & c(p_hat, Se_hat, Sp_hat) < 1)

  reason <- if (admissible) {
    "ok"
  } else if (!is.null(fit$failure_reason)) {
    as.character(fit$failure_reason)
  } else if (any(!is.finite(c(p_hat, Se_hat, Sp_hat)))) {
    "point_estimate_nonfinite"
  } else if (any(c(p_hat, Se_hat, Sp_hat) <= 0 |
                 c(p_hat, Se_hat, Sp_hat) >= 1)) {
    "point_estimate_outside_open_unit_interval"
  } else {
    "ci_or_information_failure"
  }

  list(
    admissible = admissible,
    failure_reason = reason,
    p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
    p_lower = p_ci[1], p_upper = p_ci[2],
    Se_lower = Se_ci[1], Se_upper = Se_ci[2],
    Sp_lower = Sp_ci[1], Sp_upper = Sp_ci[2],
    orientation = ifelse(
      is.null(fit$orientation), NA_character_, as.character(fit$orientation)
    )
  )
}

extract_bfgs <- function(fit, error_message = NA_character_) {
  if (is.null(fit)) {
    return(list(
      admissible = FALSE,
      failure_reason = ifelse(
        is.na(error_message), "BFGS_error", paste0("error: ", error_message)
      ),
      p_hat = NA_real_, Se_hat = NA_real_, Sp_hat = NA_real_,
      p_lower = NA_real_, p_upper = NA_real_,
      Se_lower = NA_real_, Se_upper = NA_real_,
      Sp_lower = NA_real_, Sp_upper = NA_real_,
      orientation = NA_character_
    ))
  }

  par_hat <- fit$par_hat
  p_hat <- safe_num(par_hat["p.p"])
  Se_hat <- safe_num(par_hat["Se.Se"])
  Sp_hat <- safe_num(par_hat["Sp.Sp"])

  ci <- fit$ci_95
  if (is.null(ci) || !all(c("p", "Se", "Sp") %in% rownames(ci))) {
    p_ci <- Se_ci <- Sp_ci <- c(NA_real_, NA_real_)
  } else {
    p_ci  <- safe_ci(ci["p", ])
    Se_ci <- safe_ci(ci["Se", ])
    Sp_ci <- safe_ci(ci["Sp", ])
  }

  admissible <- !isTRUE(fit$inadmissible) &&
    all(is.finite(c(p_hat, Se_hat, Sp_hat, p_ci, Se_ci, Sp_ci))) &&
    all(c(p_hat, Se_hat, Sp_hat) > 0 & c(p_hat, Se_hat, Sp_hat) < 1)

  reason <- if (admissible) {
    "ok"
  } else if (!is.null(fit$failure_reason)) {
    as.character(fit$failure_reason)
  } else {
    "unknown_BFGS_failure"
  }

  list(
    admissible = admissible,
    failure_reason = reason,
    p_hat = p_hat, Se_hat = Se_hat, Sp_hat = Sp_hat,
    p_lower = p_ci[1], p_upper = p_ci[2],
    Se_lower = Se_ci[1], Se_upper = Se_ci[2],
    Sp_lower = Sp_ci[1], Sp_upper = Sp_ci[2],
    orientation = ifelse(
      isTRUE(fit$switched), "switched", "raw"
    ),
    optimizer_convergence = safe_num(fit$convergence),
    point_estimation_success = isTRUE(fit$point_estimation_success),
    inference_success = isTRUE(fit$inference_success),
    failure_stage = ifelse(
      is.null(fit$failure_stage),
      bfgs_failure_stage(reason),
      as.character(fit$failure_stage)
    ),
    selected_start = ifelse(
      is.null(fit$selected_start),
      NA_character_,
      as.character(fit$selected_start)
    ),
    n_starts_attempted = safe_num(fit$n_starts_attempted),
    closed_form_start_available = isTRUE(fit$closed_form_start_available),
    loglik = safe_num(fit$loglik)
  )
}

run_scenario_replicates <- function(
    scenario_id, p_true, Se_true, Sp_true, N, nsim = 1000,
    base_seed = 20260723L) {

  # Separate deterministic seed for each scenario.
  set.seed(as.integer(base_seed + scenario_id))

  rows <- vector("list", 2L * nsim)

  for (s in seq_len(nsim)) {
    y <- simulate_y_3rep(
      N = N, p = p_true, Se = Se_true, Sp = Sp_true
    )

    err_cf <- NA_character_
    fit_cf <- tryCatch(
      stephen_mle_3rep(y, ci_type = "theta"),
      error = function(e) {
        err_cf <<- conditionMessage(e)
        NULL
      }
    )
    cf <- extract_closed_form(fit_cf, err_cf)

    err_bfgs <- NA_character_
    closed_form_init <- if (isTRUE(cf$admissible)) {
      c(p = cf$p_hat, Se = cf$Se_hat, Sp = cf$Sp_hat)
    } else {
      NULL
    }
    fit_bfgs <- tryCatch(
      mle_3rep_grad_logit_multistart(
        y,
        closed_form_init = closed_form_init
      ),
      error = function(e) {
        err_bfgs <<- conditionMessage(e)
        NULL
      }
    )
    bg <- extract_bfgs(fit_bfgs, err_bfgs)

    common <- list(
      scenario_id = scenario_id,
      replicate = s,
      p_true = p_true,
      Se_true = Se_true,
      Sp_true = Sp_true,
      alpha_true = 1 - Se_true,
      gamma_true = Sp_true,
      N = N,
      y0 = unname(y[1]),
      y1 = unname(y[2]),
      y2 = unname(y[3]),
      y3 = unname(y[4])
    )

    rows[[2L * s - 1L]] <- data.frame(
      common,
      method = "Closed-form",
      admissible = cf$admissible,
      failure_reason = cf$failure_reason,
      orientation = cf$orientation,
      p_hat = cf$p_hat,
      Se_hat = cf$Se_hat,
      Sp_hat = cf$Sp_hat,
      alpha_hat = 1 - cf$Se_hat,
      gamma_hat = cf$Sp_hat,
      p_lower = cf$p_lower,
      p_upper = cf$p_upper,
      Se_lower = cf$Se_lower,
      Se_upper = cf$Se_upper,
      Sp_lower = cf$Sp_lower,
      Sp_upper = cf$Sp_upper,
      alpha_lower = 1 - cf$Se_upper,
      alpha_upper = 1 - cf$Se_lower,
      gamma_lower = cf$Sp_lower,
      gamma_upper = cf$Sp_upper,
      stringsAsFactors = FALSE
    )

    rows[[2L * s]] <- data.frame(
      common,
      method = "Numerical BFGS",
      admissible = bg$admissible,
      failure_reason = bg$failure_reason,
      orientation = bg$orientation,
      p_hat = bg$p_hat,
      Se_hat = bg$Se_hat,
      Sp_hat = bg$Sp_hat,
      alpha_hat = 1 - bg$Se_hat,
      gamma_hat = bg$Sp_hat,
      p_lower = bg$p_lower,
      p_upper = bg$p_upper,
      Se_lower = bg$Se_lower,
      Se_upper = bg$Se_upper,
      Sp_lower = bg$Sp_lower,
      Sp_upper = bg$Sp_upper,
      alpha_lower = 1 - bg$Se_upper,
      alpha_upper = 1 - bg$Se_lower,
      gamma_lower = bg$Sp_lower,
      gamma_upper = bg$Sp_upper,
      optimizer_convergence = bg$optimizer_convergence,
      point_estimation_success = bg$point_estimation_success,
      inference_success = bg$inference_success,
      failure_stage = bg$failure_stage,
      selected_start = bg$selected_start,
      n_starts_attempted = bg$n_starts_attempted,
      closed_form_start_available = bg$closed_form_start_available,
      loglik = bg$loglik,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(rows)
}

# -----------------------------
# Simulation design
# -----------------------------
Se_vals <- c(0.60, 0.70, 0.80, 0.90)
Sp_vals <- c(0.90, 0.93, 0.96, 0.99)
p_vals  <- c(0.05, 0.10, 0.15, 0.20)
N_vals  <- c(100, 250, 500, 750, 1000)

param_grid <- expand.grid(
  p_true = p_vals,
  Se_true = Se_vals,
  Sp_true = Sp_vals,
  N = N_vals,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) |>
  mutate(scenario_id = row_number())

read_integer_environment <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value)) default else value
}

nsim <- read_integer_environment("BFGS_RERUN_NSIM", 1000L)
base_seed <- read_integer_environment("BFGS_RERUN_BASE_SEED", 20260723L)
scenario_limit <- read_integer_environment("BFGS_RERUN_SCENARIO_LIMIT", 0L)
if (scenario_limit > 0L) {
  param_grid <- param_grid |>
    slice_head(n = min(scenario_limit, nrow(param_grid)))
}

default_output_dir <- file.path(
  script_dir,
  "simulation_final_outputs_corrected"
)
output_dir <- Sys.getenv("BFGS_RERUN_OUTPUT_DIR", default_output_dir)
output_dir <- normalizePath(
  output_dir,
  winslash = "/",
  mustWork = FALSE
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
checkpoint_dir <- file.path(output_dir, "scenario_checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Optional fast reprocessing
# -----------------------------
reprocess_existing_results <- function(raw_results) {
  required_columns <- c(
    "scenario_id", "replicate", "method", "admissible",
    "failure_reason", "p_hat", "Se_hat", "Sp_hat",
    "y0", "y1", "y2", "y3", "N"
  )
  missing_columns <- setdiff(required_columns, names(raw_results))
  if (length(missing_columns) > 0L) {
    stop(
      "The reprocessing RDS is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (nrow(raw_results) != 640000L) {
    stop("Expected 640,000 replicate-method rows; found ", nrow(raw_results), ".")
  }

  tolerance <- sqrt(.Machine$double.eps)
  closed_form <- raw_results$method == "Closed-form"
  a1 <- (3 * raw_results$y0 + 2 * raw_results$y1 + raw_results$y2) /
    (3 * raw_results$N)
  a2 <- (3 * raw_results$y0 + raw_results$y1) / (3 * raw_results$N)
  moment_variance <- a2 - a1^2

  invalid_moment <- closed_form &
    (!is.finite(moment_variance) | moment_variance <= tolerance)
  near_boundary <- closed_form & !invalid_moment & raw_results$admissible &
    (
      raw_results$p_hat <= tolerance |
      raw_results$p_hat >= 1 - tolerance |
      raw_results$Se_hat <= tolerance |
      raw_results$Se_hat >= 1 - tolerance |
      raw_results$Sp_hat <= tolerance |
      raw_results$Sp_hat >= 1 - tolerance
    )
  affected_dataset_keys <- unique(paste(
    raw_results$scenario_id[invalid_moment | near_boundary],
    raw_results$replicate[invalid_moment | near_boundary],
    sep = ":"
  ))

  raw_results$admissible[invalid_moment | near_boundary] <- FALSE
  raw_results$failure_reason[invalid_moment] <-
    "k_den_nonpositive_or_nonfinite"
  raw_results$failure_reason[near_boundary] <-
    "point_estimate_outside_open_unit_interval"
  raw_results$orientation[invalid_moment] <- NA_character_

  estimate_and_interval_columns <- c(
    "p_hat", "Se_hat", "Sp_hat", "alpha_hat", "gamma_hat",
    "p_lower", "p_upper", "Se_lower", "Se_upper",
    "Sp_lower", "Sp_upper", "alpha_lower", "alpha_upper",
    "gamma_lower", "gamma_upper"
  )
  for (column_name in estimate_and_interval_columns) {
    raw_results[[column_name]][invalid_moment | near_boundary] <- NA_real_
  }

  message(
    "Fast reprocessing corrections: ",
    sum(invalid_moment), " nonpositive/near-zero moment-variance datasets; ",
    sum(near_boundary), " mathematically boundary datasets."
  )

  bfgs_keys <- paste(
    raw_results$scenario_id,
    raw_results$replicate,
    sep = ":"
  )
  bfgs_indexes <- which(
    raw_results$method == "Numerical BFGS" &
      bfgs_keys %in% affected_dataset_keys
  )
  message(
    "Rerunning BFGS without a closed-form start for ",
    length(bfgs_indexes), " affected datasets."
  )

  for (position in seq_along(bfgs_indexes)) {
    row_index <- bfgs_indexes[position]
    y <- as.numeric(unlist(
      raw_results[row_index, c("y0", "y1", "y2", "y3")],
      use.names = FALSE
    ))
    error_message <- NA_character_
    fit <- tryCatch(
      mle_3rep_grad_logit_multistart(y, closed_form_init = NULL),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    bg <- extract_bfgs(fit, error_message)

    raw_results$admissible[row_index] <- bg$admissible
    raw_results$failure_reason[row_index] <- bg$failure_reason
    raw_results$orientation[row_index] <- bg$orientation
    raw_results$p_hat[row_index] <- bg$p_hat
    raw_results$Se_hat[row_index] <- bg$Se_hat
    raw_results$Sp_hat[row_index] <- bg$Sp_hat
    raw_results$alpha_hat[row_index] <- 1 - bg$Se_hat
    raw_results$gamma_hat[row_index] <- bg$Sp_hat
    raw_results$p_lower[row_index] <- bg$p_lower
    raw_results$p_upper[row_index] <- bg$p_upper
    raw_results$Se_lower[row_index] <- bg$Se_lower
    raw_results$Se_upper[row_index] <- bg$Se_upper
    raw_results$Sp_lower[row_index] <- bg$Sp_lower
    raw_results$Sp_upper[row_index] <- bg$Sp_upper
    raw_results$alpha_lower[row_index] <- 1 - bg$Se_upper
    raw_results$alpha_upper[row_index] <- 1 - bg$Se_lower
    raw_results$gamma_lower[row_index] <- bg$Sp_lower
    raw_results$gamma_upper[row_index] <- bg$Sp_upper
    raw_results$optimizer_convergence[row_index] <-
      bg$optimizer_convergence
    raw_results$point_estimation_success[row_index] <-
      bg$point_estimation_success
    raw_results$inference_success[row_index] <- bg$inference_success
    raw_results$failure_stage[row_index] <- bg$failure_stage
    raw_results$selected_start[row_index] <- bg$selected_start
    raw_results$n_starts_attempted[row_index] <- bg$n_starts_attempted
    raw_results$closed_form_start_available[row_index] <- FALSE
    raw_results$loglik[row_index] <- bg$loglik

    if (position %% 100L == 0L || position == length(bfgs_indexes)) {
      message(
        "  Reprocessed affected BFGS dataset ",
        position, "/", length(bfgs_indexes), "."
      )
    }
  }

  raw_results
}

reprocess_input_rds <- Sys.getenv("BFGS_REPROCESS_INPUT_RDS", "")

# -----------------------------
# Simulation or fast reprocessing
# -----------------------------
ncores_requested <- read_integer_environment("BFGS_RERUN_CORES", 0L)
detected_cores <- suppressWarnings(parallel::detectCores())
if (!is.finite(detected_cores)) detected_cores <- 1L
ncores <- if (ncores_requested > 0L) {
  ncores_requested
} else {
  max(1L, detected_cores - 2L)
}
ncores <- min(ncores, nrow(param_grid))
message("Output directory: ", output_dir)

run_one <- function(i) {
  row <- param_grid[i, ]
  checkpoint_file <- file.path(
    checkpoint_dir,
    sprintf("scenario_%03d.rds", row$scenario_id)
  )
  if (file.exists(checkpoint_file)) {
    checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (
      !is.null(checkpoint) &&
        nrow(checkpoint) == 2L * nsim &&
        all(checkpoint$scenario_id == row$scenario_id)
    ) {
      message("Scenario ", row$scenario_id, " loaded from checkpoint.")
      return(checkpoint)
    }
  }
  message(
    "Scenario ", i, "/", nrow(param_grid),
    ": p=", row$p_true,
    ", Se=", row$Se_true,
    ", Sp=", row$Sp_true,
    ", N=", row$N
  )
  result <- run_scenario_replicates(
    scenario_id = row$scenario_id,
    p_true = row$p_true,
    Se_true = row$Se_true,
    Sp_true = row$Sp_true,
    N = row$N,
    nsim = nsim,
    base_seed = base_seed
  )
  saveRDS(result, checkpoint_file, compress = FALSE)
  result
}

if (nzchar(reprocess_input_rds)) {
  reprocess_input_rds <- normalizePath(
    reprocess_input_rds,
    winslash = "/",
    mustWork = TRUE
  )
  message("Fast reprocessing from: ", reprocess_input_rds)
  raw_results <- readRDS(reprocess_input_rds)
  raw_results <- reprocess_existing_results(raw_results)
} else {
  message("Running ", nrow(param_grid), " scenarios with ", nsim,
          " replicates each on ", ncores, " cores.")
  if (.Platform$OS.type == "windows") {
    # Safe sequential fallback. For Windows parallelism, use a PSOCK cluster.
    raw_list <- lapply(seq_len(nrow(param_grid)), run_one)
  } else {
    raw_list <- mclapply(
      seq_len(nrow(param_grid)),
      run_one,
      mc.cores = ncores,
      mc.preschedule = TRUE
    )
  }

  raw_results <- bind_rows(raw_list)
}

saveRDS(
  raw_results,
  file.path(output_dir, "simulation_replicate_results.rds"),
  compress = "xz"
)
write.csv(
  raw_results,
  file.path(output_dir, "simulation_replicate_results.csv"),
  row.names = FALSE
)

# ============================================================
# Scenario summaries: bias, SD, MSE, coverage, inadmissibility
# ============================================================

long_results <- bind_rows(
  raw_results |>
    transmute(
      scenario_id, replicate, method, admissible, failure_reason,
      p_true, alpha_true, gamma_true, Se_true, Sp_true, N,
      parameter = "p", true_value = p_true,
      estimate = p_hat, lower = p_lower, upper = p_upper
    ),
  raw_results |>
    transmute(
      scenario_id, replicate, method, admissible, failure_reason,
      p_true, alpha_true, gamma_true, Se_true, Sp_true, N,
      parameter = "alpha", true_value = alpha_true,
      estimate = alpha_hat, lower = alpha_lower, upper = alpha_upper
    ),
  raw_results |>
    transmute(
      scenario_id, replicate, method, admissible, failure_reason,
      p_true, alpha_true, gamma_true, Se_true, Sp_true, N,
      parameter = "gamma", true_value = gamma_true,
      estimate = gamma_hat, lower = gamma_lower, upper = gamma_upper
    )
) |>
  mutate(
    error = estimate - true_value,
    squared_error = error^2,
    covered = admissible & is.finite(lower) & is.finite(upper) &
      lower <= true_value & upper >= true_value
  )

scenario_summary <- long_results |>
  group_by(
    scenario_id, method, parameter,
    p_true, alpha_true, gamma_true, Se_true, Sp_true, N
  ) |>
  summarise(
    nsim = n(),
    n_admissible = sum(admissible),
    inadmissibility_rate = mean(!admissible),
    mean_estimate = if (any(admissible)) mean(estimate[admissible]) else NA_real_,
    bias = if (any(admissible)) mean(error[admissible]) else NA_real_,
    empirical_sd = if (sum(admissible) > 1L) sd(estimate[admissible]) else NA_real_,
    mse = if (any(admissible)) mean(squared_error[admissible]) else NA_real_,
    coverage = if (any(admissible)) mean(covered[admissible]) else NA_real_,
    .groups = "drop"
  )

# ============================================================
# 2x2 admissibility cross-tabs
# ============================================================

admissibility_wide <- raw_results |>
  select(
    scenario_id, replicate, p_true, alpha_true, gamma_true, N,
    method, admissible
  ) |>
  pivot_wider(names_from = method, values_from = admissible)

admissibility_crosstab <- admissibility_wide |>
  mutate(
    status = case_when(
      `Closed-form` & `Numerical BFGS` ~ "Both admissible",
      `Closed-form` & !`Numerical BFGS` ~ "Closed-form only",
      !`Closed-form` & `Numerical BFGS` ~ "BFGS only",
      TRUE ~ "Neither admissible"
    )
  ) |>
  count(
    scenario_id, p_true, alpha_true, gamma_true, N, status,
    name = "frequency"
  ) |>
  group_by(scenario_id) |>
  mutate(proportion = frequency / sum(frequency)) |>
  ungroup()

admissibility_overall <- admissibility_wide |>
  mutate(
    status = case_when(
      `Closed-form` & `Numerical BFGS` ~ "Both admissible",
      `Closed-form` & !`Numerical BFGS` ~ "Closed-form only",
      !`Closed-form` & `Numerical BFGS` ~ "BFGS only",
      TRUE ~ "Neither admissible"
    )
  ) |>
  count(status, name = "frequency") |>
  mutate(proportion = frequency / sum(frequency))

# ============================================================
# Failure reasons
# ============================================================

failure_reasons <- raw_results |>
  filter(!admissible) |>
  count(
    method, p_true, alpha_true, gamma_true, N,
    failure_reason, name = "frequency"
  ) |>
  group_by(method, p_true, alpha_true, gamma_true, N) |>
  mutate(proportion = frequency / sum(frequency)) |>
  ungroup()

failure_reasons_overall <- raw_results |>
  filter(!admissible) |>
  count(method, failure_reason, name = "frequency") |>
  group_by(method) |>
  mutate(proportion = frequency / sum(frequency)) |>
  ungroup()

bfgs_start_summary <- raw_results |>
  filter(method == "Numerical BFGS") |>
  count(
    selected_start, n_starts_attempted, admissible,
    failure_stage, failure_reason,
    name = "frequency"
  ) |>
  mutate(proportion = frequency / sum(frequency)) |>
  arrange(desc(frequency))

bfgs_stage_summary <- raw_results |>
  filter(method == "Numerical BFGS") |>
  count(
    point_estimation_success, inference_success,
    failure_stage, failure_reason,
    name = "frequency"
  ) |>
  mutate(proportion = frequency / sum(frequency)) |>
  arrange(desc(frequency))

# ============================================================
# Matched head-to-head comparisons when both methods work
# ============================================================

matched_wide <- long_results |>
  filter(admissible) |>
  select(
    scenario_id, replicate, p_true, alpha_true, gamma_true, N,
    parameter, true_value, method, estimate, lower, upper, covered
  ) |>
  pivot_wider(
    names_from = method,
    values_from = c(estimate, lower, upper, covered),
    names_sep = "__"
  ) |>
  filter(
    !is.na(.data[["estimate__Closed-form"]]),
    !is.na(.data[["estimate__Numerical BFGS"]])
  ) |>
  mutate(
    error_cf = .data[["estimate__Closed-form"]] - true_value,
    error_bfgs = .data[["estimate__Numerical BFGS"]] - true_value,
    absolute_error_cf = abs(error_cf),
    absolute_error_bfgs = abs(error_bfgs),
    squared_error_cf = error_cf^2,
    squared_error_bfgs = error_bfgs^2,
    estimate_difference =
      .data[["estimate__Closed-form"]] -
      .data[["estimate__Numerical BFGS"]],
    absolute_error_difference =
      absolute_error_cf - absolute_error_bfgs,
    squared_error_difference =
      squared_error_cf - squared_error_bfgs,
    coverage_difference =
      as.integer(.data[["covered__Closed-form"]]) -
      as.integer(.data[["covered__Numerical BFGS"]])
  )

matched_summary <- matched_wide |>
  group_by(
    scenario_id, parameter, p_true, alpha_true, gamma_true, N
  ) |>
  summarise(
    n_both_admissible = n(),
    mean_estimate_difference = mean(estimate_difference),
    bias_closed_form = mean(error_cf),
    bias_BFGS = mean(error_bfgs),
    bias_difference = bias_closed_form - bias_BFGS,
    MAE_closed_form = mean(absolute_error_cf),
    MAE_BFGS = mean(absolute_error_bfgs),
    MAE_difference = mean(absolute_error_difference),
    MSE_closed_form = mean(squared_error_cf),
    MSE_BFGS = mean(squared_error_bfgs),
    MSE_difference = mean(squared_error_difference),
    coverage_closed_form = mean(.data[["covered__Closed-form"]]),
    coverage_BFGS = mean(.data[["covered__Numerical BFGS"]]),
    coverage_difference = mean(coverage_difference),
    MCSE_estimate_difference = sd(estimate_difference) / sqrt(n()),
    MCSE_MAE_difference = sd(absolute_error_difference) / sqrt(n()),
    MCSE_MSE_difference = sd(squared_error_difference) / sqrt(n()),
    MCSE_coverage_difference = sd(coverage_difference) / sqrt(n()),
    .groups = "drop"
  )

matched_overall <- matched_wide |>
  group_by(parameter) |>
  summarise(
    n_both_admissible = n(),
    median_absolute_estimate_difference =
      median(abs(estimate_difference), na.rm = TRUE),
    mean_MAE_difference = mean(absolute_error_difference, na.rm = TRUE),
    mean_MSE_difference = mean(squared_error_difference, na.rm = TRUE),
    mean_coverage_difference = mean(coverage_difference, na.rm = TRUE),
    proportion_closed_form_lower_absolute_error =
      mean(absolute_error_difference < 0, na.rm = TRUE),
    proportion_BFGS_lower_absolute_error =
      mean(absolute_error_difference > 0, na.rm = TRUE),
    proportion_equal_absolute_error =
      mean(absolute_error_difference == 0, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Figures
# ============================================================

plot_df <- scenario_summary |>
  mutate(
    parameter = factor(
      parameter,
      levels = c("p", "alpha", "gamma"),
      labels = c("p", "alpha", "gamma")
    ),
    method = factor(
      method,
      levels = c("Closed-form", "Numerical BFGS")
    ),
    N_factor = factor(N, levels = N_vals),
    prevalence = factor(p_true, levels = p_vals),
    alpha_panel = factor(alpha_true, levels = rev(sort(unique(alpha_true)))),
    gamma_panel = factor(gamma_true, levels = sort(unique(gamma_true)))
  )

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing = grid::unit(0.35, "lines"),
    plot.title = element_text(face = "bold")
  )

parameter_labels <- c(
  p = "prevalence p",
  alpha = "false-negative probability alpha",
  gamma = "specificity gamma"
)

make_performance_plot <- function(
    metric,
    parameter_name,
    y_label,
    title_prefix,
    reference = NULL,
    percent_axis = FALSE) {
  data <- plot_df |>
    filter(parameter == parameter_name)

  plot_object <- ggplot(
    data,
    aes(
      x = N_factor,
      y = .data[[metric]],
      color = prevalence,
      group = prevalence
    )
  ) +
    geom_point(size = 1.35) +
    geom_line(linewidth = 0.45) +
    facet_grid(method ~ alpha_panel + gamma_panel) +
    labs(
      x = "Sample size N",
      y = y_label,
      color = "True prevalence p",
      title = paste(title_prefix, parameter_labels[[parameter_name]]),
      subtitle = "Columns show alpha and gamma; all prevalence values are displayed together"
    ) +
    base_theme

  if (!is.null(reference)) {
    plot_object <- plot_object +
      geom_hline(yintercept = reference, linetype = "dashed", linewidth = 0.35)
  }
  if (percent_axis) {
    plot_object <- plot_object +
      scale_y_continuous(limits = c(0, 1), labels = percent)
  }
  plot_object
}

bias_plots <- lapply(
  c("p", "alpha", "gamma"),
  function(parameter_name) {
    make_performance_plot(
      "bias", parameter_name, "Bias", "Bias for",
      reference = 0
    )
  }
)
mse_plots <- lapply(
  c("p", "alpha", "gamma"),
  function(parameter_name) {
    make_performance_plot(
      "mse", parameter_name, "Mean squared error", "MSE for"
    )
  }
)
coverage_plots <- lapply(
  c("p", "alpha", "gamma"),
  function(parameter_name) {
    make_performance_plot(
      "coverage", parameter_name,
      "Coverage of nominal 95% confidence interval",
      "Coverage for",
      reference = 0.95,
      percent_axis = TRUE
    )
  }
)

inad_plot_df <- raw_results |>
  distinct(
    scenario_id, replicate, method, admissible,
    p_true, alpha_true, gamma_true, N
  ) |>
  group_by(method, p_true, alpha_true, gamma_true, N) |>
  summarise(
    inadmissibility_rate = mean(!admissible),
    .groups = "drop"
  ) |>
  mutate(
    method = factor(method, levels = c("Closed-form", "Numerical BFGS")),
    N_factor = factor(N, levels = N_vals),
    prevalence = factor(p_true, levels = p_vals),
    alpha_panel = factor(alpha_true, levels = rev(sort(unique(alpha_true)))),
    gamma_panel = factor(gamma_true, levels = sort(unique(gamma_true)))
  )

fig_inad <- ggplot(
  inad_plot_df,
  aes(
    x = N_factor,
    y = inadmissibility_rate,
    color = prevalence,
    group = prevalence
  )
) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.35) +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  facet_grid(method ~ alpha_panel + gamma_panel) +
  labs(
    x = "Sample size N",
    y = "Inadmissibility rate",
    color = "True prevalence p",
    title = "Inadmissibility rate",
    subtitle = "Columns show alpha and gamma; all prevalence values are displayed together"
  ) +
  base_theme

save_multipage_pdf <- function(plots, filename) {
  grDevices::pdf(
    file.path(output_dir, filename),
    width = 16,
    height = 9,
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  for (plot_object in plots) print(plot_object)
}

save_multipage_pdf(bias_plots, "Figure_bias.pdf")
save_multipage_pdf(mse_plots, "Figure_MSE.pdf")
save_multipage_pdf(coverage_plots, "Figure_coverage.pdf")
save_multipage_pdf(list(fig_inad), "Figure_inadmissibility.pdf")
save_multipage_pdf(
  c(bias_plots, mse_plots, coverage_plots, list(fig_inad)),
  "Figure_all_metrics_multi_prevalence.pdf"
)

# ============================================================
# Excel outputs
# ============================================================

wb <- createWorkbook()

sheets <- list(
  scenario_summary = scenario_summary,
  admissibility_crosstab = admissibility_crosstab,
  admissibility_overall = admissibility_overall,
  failure_reasons = failure_reasons,
  failure_reasons_overall = failure_reasons_overall,
  bfgs_start_summary = bfgs_start_summary,
  bfgs_stage_summary = bfgs_stage_summary,
  matched_summary = matched_summary,
  matched_overall = matched_overall,
  parameter_grid = param_grid
)

for (nm in names(sheets)) {
  write.csv(
    sheets[[nm]],
    file.path(output_dir, paste0(nm, ".csv")),
    row.names = FALSE
  )
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]])
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
}

saveWorkbook(
  wb,
  file.path(output_dir, "simulation_complete_analysis.xlsx"),
  overwrite = TRUE
)

message("\nDONE. Main outputs:")
message("  ", file.path(output_dir, "simulation_replicate_results.rds"))
message("  ", file.path(output_dir, "simulation_replicate_results.csv"))
message("  ", file.path(output_dir, "simulation_complete_analysis.xlsx"))
message("  Five multi-prevalence PDF figure files in ", output_dir)
