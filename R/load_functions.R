## Loads only the function definitions from simulation_closed_form_vs_bfgs.R,
## so the estimators can be used without starting the full simulation.
load_estimator_functions <- function(
    path = "R/simulation_closed_form_vs_bfgs.R", envir = globalenv()) {
  exprs <- parse(path)
  is_function_def <- vapply(exprs, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
  }, logical(1))
  for (e in exprs[is_function_def]) eval(e, envir = envir)
  invisible(sum(is_function_def))
}
