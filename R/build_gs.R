#' @export
build_gs <- function(J = 2, n = c(10, 19), a = c(1, 5), r = c(Inf, 6),
                     pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                     summary = T) {

  ##### Input checking #########################################################

  check_integer_range(J, "J", c(1, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_gs_boundaries(J, n, a, r)
  check_logical(summary, "summary")

  ##### Main computations ######################################################

  if (summary){
    message("Building the design...")
  }
  des           <- list(J = J, n = n, a = a, r = r, pi0 = pi0, pi1 = pi1,
                        alpha = alpha, beta = beta)

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = NULL, J = J, pi0 = pi0,
                        pi1 = pi1, alpha = alpha, beta = beta, Nmin = NULL,
                        Nmax = NULL, futility = any(is.finite(a[1:(J - 1)])),
                        efficacy = any(is.finite(r[1:(J - 1)])),
                        optimality = NULL, point_prior = NULL,
                        beta_prior = NULL, equal_n = NULL, ensign = NULL,
                        summary = summary)
  class(output) <- "sa_des_gs"
  return(output)

}
