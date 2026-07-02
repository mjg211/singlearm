#' Build a curtailed group sequential single-arm trial design
#'
#' Constructs an object of class \code{"sa_des_curtailed"} from user-supplied
#' design parameters, without performing an optimisation search. This allows
#' bespoke curtailed group sequential designs to be passed to downstream
#' functions such as \code{\link{opchar_curtailed}} and
#' \code{\link{est_curtailed}}.
#'
#' @param J The maximal number of stages. Must be a positive integer.
#' @param n A vector of stage sample sizes of length \code{J}.
#' @param a A vector of futility boundaries of length \code{J}. Use
#'   \code{-Inf} to suppress early stopping for futility at a given stage.
#' @param r A vector of efficacy boundaries of length \code{J}. Use
#'   \code{Inf} to suppress early stopping for efficacy at a given stage.
#' @param a_curt A vector of curtailed futility boundaries of length
#'   \code{sum(n)}.
#' @param r_curt A vector of curtailed efficacy boundaries of length
#'   \code{sum(n)}.
#' @param pi0 The (undesirable) response probability used in the definition of
#'   the null hypothesis.
#' @param pi1 The (desirable) response probability at which the trial is
#'   powered.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param summary A logical variable indicating whether a summary of the
#'   function's progress should be printed to the console.
#' @return A list of class \code{"sa_des_curtailed"} containing the following
#'   elements
#'   \itemize{
#'     \item A list in the slot \code{$des} containing the design parameters.
#'     \item Each of the input variables as specified.
#'   }
#' @examples
#' des <- build_curtailed(J = 2, n = c(10, 19), a = c(1, 5), r = c(Inf, 6),
#'                        a_curt = c(rep(-Inf, 28), 5),
#'                        r_curt = c(rep(Inf, 28), 6))
#' @seealso \code{\link{des_curtailed}}, \code{\link{opchar_curtailed}},
#'   \code{\link{est_curtailed}}, and their associated \code{plot} family of
#'   functions.
#' @export
build_curtailed <- function(J = 2, n = c(10, 19), a = c(1, 5), r = c(Inf, 6),
                            a_curt, r_curt, pi0 = 0.1, pi1 = 0.3, alpha = 0.05,
                            beta = 0.2, summary = TRUE) {

  ##### Input checking #########################################################

  check_integer_range(J, "J", c(0, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_gs_boundaries(J, n, a, r)
  check_gs_boundaries(sum(n), rep(1, sum(n)), a_curt, r_curt)
  check_logical(summary, "summary")

  ##### Main computations ######################################################

  if (summary){
    message("Building the design...")
  }
  des           <- list(J = J, n = n, a = a, r = r, pi0 = pi0, pi1 = pi1,
                        alpha = alpha, beta = beta, J_curt = sum(n),
                        a_curt = a_curt, r_curt = r_curt,
                        n_curt = rep(1, sum(n)))

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
  class(output) <- "sa_des_curtailed"
  return(output)

}
