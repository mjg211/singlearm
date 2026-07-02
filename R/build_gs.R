#' Build a group sequential single-arm trial design
#'
#' Constructs an object of class \code{"sa_des_gs"} from user-supplied design
#' parameters, without performing an optimisation search. This allows bespoke
#' group sequential designs to be passed to downstream functions such as
#' \code{\link{opchar_gs}}, \code{\link{est_gs}}, \code{\link{pval_gs}}, and
#' \code{\link{ci_gs}}.
#'
#' @param J The maximal number of stages. Must be a positive integer.
#' @param n A vector of stage sample sizes of length \code{J}.
#' @param a A vector of futility boundaries of length \code{J}. Use
#'   \code{-Inf} to suppress early stopping for futility at a given stage.
#' @param r A vector of efficacy boundaries of length \code{J}. Use
#'   \code{Inf} to suppress early stopping for efficacy at a given stage.
#' @param pi0 The (undesirable) response probability used in the definition of
#'   the null hypothesis.
#' @param pi1 The (desirable) response probability at which the trial is
#'   powered.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param summary A logical variable indicating whether a summary of the
#'   function's progress should be printed to the console.
#' @return A list of class \code{"sa_des_gs"} containing the following
#'   elements
#'   \itemize{
#'     \item A list in the slot \code{$des} containing the design parameters.
#'     \item Each of the input variables as specified.
#'   }
#' @examples
#' des <- build_gs(J = 2, n = c(10, 19), a = c(1, 5), r = c(Inf, 6))
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}},
#'   \code{\link{est_gs}}, \code{\link{pval_gs}}, \code{\link{ci_gs}},
#'   and their associated \code{plot} family of functions.
#' @export
build_gs <- function(J = 2, n = c(10, 19), a = c(1, 5), r = c(Inf, 6),
                     pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                     summary = TRUE) {

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
