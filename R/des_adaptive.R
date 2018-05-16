#' Design an adaptive two-stage single-arm trial for a single binary endpoint
#'
#' Determines adaptive two-stage single-arm clinical trial designs for a single
#' binary primary endpoint.
#'
#' \code{des_adaptive()} supports the determination of adaptive two-stage
#' single-arm clinical trial designs for a single binary primary endpoint. For
#' all supported designs, the following hypotheses are tested for the response
#' probability \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}
#'
#' \ifelse{html}{\out{<center><i>H</i><sub>0</sub> : <i>&pi;</i> = <i>&pi;
#' </i><sub>0</sub>, <i>H</i><sub>1</sub> : <i>&pi;</i> = <i>&pi;</i><sub>
#' 1</sub>,</center>}}{\deqn{H_0 : \pi = \pi_0,\qquad H_1 : \pi = \pi_1,}}
#'
#' for \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\eqn{\pi_0}},
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\eqn{\pi_1}}, satisfying
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub> &lt;
#' <i>&pi;</i><sub>1</sub>}}{\eqn{\pi_0 < \pi_1}}, are specified using the
#' arguments \code{pi0} and \code{pi1}.
#'
#' In each instance, the optimal design is required to meet the following
#' operating characteristics
#'
#' \ifelse{html}{\out{<center><i>P</i>(<i>&pi;</i><sub>0</sub>) &le;
#' <i>&alpha;</i>, <i>P</i>(<i>&pi;</i><sub>1</sub>) &ge; 1 - <i>&beta;</i>,
#' </center>}}{\deqn{P(\pi_0) \le \alpha,\qquad P(\pi_1) \ge 1 - \beta,}}
#'
#' where \ifelse{html}{\out{<i>P</i>(<i>&pi;</i>)}}{\eqn{P(\pi)}} is the
#' probability of rejecting \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}
#' when the true response probability is
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}, and the values of
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}} and
#' \ifelse{html}{\out{<i>&beta;</i>}}{\eqn{\beta}} are specified using the
#' arguments \code{alpha} and \code{beta} respectively.
#'
#' An adaptive two-stage single-arm design for a single binary endpoint, is
#' then indexed by values for
#' \ifelse{html}{\out{<i>n</i><sub>1</sub>}}{\eqn{n_1}}, and two vectors:
#' \ifelse{html}{\out{<b><i>a</i></b><sub>2</sub> =
#' (<i>a</i><sub>20</sub>,&hellip;,<i>a</i><sub><i>2<i>n</i><sub>1</sub>
#' </i></sub>)}}{\eqn{\bold{a}_2=
#' (a_{20},\dots,a_{2n_1})}} and \ifelse{html}{\out{<b><i>n</i></b><sub>2</sub>
#' = (<i>n</i><sub>20</sub>,&hellip;,<i>n</i><sub><i>2<i>n</i><sub>1</sub>
#' </i></sub>)}}{\eqn{\bold{n}_2=(n_{20},\dots,n_{2n_1})}}.
#'
#' The purpose of this function is then to optimise the above parameters,
#' accounting for the chosen restrictions placed on these vectors, and the
#' chosen optimality criteria.
#'
#' The arguments \code{Nmin} and \code{Nmax} allow restrictions
#' to be placed on \ifelse{html}{\out{<i>n</i><sub>1</sub>}}{\eqn{n_1}} and
#' \ifelse{html}{\out{<b><i>n</i></b><sub>2</sub>}}{\eqn{\bold{n}_2}}.
#' Precisely, \code{Nmin} and \code{Nmax} set an inclusive range of allowed
#' values for the possible minimal and maximal trial sample sizes.
#'
#' In addition, \code{monotonic} also allows restrictions to be placed on
#' \ifelse{html}{\out{<b><i>n</i></b><sub>2</sub>}}{\eqn{\bold{n}_2}}.
#' Specifically, if \code{monotonic = TRUE}, the values in
#' \ifelse{html}{\out{<b><i>n</i></b><sub>2</sub>}}{\eqn{\bold{n}_2}} must be
#' monotonically decreasing.
#'
#' To describe the supported optimality criteria, denote the expected sample
#' size and median required sample size when the true response probability is
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} by
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i>)}}{\eqn{ESS(\pi)}} and
#' \ifelse{html}{\out{<i>Med</i>(<i>&pi;</i>)}}{\eqn{Med(\pi)}} respectively.
#' Then, the following optimality criteria are currently supported:
#' \itemize{
#' \item \code{"minimax"}: The design which minimises the maximal possible
#' sample size.
#' \item \code{"null_ess"}: The design which minimises
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i><sub>0</sub>)}}{\eqn{ESS(\pi_0)}}.
#' \item \code{"alt_ess"}: The design which minimises
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i><sub>1</sub>)}}{\eqn{ESS(\pi_1)}}.
#' }
#'
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desirable) response probability used in the definition of the
#' alternative hypothesis.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param Nmin The minimal total sample size to allow in considered
#' designs.
#' @param Nmax The maximal total sample size to allow in considered designs.
#' @param optimality Choice of optimal design criteria. Must be one of
#' \code{"null_ess"}, \code{"alt_ess"}, or \code{"minimax"}.
#' @param monotonic A logical variable indicating whether the second stage
#' sample sizes must be monotonically decreasing across their positive range.
#' @param summary A logical variable indicating a summary of the function's
#' progress should be printed to the console.
#' @return A list of class \code{"sa_des_adaptive"} containing the following
#' elements
#' \itemize{
#' \item A list in the slot \code{$des} containing details of the identified
#' optimal design.
#' \item A tibble in the slot \code{$feasible}, consisting of the
#' identified designs which met the required operating characteristics.
#' \item Each of the input variables as specified.
#' }
#' @examples
#' # The null-optimal design for the default parameters
#' null_ess <- des_adaptive()
#' @seealso \code{\link{opchar_adaptive}},  and their associated \code{plot}
#' family of functions.
#' @export
des_adaptive <- function(pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                         Nmin = 1, Nmax = 30, optimality = "null_ess",
                         monotonic = F, summary = F){

  ##### Input Checking #########################################################

  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_integer_pair_range(Nmin, Nmax, "Nmin", "Nmax", c(0, Inf))
  check_belong(optimality, "optimality", c("minimax", "null_ess", "alt_ess"),
               "1")
  check_logical(monotonic, "monotonic")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  ##### Main Computations ######################################################

  if (optimality == "alt_ess") {
    w <- c(0, 1, 0)
  } else {
    w <- c(1, 0, 0)
  }
  if (optimality != "minimax") {
    poss_n1     <- 1:(Nmax - 1)
    possible    <- matrix(NA, nrow = Nmax - 1, ncol = 14 + 3*Nmax)
    for (n1 in poss_n1) {
      possible[n1, c(1:(n1 + 2), (Nmax + 2):(Nmax + 4 + n1),
                     (4 + 2*Nmax):(4 + 2*Nmax + n1),
                     (ncol(possible) - 10):ncol(possible))] <-
        des_adaptive_n2forn1(pi0, pi1, alpha, beta, n1, max(Nmin - n1, 0),
                             Nmax - n1, monotonic, w)
      if (summary) {
        message("...optimal adaptive design search with n1 = ", n1, " completed...")
      }
    }
  } else {
    poss_n1n2max <- as.matrix(expand.grid(rep(list(1:(Nmax - 1), 2))))
    poss_n1n2max <- poss_n1n2max[which(rowSums(poss_n1n2max)  <= Nmax), ]
    possible     <- matrix(NA, nrow = nrow(poss_n1n2max), ncol = 14 + 3*Nmax)
    curr_min     <- Nmax
    for (i in 1:nrow(possible)) {
      if (curr_min >= sum(poss_n1n2max[i, ])) {
        n1    <- poss_n1n2max[i, 1]
        n2max <- poss_n1n2max[i, 2]
        possible[i, c(1:(n1 + 2), (Nmax + 2):(Nmax + 4 + n1),
                      (4 + 2*Nmax):(4 + 2*Nmax + n1),
                      (ncol(possible) - 10):ncol(possible))] <-
          des_adaptive_n2forn1(pi0, pi1, alpha, beta, n1, 0, n2max, monotonic,
                               w)
        if (possible[i, 1] != 0) {
          curr_min <- possible[i, 14 + 3*Nmax]
        }
      }
      if (all(summary, i%%25 == 0)) {
        message("...optimal adaptive design search with n1 = ", n1, " and n2max = ", n2max, " completed...")
      }
    }
    w <- c(10^-6, 0, 1)
  }
  colnames(possible) <- c("n1", paste("n2(", 0:(Nmax - 1), ")", sep = ""),
                          "a1", "r1", paste("a2(", 0:(Nmax - 1), ")", sep = ""),
                          paste("r2(", 0:(Nmax - 1), ")", sep = ""), "P(pi0)",
                          "P(pi1)", "ESS(pi0)", "ESS(pi1)", "PET(pi0)",
                          "PET(pi1)", "Med(pi0)", "Med(pi1)", "VSS(pi0)",
                          "VSS(pi1)", "max(N)")
  possible           <- tibble::as_tibble(possible)
  feasible           <- dplyr::filter(possible, `P(pi0)` > 0)
  if (nrow(feasible) > 0) {
    feasible           <- dplyr::mutate(feasible, O = w[1]*`ESS(pi0)` +
                                          w[2]*`ESS(pi1)` +
                                          w[3]*`max(N)`)
    feasible           <- dplyr::arrange(feasible, O)
    feasible[, c(1:(3*(Nmax + 1)), 14 + 3*Nmax)] <-
      suppressWarnings(dplyr::mutate_if(feasible[, c(1:(3*(Nmax + 1)),
                                                     14 + 3*Nmax)], is.double,
                       as.integer))
    max_feas_n1        <- max(feasible$n1)
    feasible           <- feasible[, c(1:(max_feas_n1 + 2),
                                       (Nmax + 2):(Nmax + 4 + max_feas_n1),
                                       (4 + 2*Nmax):(4 + 2*Nmax + max_feas_n1),
                                       (ncol(feasible) - 11):ncol(feasible))]
    opt_n1             <- as.numeric(feasible[1, 1])
    des                <- list(J = 2, n1 = opt_n1,
                               n2 = as.numeric(feasible[1, 2:(opt_n1 + 2)]),
                               a1 = as.numeric(feasible[1, max_feas_n1 + 3]),
                               r1 = as.numeric(feasible[1, max_feas_n1 + 4]),
                               a2 = as.numeric(feasible[1, (max_feas_n1 + 5):
                                                          (max_feas_n1 + 5 +
                                                             opt_n1)]),
                               r2 = as.numeric(feasible[1, (6 + 2*max_feas_n1):
                                                          (6 + 2*max_feas_n1 +
                                                             opt_n1)]),
                               pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                               opchar = feasible[1, c(1:(opt_n1 + 2),
                                                      (max_feas_n1 + 3):
                                                        (max_feas_n1 + 5 +
                                                           opt_n1),
                                                       (6 + 2*max_feas_n1):
                                                        (6 + 2*max_feas_n1 +
                                                           opt_n1),
                                                      (ncol(feasible) - 11):
                                                        ncol(feasible))])
  } else {
    feasible <- des <- NULL
    if (summary) {
      message("...no feasible designs found in range of considered maximal allowed sample size. Consider decreasing Nmin and increasing Nmax...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = feasible, pi0 = pi0, pi1 = pi1,
                        alpha = alpha, beta = beta, Nmin = Nmin, Nmax = Nmax,
                        optimality = optimality, monotonic = monotonic,
                        summary = summary)
  class(output) <- "sa_des_adaptive"
  return(output)
}
