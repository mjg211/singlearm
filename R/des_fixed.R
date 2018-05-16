#' Design a single-stage single-arm trial for a single binary endpoint
#'
#' Determines single-stage single-arm clinical trial designs
#' for a single binary primary endpoint, using either exact binomial
#' calculations or a normal approximation approach.
#'
#' \code{des_fixed()} supports the determination of single-stage single-arm
#' clinical trial designs for a single binary primary endpoint. The following
#' hypotheses are tested for the response probability
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}
#'
#' \ifelse{html}{\out{<center><i>H</i><sub>0</sub> : <i>&pi;</i> &le; <i>&pi;
#' </i><sub>0</sub>, <i>H</i><sub>1</sub> : <i>&pi;</i> > <i>&pi;</i><sub>
#' 0</sub>,</center>}}{\deqn{H_0 : \pi \le \pi_0,\qquad H_1 : \pi > \pi_0,}}
#'
#' for \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\eqn{\pi_0}},  specified
#' using the argument \code{pi0}.
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
#' arguments \code{alpha} and \code{beta} respectively. Moreover,
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\eqn{\pi_1}}, satisfying
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub> &lt;
#' <i>&pi;</i><sub>1</sub>}}{\eqn{\pi_0 < \pi_1}}, is specified using the
#' argument \code{pi1}.
#'
#' A single-stage single-arm design for a single binary endpoint is ultimately
#' indexed by three parameters: \ifelse{html}{\out{<i>a</i>}}{\eqn{a}},
#' \ifelse{html}{\out{<i>r</i>}}{\eqn{r}}, and
#' \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}.
#'
#' With these parameters, and denoting the number of responses after
#' \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} outcomes have been observed by
#' \ifelse{html}{\out{<i>s</i><sub><i>m</i></sub>}}{\eqn{s_m}}, the testing
#' rules for the trial are as follows
#'
#' \itemize{
#' \item If \ifelse{html}{\out{<i>s</i><sub><i>n</i></sub>
#' &le; <i>a</i>}}{\eqn{s_{n} \le a}}, then do not reject
#' \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' \item Else if \ifelse{html}{\out{<i>s</i><sub><i>n</i></sub>
#' &ge; <i>r</i>}}{\eqn{s_{n} \ge r}}, then reject
#' \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' }
#'
#' The purpose of this function is then to determine (optimise)
#' \ifelse{html}{\out{<i>a</i>}}{\eqn{a}},
#' \ifelse{html}{\out{<i>r</i>}}{\eqn{r}}, and
#' \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}, accounting for the chosen
#' restrictions placed on these parameters.
#'
#' The arguments \code{Nmin} and \code{Nmax} allow restrictions
#' to be placed on \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}.
#' Precisely, \code{Nmin} and \code{Nmax} set an inclusive range of allowed
#' values for \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}.
#'
#' Note that to ensure a decision is made about
#' \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}, this function always
#' enforces \ifelse{html}{\out{<i>a</i> + 1 = <i>r</i>}}{\eqn{a + 1 = r}}.
#'
#' The optimal design is then the one that minimises
#' \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}. In the case where there are multiple
#' feasible designs with the same minimal value of
#' \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}, the optimal design is the one amongst
#' these which maximises
#' \ifelse{html}{\out{<i>P</i>(<i>&pi;</i><sub>1</sub>)}}{\eqn{P(\pi_1)}}.
#'
#' If \code{exact} is set to \code{TRUE} then exact binomial probability
#' calculations are used to identify the optimal design. Otherwise, a normal
#' approximation approach is used. Note that setting \code{exact = TRUE} is
#' recommended.
#'
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desirable) response probability at which the trial is
#' powered.
#' @param alpha The desired maximal type-I error rate.
#' @param beta The desired maximal type-II error rate.
#' @param Nmin The minimal sample size to allow in considered designs.
#' @param Nmax The maximal sample size to allow in considered designs.
#' @param exact A logical variable indicating whether exact binomial
#' calculations or a normal approximation approach should be used to determine
#' the optimal design.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_des_fixed"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$des} summarising the characteristics of the
#' identified optimal design. This will be \code{NULL} if no feasible design was
#' identified in the considered range for \ifelse{html}{\out{<i>n</i>}}{\eqn{n}}.
#' \item A tibble in the slot \code{$feasible}, consisting of the identified
#' designs which met the required operating characteristics.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # The optimal design for the default parameters
#' des <- des_fixed()
#' # Find the optimal single-stage design for a 10% type-I error rate
#' des_10      <- des_fixed(alpha = 0.1)
#' @references A'Hern RP (2001) Sample size tables for exact single-stage phase
#' II designs. \emph{Statistics in Medicine} \strong{20:}859-66.
#' @references Fleming TR (1982) One-sample multiple testing procedure for phase
#' II clinical trials. \emph{Biometrics} \strong{38:}143-51.
#' @seealso \code{\link{opchar_fixed}}, \code{\link{est_fixed}},
#' \code{\link{pval_fixed}}, \code{\link{ci_fixed}}, and their
#' associated \code{plot} family of functions. Note that similar functionality
#' is available through \code{\link[clinfun]{ph2single}}.
#' @export
des_fixed <- function(pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2, Nmin = 1,
                      Nmax = 50, exact = T, summary = F) {

  ##### Input Checking #########################################################

  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_integer_pair_range(Nmin, Nmax, "Nmin", "Nmax", c(0, Inf))
  check_logical(exact, "exact")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Design of a single-stage single-arm trial for a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("\nYou have chosen to test the following hypotheses\n")
    message("     H\u2080: \u03c0 \u2264 \u03c0\u2080 = ", pi0, ", H\u2081: \u03c0 > \u03c0\u2080 = ", pi0, ".\n")
    message("with the following error constraints\n")
    message("     P(\u03c0\u2080) = P(", pi0, ") \u2264 \u03b1 = ", alpha, ", P(\u03c0\u2081) = P(", pi1, ") \u2265 1 - \u03b2 = ", 1 - beta, ".\n")
    Sys.sleep(2)
    message("You have chosen to restrict the allowed possible sample size N = n such that\n")
    message("  \u2022 N \u2265 ", Nmin, ".")
    message("  \u2022 N \u2264 ", Nmax, ".\n")
    Sys.sleep(2)
    if (exact) {
      message("You have chosen to compute the optimal design using exact binomial calculations.\n")
    } else {
      message("You have chosen to compute the optimal design using a normal approximation approach.\n")
    }
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  pmf         <- list()
  for (n in 1:Nmax) {
    pmf[[n]] <- pmf_fixed(c(pi0, pi1), n)
  }
  pmf         <- tibble::as_tibble(plyr::rbind.fill(pmf))
  if (exact) {
    possible  <- as.matrix(expand.grid(rep(list(0:Nmax), 2)))[, 2:1]
  } else {
    possible  <- Nmin:Nmax
    possible  <- cbind(possible, round(possible*pi0 + stats::qnorm(1 - alpha)*
                                         sqrt(possible*pi0*(1 - pi0))))
  }
  possible    <- possible[which(possible[, 1] > possible[, 2] &
                                  possible[, 2] >= 0), ]
  possible   <- tibble::tibble(n = as.integer(possible[, 1]),
                               a = as.integer(possible[, 2]),
                               r = a + 1L, `P(pi0)` = NA, `P(pi1)` = NA)
  for (i in 1:nrow(possible)) {
    possible[i, 4:5] <-
      c(sum(dplyr::filter(pmf, pi == pi0 & m == as.integer(possible[i, 1]) &
                            s > as.integer(possible[i, 2]))$`f(s,m|pi)`),
        sum(dplyr::filter(pmf, pi == pi1 & m == as.integer(possible[i, 1]) &
                            s > as.integer(possible[i, 2]))$`f(s,m|pi)`))
    if (all(summary, i%%1000 == 0)) {
      message("...", i, " designs evaluated...")
    }
  }
  feasible           <- dplyr::filter(possible, `P(pi0)` <= alpha &
                                        `P(pi1)` >= 1 - beta)
  if (nrow(feasible) > 0) {
    if (summary) {
      message("...feasible designs identified in range of considered N...")
      Sys.sleep(2)
      message("...now identifying the optimal design...")
    }
    feasible         <- dplyr::mutate(feasible, `ESS(pi0)` = n,
                                      `ESS(pi1)` = n, `VSS(pi0)` = 0,
                                      `VSS(pi1)` = 0, `Med(pi0)` = n,
                                      `Med(pi1)` = n, `A1(pi0)` = 1 - `P(pi0)`,
                                      `R1(pi0)` = `P(pi0)`,
                                      `A1(pi1)` = 1 - `P(pi1)`,
                                      `R1(pi1)` = `P(pi1)`, `S1(pi0)` = 1,
                                      `S1(pi1)` = 1, `cum(S1(pi0))` = 1,
                                      `cum(S1(pi1))` = 1)
    feasible         <- dplyr::arrange(feasible, n, dplyr::desc(`P(pi1)`))
    feasible[, 4:19] <- dplyr::mutate_if(feasible[, 4:19], is.integer,
                                           as.double)
    des                <- list(n = as.numeric(feasible[1, 1]),
                               a = as.numeric(feasible[1, 2]),
                               r = as.numeric(feasible[1, 3]), pi0 = pi0,
                               pi1 = pi1, alpha = alpha, beta = beta,
                               opchar = feasible[1, ])
  } else {
    feasible <- des <- NULL
    if (summary) {
      message("...no feasible designs identified in range of considered N. Consider decreasing Nmin and increasing Nmax.")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = feasible, pi0 = pi0, pi1 = pi1,
                        alpha = alpha, beta = beta, Nmin = Nmin, Nmax = Nmax,
                        exact = exact, summary = summary)
  class(output) <- "sa_des_fixed"
  return(output)
}
