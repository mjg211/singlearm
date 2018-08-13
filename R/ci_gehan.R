#' Determine confidence intervals in a group sequential single-arm trial design
#' for a single binary endpoint
#'
#' Determines possible confidence intervals at the end of a
#' group sequential single-arm trial for a single binary endpoint, as determined
#' using \code{des_gs()}. Support is available to compute confidence intervals
#' using the naive (\code{"naive"}), exact (\code{"exact"}), and Mid-p
#' (\code{"mid_p"}) approaches.
#'
#' In addition, the performance of the chosen confidence interval determination
#' procedures (including their coverage and expected length) for each value of
#' \ifelse{html}{\out{<i>pi</i>)}}{\eqn{\pi}} in the supplied vector
#' \ifelse{html}{\out{<b><i>pi</i></b>)}}{\eqn{\bold{\pi}}}, will also be
#' evaluated.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_gs"}, as returned by
#' \code{des_gs()}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and the maximum number of possible
#' stages.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the point estimation procedures at. This will internally
#' default to be the \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}
#' and \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from
#' \code{des} if it is left unspecified.
#' @param alpha Level to use in confidence interval construction. Defaults to
#' the value of \ifelse{html}{\out{<i>&alpha;</i>}}{\deqn{\alpha}} used in the
#' construction of \code{des} (i.e., \code{des$alpha}).
#' @param method A vector of methods to use to construct confidence intervals.
#' Currently, support is available to use the naive (\code{"naive"}), exact
#' (\code{"exact"}), and Mid-p (\code{"mid_p"}) approaches. Defaults to all
#' available methods.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_ci_gs"} containing the following elements
#' \itemize{
#' \item A tibble in the slot \code{$ci} summarising the possible confidence
#' intervals at the end of the trial for the supplied design, according to the
#' chosen methods.
#' \item A tibble in the slot \code{$perf} summarising the performance of the
#' chosen confidence interval determination procedures for each specified value
#' of \ifelse{html}{\out{<i>&pi;</i>}}{\deqn{\pi}}.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal two-stage design for the default parameters
#' des <- des_gs()
#' # Determine the performance of all supported confidence interval
#' # determination procedures for a range of possible response probabilities
#' ci  <- ci_gs(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}}, \code{\link{est_gs}},
#' \code{\link{pval_gs}}, and their associated \code{plot} family of functions.
#' @export
ci_gehan <- function(des, k, pi, alpha = des$alpha, summary = F) {

  ##### Input Checking #########################################################

  #check_sa_des_gehan(des, "des")
  if (!missing(k)) {
    check_k(k, des, NULL)
  } else {
    k <- 1:des$des$J
  }
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Confidence interval determination for group-sequential single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k \u2208 {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("The following methods will be used to construct confidence intervals\n")
    message("  \u2022 Naive.")
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  a1 <- des$des$a1; r1 <- des$des$r1; n1 <- des$des$n1; n2 <- des$des$n2
  terminal   <- terminal_states_adaptive(a1, r1, n1, n2, k)
  ci         <- tibble::tibble(s = terminal$s, m = terminal$m,
                               k = factor(terminal$k, k), `clow(s,m)` = NA,
                               `cupp(s,m)` = NA)
  for (i in 1:nrow(ci)) {
    ci[i, 4:5] <- ci_fixed_clopper_pearson(ci$s[i], ci$m[i], alpha)
    if (all(summary, i%%100 == 0)) {
      message("... ", i, " confidence intervals determined...")
    }
  }
  ci     <- dplyr::mutate(ci, `l(s,m)` = `cupp(s,m)` - `clow(s,m)`)
  pmf    <- pmf_adaptive(pi, a1, r1, n1, n2, k)
  len_pi <- length(pi)
  perf   <- tibble::tibble(pi = pi, `bar(L)` = NA, `max(L)` = NA,
                           `E(L|pi)` = NA, `Var(L|pi)` = NA, `Cover(C|pi)` = 0)
  for (i in 1:nrow(perf)) {
    pmf_i                   <- dplyr::filter(pmf, pi == perf$pi[i])
    perf$`bar(L)`[i]        <- mean(ci$`l(s,m)`)
    perf$`max(L)`[i]        <- max(ci$`l(s,m)`)
    perf$`E(L|pi)`[i]       <- sum(pmf_i$`f(s1,s,m|pi)`*ci$`l(s,m)`)
    perf$`Var(L|pi)`[i]     <- sum(pmf_i$`f(s1,s,m|pi)`*ci$`l(s,m)`^2) -
                                 perf$`E(L|pi)`[i]^2
    coverage                <- dplyr::filter(ci, `clow(s,m)` <= perf$pi[i] &
                                                   `cupp(s,m)` >= perf$pi[i])
    if (nrow(coverage) > 0) {
      perf$`Cover(C|pi)`[i] <- sum(dplyr::filter(pmf_i,
                                                 s %in% coverage$s)$`f(s1,s,m|pi)`)
    }
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " values of pi evaluated...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(ci = ci, perf = perf, des = des, k = k, pi = pi,
                        alpha = alpha, summary = summary)
  class(output) <- "sa_ci_gehan"
  return(output)
}
