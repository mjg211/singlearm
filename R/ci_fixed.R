#' Determine confidence intervals in a single-stage single-arm trial design for
#' a single binary endpoint
#'
#' Determines all possible confidence intervals at the end of a single-stage
#' single-arm trial for a single binary endpoint, as determined
#' using \code{des_fixed()}. Support is available to compute confidence
#' intervals using the Agresti-Coull, Clopper-Pearson, Jeffreys, Mid-p, Wald,
#' and Wilson Score approaches.
#'
#' In addition, the performance of the chosen confidence interval determination
#' procedures (including their coverage and expected length) for each value of
#' \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in the supplied vector
#' \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}}, will also be
#' evaluated.
#'
#' @param des An object of class \code{"sa_des_fixed"}, as returned by
#' \code{des_fixed()}.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the confidence interval calculation procedure at. This will
#' internally default to be the
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from \code{des} if
#' it is left unspecified.
#' @param alpha \ifelse{html}{\out{<i>&alpha;</i>}}{\deqn{\alpha}}-level to use
#' in confidence interval construction. Defaults to the value of
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\deqn{\alpha}} used in the construction
#' of \code{des} (i.e., \code{des$alpha}).
#' @param method A vector of methods to use to construct p-values. Currently,
#' support is available to use the Agresti-Coull (\code{"agresti_coull"}),
#' Clopper-Pearson (\code{"clopper_pearson"}), Jeffreys (\code{"jeffreys"}),
#' Mid-p (\code{"mid_p"}), Wald (\code{"wald"}), and Wilson Score
#' (\code{"wilson"}) approaches to confidence interval determination.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_ci_fixed"} containing the following
#' elements
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
#' # Find the optimal single-stage design for the default parameters
#' des <- des_fixed()
#' # Determine the performance of all supported confidence interval
#' # determination procedures for a range of possible response probabilities
#' ci  <- ci_fixed(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, and their associated
#' \code{plot} family of functions.
#' @export
ci_fixed <- function(des, pi, alpha = des$alpha,
                     method = c("agresti_coull", "clopper_pearson", "jeffreys",
                                "mid_p", "wald", "wilson"), summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_fixed(des, "des")
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  check_belong(method, "method", c("agresti_coull", "clopper_pearson",
                                   "jeffreys", "mid_p", "wald", "wilson"),
               "any")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Confidence interval determination for single-stage single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to use the following methods to construct confidence intervals\n")
    if ("aggresti_coull" %in% method) {
      message("  • Aggresti-Coull.")
    }
    if ("clopper_pearson" %in% method) {
      message("  • Clopper-Pearson.")
    }
    if ("jeffreys" %in% method) {
      message("  • Jeffreys.")
    }
    if ("mid_p" %in% method) {
      message("  • Mid-p.")
    }
    if ("wald" %in% method) {
      message("  • Wald.")
    }
    if ("wilson" %in% method) {
      message("  • Wilson.")
    }
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  ci <- tibble::tibble(s = rep(0:des$des$n, each = length(method)),
                       m = rep(des$des$n, (des$des$n + 1)*length(method)),
                       method = factor(rep(method, des$des$n + 1),
                                       levels = method), `clow(s,m)` = NA,
                       `cupp(s,m)` = NA)
  for (i in 1:nrow(ci)) {
    ci[i, 4:5] <- switch(as.character(ci$method[i]),
                         agresti_coull   = ci_fixed_agresti_coull(ci$s[i],
                                                                  ci$m[i],
                                                                  alpha),
                         clopper_pearson = ci_fixed_clopper_pearson(ci$s[i],
                                                                    ci$m[i],
                                                                    alpha),
                         jeffreys        = ci_fixed_jeffreys(ci$s[i], ci$m[i],
                                                             alpha),
                         mid_p           = ci_fixed_mid_p(ci$s[i], ci$m[i],
                                                          alpha),
                         wald            = ci_fixed_wald(ci$s[i], ci$m[i],
                                                         alpha),
                         wilson          = ci_fixed_wilson(ci$s[i], ci$m[i],
                                                           alpha))
    if (all(summary, i%%100 == 0)) {
      message("...", i, " confidence intervals determined...")
    }
  }
  ci     <- dplyr::mutate(ci, `l(s,m)` = `cupp(s,m)` - `clow(s,m)`)
  len_pi <- length(pi)
  pmf    <- pmf_fixed(pi, des$des$n)
  perf   <- tibble::tibble(pi = rep(pi, length(method)),
                           method = factor(rep(method, each = len_pi), method),
                           `bar(L)` = NA, `max(L)` = NA, `E(L|pi)` = NA,
                           `Var(L|pi)` = NA, `Cover(C|pi)` = 0)
  for (i in 1:nrow(perf)) {
    ci_i                    <- dplyr::filter(ci, method == perf$method[i])
    pmf_i                   <- dplyr::filter(pmf, pi == perf$pi[i])
    perf$`E(L|pi)`[i]       <- sum(pmf_i$`f(s,m|pi)`*ci_i$`l(s,m)`)
    perf$`Var(L|pi)`[i]     <- sum(pmf_i$`f(s,m|pi)`*ci_i$`l(s,m)`^2) -
                                 perf$`E(L|pi)`[i]^2
    perf$`bar(L)`[i]        <- mean(ci_i$`l(s,m)`)
    perf$`max(L)`[i]        <- max(ci_i$`l(s,m)`)
    coverage                 <- dplyr::filter(ci_i, `clow(s,m)` <= perf$pi[i] &
                                               `cupp(s,m)` >= perf$pi[i])
    if (nrow(coverage) > 0) {
      perf$`Cover(C|pi)`[i] <- sum(dplyr::filter(pmf_i,
                                                 s %in% coverage$s)$`f(s,m|pi)`)
    }
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " pi-method combinations evaluated...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(ci = ci, perf = perf, des = des, pi = pi, alpha = alpha,
                        summary = summary)
  class(output) <- "sa_ci_fixed"
  return(output)
}
