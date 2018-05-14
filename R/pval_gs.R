#' Determine p-values in a group sequential single-arm trial design for a single
#' binary endpoint
#'
#' Determines possible p-values at the end of a
#' group sequential single-arm trial for a single binary endpoint, as determined
#' using \code{des_gs()}. Support is available to compute point estimates using
#' the naive (\code{"naive"}), MLE-ordering (\code{"mle"}), UMVUE-ordering
#' (\code{"umvue"}), and conditional (\code{"conditional"}) approaches.
#'
#' In addition, the performance of the chosen p-value calculation procedures
#' (including their expected value and variance) for each value of
#' \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in the supplied vector
#' \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}}, will also be
#' evaluated.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_gs"}, as returned by
#' \code{des_gs()}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and the maximum number of allowed stages
#' across the supplied designs. If left unspecified, it will internally default
#' to all possible stages.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the point estimation procedures at. This will internally
#' default to be the \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}
#' and \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from
#' \code{des} if it is left unspecified.
#' @param method A vector of methods to use to construct point estimates.
#' Currently, support is available to use the naive (\code{"naive"}),
#' MLE-ordering (\code{"mle"}), UMVUE-ordering (\code{"umvue"}), and conditional
#' (\code{"conditional"}), approaches. Defaults to all available methods.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_pval_gs"} containing the following elements
#' \itemize{
#' \item A tibble in the slot \code{$pval} summarising the possible p-values at
#' the end of the trial for the supplied design, according to the chosen
#' methods.
#' \item A tibble in the slot \code{$perf} summarising the performance of the
#' chosen p-value calculation procedures for each specified value of
#' \ifelse{html}{\out{<i>&pi;</i>}}{\deqn{\pi}}.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal two-stage design for the default parameters
#' des  <- des_gs()
#' # Determine the performance of all supported p-value calculation procedures
#' # for a range of possible response probabilities
#' pval <- pval_gs(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}}, \code{\link{est_gs}},
#' \code{\link{ci_gs}}, and their associated \code{plot} family of functions.
#' @export
pval_gs <- function(des, k, pi, method = c("conditional", "mle", "naive",
                                           "umvue"), summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_gs(des, "des")
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
  check_belong(method, "method", c("naive", "mle", "umvue", "conditional"),
               "any")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("p-value determination for group-sequential single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k \u2208 {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("You have chosen to use the following methods to construct confidence intervals\n")
    if ("naive" %in% method) {
      message("  • Naive.")
    }
    if ("mle" %in% method) {
      message("  • MLE-ordering.")
    }
    if ("umvue" %in% method) {
      message("  • UMVUE-ordering.")
    }
    if ("conditional" %in% method) {
      message("  • Conditional.")
    }
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  J          <- des$des$J; a <- des$des$a; r <- des$des$r; n <- des$des$n
  pi0        <- des$des$pi0
  terminal   <- terminal_states_gs(J, a, r, n, k)
  len_method <- length(method)
  pval     <- tibble::tibble(s = rep(terminal$s, each = len_method),
                             m = rep(terminal$m, each = len_method),
                             k = factor(rep(terminal$k, each = len_method), k),
                             method = factor(rep(method, nrow(terminal)),
                                             levels = method), `p(s,m)` = NA)
  for (i in 1:nrow(pval)) {
    pval$`p(s,m)`[i] <-
      switch(as.character(pval$method[i]),
             naive       = pval_fixed_exact(pval$s[i], pval$m[i], pi0),
             mle         = pval_gs_mle(pi0, pval$s[i], pval$m[i], J, a, r, n),
             umvue       = pval_gs_umvue(pi0, pval$s[i], pval$m[i],
                                         as.numeric(pval$k[i]), J, a, r, n),
             conditional = pval_gs_cond(pi0, pval$s[i], pval$m[i],
                                        as.numeric(pval$k[i]), J, a, r, n))
    if (all(summary, i%%100 == 0)) {
      message("... ", i, " p-values determined...")
    }
  }
  pmf    <- pmf_gs(pi, J, a, r, n, k)
  len_pi <- length(pi)
  perf   <- tibble::tibble(pi = rep(pi, len_method),
                           method = factor(rep(method, each = len_pi), method),
                           `E(p|pi)` = NA, `Var(p|pi)` = NA)
  for (i in 1:nrow(perf)) {
    pmf_i               <- dplyr::filter(pmf, pi == perf$pi[i])
    pval_i              <- dplyr::filter(pval, method == perf$method[i])
    perf$`E(p|pi)`[i]   <- sum(pmf_i$`f(s,m|pi)`*pval_i$`p(s,m)`)
    perf$`Var(p|pi)`[i] <- sum(pmf_i$`f(s,m|pi)`*pval_i$`p(s,m)`^2) -
                             perf$`E(p|pi)`[i]^2
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " pi-method combinations evaluated...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(pval = pval, perf = perf, des = des, k = k, pi = pi,
                        method = method, summary = summary)
  class(output) <- "sa_pval_gs"
  return(output)
}
