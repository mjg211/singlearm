#' Determine p-values in a single-stage single-arm trial design for a single
#' binary endpoint
#'
#' Determines all possible p-values at the end of a single-stage single-arm
#' trial for a single binary endpoint, as determined using \code{des_fixed()}.
#' Support is available to compute p-values using exact binomial tail
#' probabilities, or using a normal approximation.
#'
#' In addition, the performance of the chosen p-value calculation procedures
#' (including their expected value and variance) for each value of
#' \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in the supplied vector
#' \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}}, will also be
#' evaluated.
#'
#' @param des An object of class \code{"sa_des_fixed"}, as returned by
#' \code{des_fixed()}.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the p-value calculation procedure at. This will internally
#' default to be the \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}
#' and \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from
#' \code{des} if it is left unspecified.
#' @param method A vector of methods to use to construct p-values. Currently,
#' support is available to use exact binomial tail probabilities
#' (\code{"exact"}) and to use a normal approximation (\code{"normal"})
#' approach.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_pval_fixed"} containing the following
#' elements
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
#' # Find the optimal single-stage design for the default parameters
#' des  <- des_fixed()
#' # Determine the performance of both supported p-value calculation procedures
#' # for a range of possible response probabilities
#' pval <- pval_fixed(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{ci_fixed}}, and their associated
#' \code{plot} family of functions.
#' @export
pval_fixed <- function(des, pi, method = c("exact", "normal"), summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_fixed(des, "des")
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  check_belong(method, "method", c("exact", "normal"), "any")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("p-value determination for single-stage single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to use the following methods to construct p-values\n")
    if ("exact" %in% method) {
      message("  • Exact.")
    }
    if ("normal" %in% method) {
      message("  • Normal.")
    }
    Sys.sleep(2)
    message("\nBeginning the required calculations...")
  }

  ##### Main Computations ######################################################

  s    <- as.integer(0:des$des$n)
  n    <- as.integer(rep(des$des$n, des$des$n + 1))
  pval <- tibble::tibble(s = rep(0:des$des$n, each = length(method)),
                         m = rep(des$des$n, (des$des$n + 1)*length(method)),
                         method = factor(rep(method, des$des$n + 1),
                                         method), `p(s,m)` = NA)
  for (i in 1:nrow(pval)) {
    pval$`p(s,m)`[i] <- switch(as.character(pval$method[i]),
                               exact  = pval_fixed_exact(pval$s[i], pval$m[i],
                                                         des$des$pi0),
                               normal = pval_fixed_normal(pval$s[i], pval$m[i],
                                                          des$des$pi0))
    if (all(summary, i%%100 == 0)) {
      message("...", i, " p-values determined...")
    }
  }
  len_pi <- length(pi)
  pmf    <- pmf_fixed(pi, des$des$n)
  perf   <- tibble::tibble(pi = rep(pi, length(method)),
                           method = factor(rep(method, each = len_pi), method),
                           `E(p|pi)` = NA, `Var(p|pi)` = NA)
  for (i in 1:nrow(perf)) {
    pval_i              <- dplyr::filter(pval, method == perf$method[i])
    pmf_i               <- dplyr::filter(pmf, pi == perf$pi[i])
    perf$`E(p|pi)`[i]   <- sum(pmf_i$`f(s,m|pi)`*pval_i$`p(s,m)`)
    perf$`Var(p|pi)`[i] <- sum(pmf_i$`f(s,m|pi)`*(pval_i$`p(s,m)`)^2) -
                             perf$`E(p|pi)`[i]^2
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " pi-method combinations evaluated...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(pval = pval, perf = perf, des = des, pi = pi,
                        summary = summary)
  class(output) <- "sa_pval_fixed"
  return(output)
}
