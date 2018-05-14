#' Determine point estimates in a group sequential single-arm trial design for
#' a single binary endpoint
#'
#' Determines possible point estimates at the end of a
#' group sequential single-arm trial for a single binary endpoint, as determined
#' using \code{des_gs()}. Support is available to compute point estimates using
#' the naive (\code{"naive"}), bias-adjusted (\code{"bias_adj"}),
#' bias-subtracted (\code{"bias_sub"}), conditional (\code{"conditional"}),
#' median unbiased (\code{"mue"}), UMVUE (\code{"umvue"}), and UMVCUE
#' (\code{"umvcue"}) approaches.
#'
#' In addition, the performance of the chosen point estimate procedures
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
#' integers, with elements between one and the maximum number of possible
#' stages. If left unspecified, it will internally default to all possible
#' stages.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the point estimation procedures at. This will internally
#' default to be the \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}
#' and \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from
#' \code{des} if it is left unspecified.
#' @param method A vector of methods to use to construct point estimates.
#' Currently, support is available to use the naive (\code{"naive"}),
#' bias-adjusted (\code{"bias_adj"}), bias-subtracted (\code{"bias_sub"}),
#' conditional (\code{"conditional"}), median unbiased (\code{"mue"}), UMVUE
#' (\code{"umvue"}), and UMVCUE (\code{"umvcue"}) approaches. Defaults to all
#' available methods.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_est_gs"} containing the following elements
#' \itemize{
#' \item A tibble in the slot \code{$est} summarising the possible point
#' estimates at the end of the trial for the supplied design, according to the
#' chosen methods.
#' \item A tibble in the slot \code{$perf} summarising the performance of the
#' chosen point estimation procedures for each specified value of
#' \ifelse{html}{\out{<i>&pi;</i>}}{\deqn{\pi}}.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal two-stage design for the default parameters
#' des <- des_gs()
#' # Determine the performance of all supported point estimation procedures for
#' # a range of possible response probabilities
#' est <- est_gs(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}},
#' \code{\link{pval_gs}}, \code{\link{ci_gs}}, and their associated \code{plot}
#' family of functions.
#' @export
est_gs <- function(des, k, pi, method = c("bias_adj", "bias_sub", "conditional",
                                          "naive", "mue", "umvcue", "umvue"),
                   summary = F) {

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
  check_belong(method, "method", c("naive", "umvue", "umvcue", "mue",
                                   "bias_adj", "bias_sub", "conditional"),
               "any")
  if (all(des$des$J > 2, "umvcue" %in% method)) {
    stop("umvcue is not currently supported for designs with J > 2")
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Point estimation for group-sequential single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k ∈ {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("You have chosen to use the following methods to construct point estimates\n")
    if ("naive" %in% method) {
      message("  • Naive.")
    }
    if ("bias_adj" %in% method) {
      message("  • Bias-adjusted.")
    }
    if ("bias_sub" %in% method) {
      message("  • Bias-subtracted.")
    }
    if ("conditional" %in% method) {
      message("  • Conditional.")
    }
    if ("mue" %in% method) {
      message("  • MUE.")
    }
    if ("umvue" %in% method) {
      message("  • UMVUE.")
    }
    if ("umvcue" %in% method) {
      message("  • UMVCUE.")
    }
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  J        <- des$des$J; a <- des$des$a; r <- des$des$r; n <- des$des$n
  terminal <- terminal_states_gs(J, a, r, n, k)
  est      <- tibble::tibble(s = rep(terminal$s, each = length(method)),
                             m = rep(terminal$m, each = length(method)),
                             k = factor(rep(terminal$k,
                                            each = length(method)), k),
                             method = factor(rep(method, nrow(terminal)),
                                             method), `hat(pi)(s,m)` = NA)
  for (i in 1:nrow(est)) {
    est$`hat(pi)(s,m)`[i] <-
      switch(as.character(est$method[i]),
             naive       = est$s[i]/est$m[i],
             umvue       = est_gs_umvue(est$s[i], est$m[i],
                                        as.numeric(as.character(est$k[i])), a,
                                        r, n),
             umvcue      = est_gs_umvcue(est$s[i], est$m[i],
                                         as.numeric(as.character(est$k[i])), a,
                                         r, n),
             mue         = est_gs_mue(est$s[i], est$m[i],
                                      as.numeric(as.character(est$k[i])), J, a,
                                      r, n),
             bias_adj    = est_gs_bias_adj(est$s[i], est$m[i], J, a, r, n),
             bias_sub    = est_gs_bias_sub(est$s[i], est$m[i], J, a, r, n),
             conditional = est_gs_cond_mle(est$s[i], est$m[i], est$k[i], a, r,
                                           n))
    if (all(summary, i%%100 == 0)) {
      message("... ", i, " point estimates determined...")
    }
  }
  len_pi <- length(pi)
  pmf    <- pmf_gs(pi, J, a, r, n, k)
  perf   <- tibble::tibble(pi = rep(pi, length(method)),
                           method = factor(rep(method, each = len_pi),
                                               levels = method),
                           `E(hat(pi)|pi)` = NA, `Var(hat(pi)|pi)` = NA,
                           `Bias(hat(pi)|pi)` = NA, `RMSE(hat(pi)|pi)` = NA)
  for (i in 1:nrow(perf)) {
    pmf_i                      <- dplyr::filter(pmf, pi == perf$pi[i])
    est_i                      <- dplyr::filter(est, method == perf$method[i])
    perf$`E(hat(pi)|pi)`[i]    <- sum(pmf_i$`f(s,m|pi)`*est_i$`hat(pi)(s,m)`)
    perf$`Var(hat(pi)|pi)`[i]  <- sum(pmf_i$`f(s,m|pi)`*
                                        (est_i$`hat(pi)(s,m)`)^2) -
                                    perf$`E(hat(pi)|pi)`[i]^2
    perf$`Bias(hat(pi)|pi)`[i] <- perf$`E(hat(pi)|pi)`[i] - perf$pi[i]
    perf$`RMSE(hat(pi)|pi)`[i] <- sqrt(perf$`Var(hat(pi)|pi)`[i] +
                                         perf$`Bias(hat(pi)|pi)`[i]^2)
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " pi-method combinations evaluated...")
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(est = est, perf = perf, des = des, k = k, pi = pi,
                        method = method, summary = summary)
  class(output) <- "sa_est_gs"
  return(output)
}
