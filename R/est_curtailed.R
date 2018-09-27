#' Determine point estimates in a curtailed group sequential single-arm trial
#' design for a single binary endpoint
#'
#' Determines possible point estimates at the end of a curtailed group
#' sequential single-arm trial for a single binary endpoint, as determined using
#' \code{des_curtailed()}. Support is available to compute point estimates using
#' the naive (\code{"naive"}), bias-adjusted (\code{"bias_adj"}),
#' bias-subtracted (\code{"bias_sub"}), conditional (\code{"conditional"}),
#' median unbiased (\code{"mue"}), and UMVUE (\code{"umvue"}) approaches.
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
#' @param des An object of class \code{"sa_des_curtailed"}, as returned by
#' \code{des_curtailed()}.
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
#' conditional (\code{"conditional"}), median unbiased (\code{"mue"}), and UMVUE
#' (\code{"umvue"}) approaches. Defaults to all available methods.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_est_curtailed"} containing the following
#' elements
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
#' # Find the optimal non-stochastically curtailed two-stage design for the
#' default parameters
#' des <- des_curtailed()
#' # Determine the performance of all supported point estimation procedures for
#' # a range of possible response probabilities
#' est <- est_gs(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_curtailed}}, \code{\link{opchar_curtailed}},
#' \code{\link{pval_curtailed}}, \code{\link{ci_curtailed}}, and their
#' associated \code{plot} family of functions.
#' @export
est_curtailed <- function(des, k, pi,
                          method = c("bias_adj", "bias_sub", "conditional",
                                     "naive", "mue", "umvue"), summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_curtailed(des, "des")
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
  check_belong(method, "method", c("naive", "umvue", "mue", "bias_adj",
                                   "bias_sub", "conditional"), "any")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Point estimation for curtailed group-sequential single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k \u2208 {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("You have chosen to use the following methods to construct point estimates\n")
    if ("naive" %in% method) {
      message("  \u2022 Naive.")
    }
    if ("bias_adj" %in% method) {
      message("  \u2022 Bias-adjusted.")
    }
    if ("bias_sub" %in% method) {
      message("  \u2022 Bias-subtracted.")
    }
    if ("conditional" %in% method) {
      message("  \u2022 Conditional.")
    }
    if ("mue" %in% method) {
      message("  \u2022 MUE.")
    }
    if ("umvue" %in% method) {
      message("  \u2022 UMVUE.")
    }
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  J          <- des$des$J; J_curt <- des$des$J_curt; a_curt <- des$des$a_curt;
  r_curt     <- des$des$r_curt; n <- des$des$n; n_curt <- des$des$n_curt
  N          <- c(0, cumsum(n))
  k_curt     <- NULL
  for (i in 1:length(k)) {
    k_curt   <- c(k_curt, (N[k[i]] + 1):N[k[i] + 1])
  }
  terminal   <- terminal_states_gs(J_curt, a_curt, r_curt, n_curt, k_curt)
  est        <- tibble::tibble(s = rep(terminal$s, each = length(method)),
                               m = rep(terminal$m, each = length(method)),
                               k = factor(rep(terminal$k,
                                              each = length(method)), k_curt),
                               method = factor(rep(method, nrow(terminal)),
                                               method), `hat(pi)(s,m)` = NA)
  for (i in 1:length(method)) {
    range                     <- which(est$method == method[i])
    est$`hat(pi)(s,m)`[range] <-
      switch(as.character(est$method[i]),
             naive       = est$s[range]/est$m[range],
             umvue       = est_gs_umvue(est$s[range], est$m[range],
                                        as.numeric(as.character(est$k[range])),
                                        a_curt, r_curt, n_curt),
             mue         = est_gs_mue(est$s[range], est$m[range],
                                      as.numeric(as.character(est$k[range])),
                                      J_curt, a_curt, r_curt, n_curt),
             bias_adj    = est_gs_bias_adj(est$s[range], est$m[range], J_curt,
                                           a_curt, r_curt, n_curt),
             bias_sub    = est_gs_bias_sub(est$s[range], est$m[range], J_curt,
                                           a_curt, r_curt, n_curt),
             conditional =
               est_gs_cond_mle(est$s[range], est$m[range],
                               as.numeric(as.character(est$k[range])), a_curt,
                               r_curt, n_curt))
    if (summary) {
      message("...point estimates for method \"", method[i], "\" determined...")
    }
  }
  terminal$k <- as.numeric(as.character(terminal$k))
  og_k       <- numeric(nrow(terminal))
  for (i in 1:length(og_k)) {
    og_k[i]  <- which(N[1:J] < terminal$k[i] & N[2:(J + 1)] >= terminal$k[i])
  }
  est$k      <- factor(rep(og_k, each = length(method)), unique(og_k))
  len_pi     <- length(pi)
  pmf        <- pmf_gs(pi, J_curt, a_curt, r_curt, n_curt, k_curt)
  perf       <- tibble::tibble(pi = rep(pi, length(method)),
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
  #Bias and rmse should be in a mutate afterward

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(est = est, perf = perf, des = des, k = k, pi = pi,
                        method = method, summary = summary)
  class(output) <- "sa_est_curtailed"
  return(output)
}
