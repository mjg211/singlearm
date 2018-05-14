#' Design a two-stage Gehan single-arm trial for a single binary endpoint
#'
#' Determines two-stage Gehan single-arm clinical trial designs for a single
#' binary primary endpoint.
#'
#' \code{des_gehan()} supports the determination of two-stage Gehan
#' single-arm clinical trial designs for a single binary primary endpoint.
#'
#' @param pi1 The (desirable) response probability used in computing the first
#' stage sample size.
#' @param beta1 The desired maximal type-II error-rate for stage one.
#' @param gamma The desired standard error for the estimate of the response
#' probability by the end of stage two.
#' @param alpha The confidence level to use in the formula for computing the
#' second stage sample sizes.
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
#' # The default design
#' gehan <- des_gehan()
#' @export
des_gehan <- function(pi1 = 0.3, beta1 = 0.1, gamma = 0.05, alpha = 1,
                      summary = F) {

  ##### Input Checking #########################################################

  check_real_range_strict(pi1, "pi1", c(0, 1), "1")
  check_real_range_strict(beta1, "beta1", c(0, 1), "1")
  check_real_range_strict(gamma, "gamma", c(0, 1), "1")
  check_real_range(alpha, "alpha", c(0, 1))
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################


  ##### Main Computations ######################################################

  n1          <- 1
  rej_error   <- dbinom(0, n1, pi1)
  while (rej_error > beta1) {
    n1        <- n1 + 1
    rej_error <- dbinom(0, n1, pi1)
  }
  s1          <- 0:n1
  n2          <- numeric(n1 + 1)
  for (s1 in 1:n1) {
    poss_pi_hat  <- c(ci_fixed_wald(s1, n1, alpha), s1/n1)
    pi_hat       <- poss_pi_hat[which.min(abs(poss_pi_hat - 0.5))]
    while (sqrt(pi_hat*(1 - pi_hat)/(n1 + n2[s1 + 1])) > gamma) {
      n2[s1 + 1] <- n2[s1 + 1] + 1
    }
  }
  a1   <- tail(which(cumsum(n2) == 0), n = 1) - 1
  if (n2[n1 + 1] == 0) {
    r1 <- max(which(n2 > 0))
  } else {
    r1 <- Inf
  }
  opchar <- tibble::as_tibble(matrix(c(n1, n2, a1, r1,
                                       unlist(int_opchar_gehan(pi1, a1, r1, n1,
                                                               n2, 1:2))[-1],
                                       n1 + max(n2)), nrow = 1))
  colnames(opchar) <- c("n1", paste("n2(", 0:n1, ")", sep = ""),
                        "a1", "r1", "ESS(pi1)", "VSS(pi1)", "Med(pi1)", "S1(pi1)",
                        "S2(pi1)", "cumS1(pi1)", "cumS2(pi1)", "max(N)")
  des                <- list(J = 2, n1 = n1, n2 = n2, a1 = a1, r1 = r1,
                             pi1 = pi1, beta1 = beta1, gamma = gamma,
                             alpha = alpha, opchar = opchar)

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, pi1 = pi1, beta1 = beta1, gamma = gamma,
                        alpha = alpha, summary = summary)
  class(output) <- "sa_des_gehan"
  return(output)
}
