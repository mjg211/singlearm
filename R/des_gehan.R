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
des_gehan <- function(pi1 = 0.3, beta1 = 0.1, gamma = 0.05, alpha_pi_hat = 0.25,
                      conservative = F, alpha = 0.05, pi0 = 0.1, find_D = F,
                      summary = F) {

  ##### Input Checking #########################################################

  check_logical(find_D, "find_D")
  if (find_D) {
    check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
    check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  } else {
    check_real_range_strict(pi1, "pi1", c(0, 1), "1")
  }
  check_real_range_strict(alpha_pi_hat, "alpha_pi_hat", c(0, 1), 1)
  check_logical(conservative, "conservative")
  check_real_range_strict(beta1, "beta1", c(0, 1), 1)
  check_real_range_strict(gamma, "gamma", c(0, 1), "1")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################


  ##### Main Computations ######################################################

  n1              <- 1
  rej_error       <- stats::dbinom(0, n1, pi1)
  while (rej_error > beta1) {
    n1            <- n1 + 1
    rej_error     <- stats::dbinom(0, n1, pi1)
  }
  s1              <- 0:n1
  n2              <- numeric(n1 + 1)
  for (s1 in 1:n1) {
    if (!conservative) {
      pi_hat      <- ci_fixed_wald(s1, n1, alpha_pi_hat)[2]
    } else {
      poss_pi_hat <- c(ci_fixed_clopper_pearson(s1, n1, alpha_pi_hat), s1/n1)
      pi_hat      <- poss_pi_hat[which.min(abs(poss_pi_hat - 0.5))]
    }
    while (sqrt(pi_hat*(1 - pi_hat)/(n1 + n2[s1 + 1])) > gamma) {
      n2[s1 + 1]  <- n2[s1 + 1] + 1
    }
  }
  if (find_D) {
    dbinom_pi0          <- stats::dbinom(0:n1, n1, pi0)
    dbinom_pi1          <- stats::dbinom(0:n1, n1, pi1)
    dc_ef               <- dc_pf <- list()
    length_dc_ef        <- numeric(n1 + 1)
    for (s1 in 1:(n1 + 1)) {
      if (n2[s1] == 0) {
        if (s1 < n1 + 1) {
          if (any(n2[(s1 + 1):(n1 + 1)] > 0)) {
            dc_ef[[s1]] <- dc_pf[[s1]] <- 0
          } else {
            dc_ef[[s1]] <- dc_pf[[s1]] <- 1
          }
        } else {
          dc_ef[[s1]]   <- dc_pf[[s1]] <- 1
        }
      } else {
        dc_ef[[s1]]     <- pbinom((n2[s1] - 1):0, n2[s1], pi0, lower.tail = F)
        dc_pf[[s1]]     <- pbinom((n2[s1] - 1):0, n2[s1], pi1, lower.tail = F)
        dc_pf[[s1]]     <- dc_pf[[s1]][which(dc_ef[[s1]]*dbinom_pi0[s1] <=
                                               alpha)]
        dc_ef[[s1]]     <- dc_ef[[s1]][which(dc_ef[[s1]]*dbinom_pi0[s1] <=
                                               alpha)]
      }
      length_dc_ef[s1]  <- length(dc_ef[[s1]])
    }
    optimal_D           <- gehan_dc_ef(pi0, pi1, alpha, n1, n2, dc_ef, dc_pf,
                                       length_dc_ef)
    D                   <- optimal_D$D
    a1                  <- optimal_D$a1
    r1                  <- optimal_D$r1
    a2                  <- optimal_D$a2
    r2                  <- optimal_D$r2
    opchar              <- int_opchar_adaptive(c(pi0, pi1), a1, r1, a2, r2, n1,
                                               n2, 1:2)
    des                 <- list(J = 2, n1 = n1, n2 = n2, a1 = a1, r1 = r1,
                                a2 = a2, r2 = r2, D = D, pi0 = pi0, pi1 = pi1,
                                beta1 = beta1, gamma = gamma,
                                alpha_pi_hat = alpha_pi_hat, alpha = alpha,
                                opchar = opchar)
  } else {
    a1     <- utils::tail(which(cumsum(n2) == 0), n = 1) - 1
    if (n2[n1 + 1] == 0) {
      r1   <- max(which(n2 > 0))
    } else {
      r1   <- Inf
    }
    opchar <- int_opchar_gehan(pi1, a1, r1, n1, n2, 1:2)
    des    <- list(J = 2, n1 = n1, n2 = n2, a1 = a1, r1 = r1, pi1 = pi1,
                   beta1 = beta1, gamma = gamma, alpha_pi_hat = alpha_pi_hat,
                   opchar = opchar)
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, pi1 = pi1, beta1 = beta1, gamma = gamma,
                        alpha_pi_hat = alpha_pi_hat, alpha = alpha, pi0 = pi0,
                        find_D = find_D, summary = summary)
  class(output) <- "sa_des_gehan"
  return(output)
}
