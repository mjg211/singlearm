#' Design a Bayesian-frequentist single-arm trial for a single binary endpoint
#'
#' Determines optimised single- and two-stage Bayesian-frequentst single-arm
#' clinical trial designs for a single binary primary endpoint, using exact
#' calculations.
#'
#' Designs controlling Bayesian, frequentist, or Bayesian and frequentist
#' operating characteristics can be determining, which optimise either the
#' Bayesian expected sample size or the maximal sample size.
#'
#' @param J The maximal number of stages to allow.
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desiable) response probability used in the definition of the
#' alternative hypothesis.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param mu The first shape parameter of the Beta distribution.
#' @param nu The second shape parameter of the Beta distribution.
#' @param Nmin The minimal total sample size to allow in considered
#' designs.
#' @param Nmax The maximal total sample size to allow in considered designs.
#' @param optimality Choice of optimal design criteria. Must be one of
#' \code{"ess"} or \code{"minimax"}.
#' @param control Error-rates to control. Should be a vector containing elements
#' chosen from "frequentist" and "bayesian".
#' @param equal_n A logical variable indicating that the sample size of each
#' stage should be equal.
#' @param PL Predictive probability value used in determining when to stop the
#' trial early for futility.
#' @param PU Predictive probability value used in determining when to stop the
#' trial early for futility.
#' @param PT Terminal predictie probability value used in determining when the
#' trial is a success.
#' @param summary A logical variable indicating a summary of the function's
#' progress should be printed to the console.
#' @return A list of class \code{"sa_des_bayesfreq"} containing the following
#' elements
#' \itemize{
#' \item A list in the slot \code{$des} containing details of the
#' identified optimal design.
#' \item A tibble in the slot \code{$feasible}, consisting of the
#' identified designs which met the required operating characteristics.
#' \item Each of the input variables as specified.
#' }
#' @examples
#' # The ESS-optimal design for the default parameters
#' ess_optimal <- des_bayesfreq()
#' # The corresponding minimax design
#' minimax     <- des_bayesfreq(optimality = "minimax")
#' @seealso \code{\link{opchar_bayesfreq}}, and their associated \code{plot}
#' family of functions.
#' @export
des_bayesfreq <- function(J = 2, pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                          mu = 0.1, nu = 0.9, Nmin = 1, Nmax = 30,
                          optimality = "ess", control = c("frequentist",
                                                          "bayesian"),
                          equal_n = F, PL = 0.5, PU = 0.9, PT = 0.95,
                          summary = F) {

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(0, 3))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_real_range_strict(mu, "mu", c(0, Inf), 1)
  check_real_range_strict(nu, "nu", c(0, Inf), 1)
  check_integer_pair_range(Nmin, Nmax, "Nmin", "Nmax", c(0, Inf))
  check_belong(optimality, "optimality", c("minimax", "ess"), "1")
  check_belong(control, "control", c("frequentist", "bayesian"), "any")
  check_logical(equal_n, "equal_n")
  check_real_range_strict(PL, "PL", c(-Inf, Inf), "1")
  check_real_range_strict(PU, "PU", c(-Inf, Inf), "1")
  check_real_range_strict(PT, "PT", c(-Inf, Inf), "1")
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Design of ", J, "-stage Bayesian-frequentist group sequential single-arm trials with a single binary endpoint")
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
    if (equal_n) {
      Sys.sleep(2)
      if (J == 2) {
        message("You have chosen to restrict the allowed values of the n\u2c7c, j = 1,2, such that\n")
        message("  \u2022 n\u2081 = n\u2082.\n")
      }
    }
    Sys.sleep(2)
    message("You have chosen to restrict the allowed values in a and r such that\n")
    message("  \u2022 a", sub_num(J), " + 1 = r", sub_num(J), ".")
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  prob_s1_n1   <- matrix(0, nrow = Nmax + 1, ncol = Nmax)
  Beta         <- beta(mu, nu)
  for (n1 in 1:Nmax) {
    prob_s1_n1[1:(n1 + 1), n1] <- choose(n1, 0:n1)*beta(mu + 0:n1,
                                                        nu + n1 - 0:n1)/Beta
  }
  dbinomial_pi0 <- matrix(0, Nmax + 1, Nmax)
  dbinomial_pi1 <- matrix(0, Nmax + 1, Nmax)
  for (n in 1:Nmax) {
    dbinomial_pi0[1:(n + 1), n] <- stats::dbinom(0:n, n, pi0)
    dbinomial_pi1[1:(n + 1), n] <- stats::dbinom(0:n, n, pi1)
  }
  if (J == 2) {
    prob_s2_s1n1n2 <- array(0, c(Nmax + 1, Nmax + 1, Nmax, Nmax))
    for (n1 in 1:Nmax) {
      for (s1 in 0:n1) {
        for (n2 in 1:Nmax) {
          prob_s2_s1n1n2[1:(n2 + 1), s1 + 1, n1, n2] <- choose(n2, 0:n2)*
            beta(mu + s1 + 0:n2,
                 nu + n1 - s1 + n2 - 0:n2)/
            beta(mu + s1, nu + n1 - s1)
        }
      }
    }
    PP_TS_s1n1n2 <- array(0, c(Nmax + 1, Nmax, Nmax))
    for (n1 in 1:Nmax) {
      for (n2 in 1:Nmax) {
        for (s1 in 0:n1) {
          cond_s2 <- prob_s2_s1n1n2[1:(n2 + 1), s1 + 1, n1, n2]
          prob    <- stats::pbeta(pi0, mu + s1 + 0:n2, nu + n1 + n2 - s1 - 0:n2,
                           lower.tail = F)
          PP_TS_s1n1n2[s1 + 1, n1, n2] <- sum(cond_s2[which(prob > PT)])
        }
      }
    }
  }
  feasible <- matrix(0, nrow = 10^7, ncol = 7 + 10*(J == 2))
  counter  <- 1
  if (J == 1) {
    AB_pi1 <- 0
    RB_pi0 <- 0
    AF_pi1 <- 0
    RF_pi0 <- 0
    for (n in Nmin:Nmax) {
      for (a in 0:(n - 1)) {
        r      <- a + 1
        RB_pi0 <- sum(prob_s1_n1[(r + 1):(n + 1), n]*stats::pbeta(pi0, mu + r:n,
                                                           nu + n - r:n))
        RF_pi0 <- sum(dbinomial_pi0[(r + 1):(n + 1), n])
        AB_pi1 <- sum(prob_s1_n1[1:(a + 1), n]*stats::pbeta(pi1, mu + 0:a,
                                                     nu + n - 0:a,
                                                     lower.tail = F))
        AF_pi1 <- sum(dbinomial_pi1[1:(a + 1), n])
        if (length(control) == 1) {
          if (control == "frequentist") {
            if (all(RF_pi0 <= alpha, AF_pi1 <= beta)) {
              feasible[counter, ] <- c(n, a, r, sum(RB_pi0), 1 - sum(AB_pi1),
                                       sum(RF_pi0), 1 - sum(AF_pi1))
              counter             <- counter + 1
            }
          } else {
            if (all(RB_pi0 <= alpha, AB_pi1 <= beta)) {
              feasible[counter, ] <- c(n, a, r, sum(RB_pi0), 1 - sum(AB_pi1),
                                       sum(RF_pi0), 1 - sum(AF_pi1))
              counter             <- counter + 1
            }
          }
        } else {
          if (all(RB_pi0 <= alpha, RF_pi0 <= alpha, AB_pi1 <= beta,
                  AF_pi1 <= beta)) {
            feasible[counter, ] <- c(n, a, r, sum(RB_pi0), 1 - sum(AB_pi1),
                                     sum(RF_pi0), 1 - sum(AF_pi1))
            counter             <- counter + 1
          }
        }
      }
      if (all(summary, n%%10 == 0)) {
        message("...completed evaluation of designs with n = ", n, "...")
      }
    }
  } else {
    for (n1 in switch(equal_n + 1, 1:(Nmax - 1),
                      ceiling(Nmin/2):floor(Nmax/2))) {
      for (n2 in switch(equal_n + 1, max(1, Nmin - n1):(Nmax - n1), n1)) {
        n      <- n1 + n2
        AB_pi1 <- numeric(2)
        RB_pi0 <- numeric(2)
        AF_pi1 <- numeric(2)
        RF_pi0 <- numeric(2)
        PP_TS  <- PP_TS_s1n1n2[1:(n1 + 1), n1, n2]
        a1     <- suppressWarnings(max(which(PP_TS <= PL)) - 1)
        r1     <- sign(suppressWarnings(min(which(PP_TS >= PU)) - 1))*
                    suppressWarnings(min(which(PP_TS >= PU)) - 1)
        if (r1 < Inf) {
          RB_pi0[1]  <- sum(prob_s1_n1[(r1 + 1):(n1 + 1), n1]*
                              stats::pbeta(pi0, mu + r1:n1, nu + n1 - r1:n1))/
                          sum(prob_s1_n1[(r1 + 1):(n1 + 1), n1])
          RF_pi0[1]  <- sum(dbinomial_pi0[(r1 + 1):(n1 + 1), n1])
        }
        if (a1 > -Inf) {
          AB_pi1[1]  <- sum(prob_s1_n1[1:(a1 + 1), n1]*stats::pbeta(pi1, mu + 0:a1,
                                                             nu + n1 - 0:a1,
                                                             lower.tail = F))/
                          sum(prob_s1_n1[1:(a1 + 1), n1])
          AF_pi1[1]  <- sum(dbinomial_pi1[1:(a1 + 1), n1])
        }
        if (all(RB_pi0[1] <= alpha, RF_pi0[1] <= alpha, AB_pi1[1] <= beta,
                AF_pi1[1] <= beta)) {
          PETB       <- 0
          PETF_pi0   <- 0
          PETF_pi1   <- 0
          if (a1 > -Inf) {
            PETB     <- sum(prob_s1_n1[1:(a1 + 1), n1])
            PETF_pi0 <- sum(dbinomial_pi0[1:(a1 + 1), n1])
            PETF_pi1 <- sum(dbinomial_pi1[1:(a1 + 1), n1])
          }
          if (r1 < Inf) {
            PETB     <- PETB + sum(prob_s1_n1[(r1 + 1):(n1 + 1), n1])
            PETF_pi0 <- PETF_pi0 + sum(dbinomial_pi0[(r1 + 1):(n1 + 1), n1])
            PETF_pi1 <- PETF_pi1 + sum(dbinomial_pi1[(r1 + 1):(n1 + 1), n1])
          }
          ESSB       <- n1 + (1 - PETB)*n2
          ESSF_pi0   <- n1 + (1 - PETF_pi0)*n2
          ESSF_pi1   <- n1 + (1 - PETF_pi1)*n2
          if (a1 < r1 - 1) {
            for (a2 in max(0, a1 + 1):min(r1 + n2 - 2, n1 + n2 - 1)) {
              r2    <- a2 + 1
              numer <- 0
              denom <- 0
              freq  <- 0
              for (s1 in max(0, a1 + 1):min(r1 - 1, n1)) {
                if (r2 - s1 <= n2) {
                  for (s2 in max(r2 - s1, 0):n2) {
                    s     <- s1 + s2
                    numer <- numer + sum(choose(n1, s1)*choose(n2, s2)*
                                           beta(mu + s, nu + n - s)*
                                           stats::pbeta(pi0, mu + s, nu + n - s)/Beta)
                    denom <- denom + sum(choose(n1, s1)*choose(n2, s2)*
                                           beta(mu + s, nu + n - s)/Beta)
                    freq  <- freq + dbinomial_pi0[s1 + 1, n1]*
                      dbinomial_pi0[s2 + 1, n2]
                  }
                }
              }
              RB_pi0[2]   <- numer/denom
              RF_pi0[2]   <- freq
              numer       <- 0
              denom       <- 0
              freq        <- 0
              for (s1 in max(0, a1 + 1):min(r1 - 1, n1)) {
                if (s1 <= a2) {
                  for (s2 in 0:(a2 - s1)) {
                    s     <- s1 + s2
                    numer <- numer + sum(choose(n1, s1)*choose(n2, s2)*
                                           beta(mu + s, nu + n - s)*
                                           stats::pbeta(pi1, mu + s, nu + n - s,
                                                 lower.tail = F)/Beta)
                    denom <- denom + sum(choose(n1, s1)*choose(n2, s2)*
                                           beta(mu + s, nu + n - s)/Beta)
                    freq  <- freq + dbinomial_pi1[s1 + 1, n1]*
                      dbinomial_pi1[s2 + 1, n2]
                  }
                }
              }
              AB_pi1[2] <- numer/denom
              AF_pi1[2] <- freq
              if (length(control) == 1) {
                if (control == "frequentist") {
                  if (all(sum(RF_pi0) <= alpha, sum(AF_pi1) <= beta)) {
                    feasible[counter, ] <- c(n1, n2, a1, a2, r1, r2, sum(RB_pi0),
                                             1 - sum(AB_pi1), sum(RF_pi0),
                                             1 - sum(AF_pi1), ESSB, ESSF_pi0,
                                             ESSF_pi1, PETB, PETF_pi0, PETF_pi1,
                                             n1 + n2)
                    counter             <- counter + 1
                  }
                } else {
                  if (all(sum(RB_pi0) <= alpha, sum(AB_pi1) <= beta)) {
                    feasible[counter, ] <- c(n1, n2, a1, a2, r1, r2, sum(RB_pi0),
                                             1 - sum(AB_pi1), sum(RF_pi0),
                                             1 - sum(AF_pi1), ESSB, ESSF_pi0,
                                             ESSF_pi1, PETB, PETF_pi0, PETF_pi1,
                                             n1 + n2)
                    counter             <- counter + 1
                  }
                }
              } else {
                if (all(sum(RB_pi0) <= alpha, sum(RF_pi0) <= alpha,
                        sum(AB_pi1) <= beta, sum(AF_pi1) <= beta)) {
                  feasible[counter, ] <- c(n1, n2, a1, a2, r1, r2, sum(RB_pi0),
                                           1 - sum(AB_pi1), sum(RF_pi0),
                                           1 - sum(AF_pi1), ESSB, ESSF_pi0,
                                           ESSF_pi1, PETB, PETF_pi0, PETF_pi1,
                                           n1 + n2)
                  counter             <- counter + 1
                }
              }
            }
          }
        }
      }
      if (all(summary, n1%%10 == 0)) {
        message("...completed evaluation of designs with n1 = ", n1, "...")
      }
    }
  }
  if (counter > 1) {
    feasible             <- tibble::as_tibble(feasible)
    feasible             <- feasible[1:(counter - 1), ]
    if (J == 1) {
      colnames(feasible) <- c("n", "a", "r", "PredPB(pi0)", "PB(pi1)", "PF(pi0)",
                              "PF(pi1)")
      feasible           <- dplyr::arrange(feasible, n, dplyr::desc(`PB(pi1)`))
      des                <- list(J = J, n = as.numeric(feasible$n[1]),
                                 a = as.numeric(feasible$a[1]),
                                 r = as.numeric(feasible$r[1]), pi0 = pi0,
                                 pi1 = pi1, alpha = alpha, beta = beta, mu = mu,
                                 nu = nu, opchar = feasible[1, ])
      feasible[, 1:3]    <- dplyr::mutate_if(feasible[, 1:3], is.double,
                                             as.integer)
    } else {
      colnames(feasible) <- c("n1", "n2", "a1", "a2", "r1", "r2", "margPB(pi0)",
                              "margPB(pi1)", "PF(pi0)", "PF(pi1)", "ESSB",
                              "ESSF(pi0)", "ESSF(pi1)", "PETB", "PETF(pi0)",
                              "PETF(pi1)", "max(N)")
      if (optimality == "ess") {
        feasible         <- dplyr::arrange(feasible, ESSB, `max(N)`,
                                           dplyr::desc(`margPB(pi1)`))
      } else {
        feasible         <- dplyr::arrange(feasible, `max(N)`, ESSB,
                                           dplyr::desc(`margPB(pi1)`))
      }
      des                <- list(J = J, n = as.numeric(feasible[1, 1:2]),
                                 a = as.numeric(feasible[1, 3:4]),
                                 r = as.numeric(feasible[1, 5:6]),
                                 pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                                 mu = mu, nu = nu,
                                 opchar = feasible[1, ])
      feasible[, c(1:6, 17)] <- dplyr::mutate_if(feasible[, c(1:6, 17)],
                                                 is.double, as.integer)
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = feasible, J = J, pi0 = pi0,
                        pi1 = pi1, alpha = alpha, beta = beta, mu = mu, nu = nu,
                        Nmin = Nmin, Nmax = Nmax, optimality = optimality,
                        equal_n = equal_n, PL = PL, PU = PU, PT = PT,
                        summary = summary)
  class(output) <- "sa_des_bayesfreq"
  return(output)
}
