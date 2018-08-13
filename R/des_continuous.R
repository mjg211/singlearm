des_continuous <- function(J = 2, delta0 = 10, delta1 = 0, sigma = 10,
                           alpha = 0.05, beta = 0.2, efficacy = F, futility = T,
                           optimality = "null_ess", equal_n = F, summary = F,
                           parallel = F, popSize = 100, maxiter = 1000,
                           run = 100, seed = Sys.time()) {

  # Need a sigma_known option - that only works with J=1,2 - not sure I want to
  # use James approach - think treat all of them as continuous

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(0, Inf))
  check_real_pair_range_strict(delta0, delta1, "delta0", "delta1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), "1")
  check_real_range_strict(beta, "beta", c(0, 1), "1")
  check_stopping(futility, efficacy)
  check_belong(optimality, "optimality", c("null_ess", "alt_ess",
                                           "delta_minimax"), "1")
  if (all(!futility | !efficacy, optimality == "delta_minimax")) {
    stop("The delta_minimax design is only relevant when efficacy == futility == TRUE")
  }
  check_logical(equal_n, "equal_n")
  check_logical(summary, "summary")
  # Need a check on parallel
  check_integer_range(popSize, "popSize", c(0, Inf))
  check_integer_range(maxiter, "maxiter", c(0, Inf))
  check_integer_range(run, "run", c(0, Inf))

  ##### Print Summary ##########################################################

  ##### Main Computations ######################################################

  if (J == 1) {

  } else {

  }

  penalty    <- ((stats::qnorm(1 - alpha) +
                    stats::qnorm(1 - beta))*sigma/(delta1 - delta0))^2
  if (all(efficacy, futility)) {
    fitness  <- function(...){
      -obj_fn_gs_continuous_ef(...)
    }
    sink("NULL")
    opt_des  <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                       delta0 = delta0, delta1 = delta1, sigma = sigma,
                       alpha = alpha, beta = beta, optimality = optimality,
                       equal_n = equal_n, penalty = penalty,
                       min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                               rep(-20, J), rep(0, J - 1)),
                       max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                               rep(20, J), rep(20, J - 1)), popSize = popSize,
                       maxiter = maxiter, run = run, parallel = parallel,
                       seed = seed)@solution
    sink()
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      f     <- opt_des[2:(J + 1)]
      e     <- c(opt_des[2:J] + opt_des[(J + 2):(2*J)], opt_des[J + 1])
    } else {
      n     <- opt_des[1:J]
      f     <- opt_des[(J + 1):(2*J)]
      e     <- c(opt_des[(J + 1):(2*J - 1)] + opt_des[(2*J + 1):(3*J - 1)],
                 opt_des[2*J])
    }
    N         <- cumsum(n)
    I         <- N/sigma^2
    Sigma     <- cov_gs_continuous(J, I)
    opchar_H0 <- opchar_ef(delta0, J, e, f, I, N, Sigma)
    opchar_H1 <- opchar_ef(delta1, J, e, f, I, N, Sigma)
  } else if (efficacy) {
    fitness <- function(...){
      -obj_fn_gs_continuous_e(...)
    }
    sink("NULL")
    opt_des <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                      delta0 = delta0, delta1 = delta1, sigma = sigma,
                      optimality = optimality, equal_n = equal_n,
                      penalty = penalty, min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                                                 rep(-20, J)),
                      max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                              rep(20, J)), popSize = popSize, maxiter = maxiter,
                      run = run, parallel = parallel, seed = seed)@solution
    sink()
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      f     <- c(rep(-Inf, J - 1), opt_des[J + 1])
      e     <- opt_des[2:(J + 1)]
    } else {
      n     <- opt_des[1:J]
      f     <- c(rep(-Inf, J - 1), opt_des[2*J])
      e     <- opt_des[(J + 1):(2*J)]
    }
    N       <- cumsum(n)
    I       <- N/sigma^2
    Sigma   <- cov_gs_continuous(J, I)
    opchar_H0 <- opchar_e(delta0, J, e, I, N, Sigma)
    opchar_H1 <- opchar_e(delta1, J, e, I, N, Sigma)
  } else {
    fitness <- function(...){
      -obj_fn_gs_continuous_f(...)
    }
    sink("NULL")
    opt_des <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                      delta0 = delta0, delta1 = delta1, sigma = sigma,
                      optimality = optimality, equal_n = equal_n,
                      penalty = penalty, min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                                                 rep(-20, J)),
                      max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                              rep(20, J)), popSize = popSize, maxiter = maxiter,
                      run = run, parallel = parallel, seed = seed)@solution
    sink()
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      f     <- opt_des[2:(J + 1)]
      e     <- c(rep(Inf, J - 1), opt_des[J + 1])
    } else {
      n     <- opt_des[1:J]
      f     <- opt_des[(J + 1):(2*J)]
      e     <- c(rep(Inf, J - 1), opt_des[2*J])
    }
    N       <- cumsum(n)
    I       <- N/sigma^2
    Sigma   <- cov_gs_continuous(J, I)
    opchar_H0 <- opchar_f(delta0, J, f, I, N, Sigma)
    opchar_H1 <- opchar_f(delta1, J, f, I, N, Sigma)
  }
  des <- list(J = J, f = f, e = e, n = n, opchar_H0 = opchar_H0,
              opchar_H1 = opchar_H1)

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, delta0 = delta0, delta1 = delta1,
                        alpha = alpha, beta = beta, sigma = sigma,
                        efficacy = efficacy, futility = futility,
                        optimality = optimality, equal_n = equal_n,
                        summary = summary, parallel = parallel,
                        popSize = popSize, maxiter = maxiter, run = run,
                        seed = seed)
  class(output) <- "sa_des_continuous"
  return(output)
}
