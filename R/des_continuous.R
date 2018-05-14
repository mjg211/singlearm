des_continuous <- function(J = 2, pi0 = 0.1, pi1 = 0.3, alpha = 0.05,
                           beta = 0.2, d = 5, sigma = 20, efficacy = F,
                           futility = T, optimality = "null-ess", equal_n = F,
                           summary = F, parallel = F, cpus = 1, popSize = 100,
                           maxiter = 1000, run = 100, seed = Sys.time()) {

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(0, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), "1")
  check_real_range_strict(beta, "beta", c(0, 1), "1")
  check_real_range_strict(d, "d", c(0, Inf), "1")
  check_real_range_strict(sigma, "sigma", c(0, Inf), "1")
  check_stopping(futility, efficacy)
  check_belong(optimality, "optimality", c("null-ess", "alt-ess",
                                           "delta-minimax"), "1")
  check_logical(equal_n, "equal_n")
  check_logical(summary, "summary")
  check_logical(parallel, "parallel")
  check_integer_range(cpus, "cpus", c(0, Inf))
  check_integer_range(popSize, "popSize", c(0, Inf))
  check_integer_range(maxiter, "maxiter", c(0, Inf))
  check_integer_range(run, "run", c(0, Inf))

  ##### Print Summary ##########################################################

  ##### Main Computations ######################################################

  delta1     <- d - sigma*stats::qnorm(1 - pi1)
  delta0     <- d - sigma*stats::qnorm(1 - pi0)
  penalty    <- ((stats::qnorm(1 - alpha) +
                    stats::qnorm(1 - beta))*sigma/(delta1 - delta0))^2
  if (parallel) {
    parallel <- cpus
  }
  if (all(efficacy, futility)) {
    fitness  <- function(...){
      -obj_fn_gs_continuous_ef(...)
    }
    if (!summary) {
      sink("NULL")
    }
    opt_des <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                      delta = delta1 - delta0, alpha = alpha, beta = beta,
                      sigma = sigma, optimality = optimality, equal_n = equal_n,
                      penalty = penalty, min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                                                 rep(-20, J), rep(0, J - 1)),
                      max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                              rep(20, J), rep(20, J - 1)), popSize = popSize,
                      maxiter = maxiter, run = run, parallel = parallel,
                      seed = seed)@solution
    if (!summary) {
      sink()
    }
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      a     <- opt_des[2:(J + 1)]
      r     <- c(opt_des[2:J] + opt_des[(J + 2):(2*J)], opt_des[J + 1])
    } else {
      n     <- opt_des[1:J]
      a     <- opt_des[(J + 1):(2*J)]
      r     <- c(opt_des[(J + 1):(2*J - 1)] + opt_des[(2*J + 1):(3*J - 1)],
                 opt_des[2*J])
    }
    N       <- cumsum(n)
    I       <- N/(sigma^2)
    Sigma   <- cov_gs_continuous(J, I)
    perf_H0 <- opchar_ef(0, J, r, a, I, N, Sigma)
    perf_H1 <- opchar_ef(delta1 - delta0, J, r, a, I, N, Sigma)
  } else if (efficacy) {
    fitness <- function(...){
      -obj_fn_gs_continuous_e(...)
    }
    if (!summary) {
      sink("NULL")
    }
    opt_des <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                      delta = delta1 - delta0, alpha = alpha, beta = beta,
                      sigma = sigma, optimality = optimality, equal_n = equal_n,
                      penalty = penalty, min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                                                 rep(-20, J)),
                      max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                              rep(20, J)), popSize = popSize, maxiter = maxiter,
                      run = run, parallel = parallel, seed = seed)@solution
    if (!summary) {
      sink()
    }
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      a     <- c(rep(-Inf, J - 1), opt_des[J + 1])
      r     <- opt_des[2:(J + 1)]
    } else {
      n     <- opt_des[1:J]
      a     <- c(rep(-Inf, J - 1), opt_des[2*J])
      r     <- opt_des[(J + 1):(2*J)]
    }
    N       <- cumsum(n)
    I       <- N/(sigma^2)
    Sigma   <- cov_gs_continuous(J, I)
    perf_H0 <- opchar_e(0, J, r, I, N, Sigma)
    perf_H1 <- opchar_e(delta1 - delta0, J, r, I, N, Sigma)
  } else {
    fitness <- function(...){
      -obj_fn_gs_continuous_f(...)
    }
    if (!summary) {
      sink("NULL")
    }
    opt_des <- GA::ga(type = "real-valued", fitness = fitness, J = J,
                      delta = delta1 - delta0, alpha = alpha, beta = beta,
                      sigma = sigma, optimality = optimality, equal_n = equal_n,
                      penalty = penalty, min = c(rep(0, 1 + (J - 1)*(!equal_n)),
                                                 rep(-20, J)),
                      max = c(rep(2*penalty, 1 + (J - 1)*(!equal_n)),
                              rep(20, J)), popSize = popSize, maxiter = maxiter,
                      run = run, parallel = parallel, seed = seed)@solution
    if (!summary) {
      sink()
    }
    if (equal_n) {
      n     <- rep(opt_des[1], J)
      a     <- opt_des[2:(J + 1)]
      r     <- c(rep(Inf, J - 1), opt_des[J + 1])
    } else {
      n     <- opt_des[1:J]
      a     <- opt_des[(J + 1):(2*J)]
      r     <- c(rep(Inf, J - 1), opt_des[2*J])
    }
    N       <- cumsum(n)
    I       <- N/(sigma^2)
    Sigma   <- cov_gs_continuous(J, I)
    perf_H0 <- opchar_f(0, J, a, I, N, Sigma)
    perf_H1 <- opchar_f(delta1 - delta0, J, a, I, N, Sigma)
  }

  des <- list(J = J, a = a, r = r, n = n, perf_H0 = perf_H0, perf_H1 = perf_H1)

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, pi0 = pi0, pi1 = pi1,
                        alpha = alpha, beta = beta, d = d, sigma = sigma,
                        efficacy = efficacy, futility = futility,
                        optimality = optimality, equal_n = equal_n,
                        summary = summary, parallel = parallel, cpus = cpus,
                        popSize = popSize, maxiter = maxiter, run = run,
                        seed = seed)
  class(output) <- "sa_des_continuous"
  return(output)

}
