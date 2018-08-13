# Function to return covariance matrix in a group sequential design with a
# continuous outcome
cov_gs_continuous <- function(J, I) {
  Sigma           <- diag(1, J, J)
  for (i in 2:J) {
    for (j in 1:(i - 1)) {
      Sigma[j, i] <- Sigma[i, j] <- sqrt(I[j]/I[i])
    }
  }
  return(Sigma)
}

# Function to return objective function in a group sequential design with
# efficacy and futility stopping boundaries
obj_fn_gs_continuous_ef <- function(parameters, J, delta0, delta1, sigma, alpha,
                                    beta, optimality, equal_n, penalty) {
  if (equal_n) {
    n <- rep(parameters[1], J)
    f <- parameters[2:(J + 1)]
    e <- c(parameters[2:J] + parameters[(J + 2):(2*J)], parameters[J + 1])
  } else {
    n <- parameters[1:J]
    f <- parameters[(J + 1):(2*J)]
    e <- c(parameters[(J + 1):(2*J - 1)] + parameters[(2*J + 1):(3*J - 1)],
           parameters[2*J])
  }
  N         <- cumsum(n)
  I         <- N/sigma^2
  Sigma     <- cov_gs_continuous(J, I)
  opchar_H0 <- opchar_ef(delta0, J, e, f, I, N, Sigma)
  opchar_H1 <- opchar_ef(delta1, J, e, f, I, N, Sigma)
  O         <- penalty*((alpha < opchar_H0$P)*(opchar_H0$P - alpha)/alpha +
                          (beta < 1 - opchar_H1$P)*
                          ((1 - opchar_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + opchar_H0$ESS)
  } else if (optimality == "alt-ess") {
    return(O + opchar_H1$ESS)
  } else {
    max_ess <- -suppressWarnings(stats::optim(par = 0.5*delta, fn = max_ess_ef,
                                              J = J, e = e, f = f, I = I, N = N,
                                              Sigma = Sigma)$value)
    return(O + max_ess)
  }
}

# Function to return objective function in a group sequential design with
# efficacy stopping boundaries only
obj_fn_gs_continuous_e <- function(parameters, J, delta0, delta1, sigma, alpha,
                                   beta, optimality, equal_n, penalty) {
  if (equal_n) {
    n <- rep(parameters[1], J)
    e <- parameters[2:(J + 1)]
  } else {
    n <- parameters[1:J]
    e <- parameters[(J + 1):(2*J)]
  }
  N         <- cumsum(n)
  I         <- N/sigma^2
  Sigma     <- cov_gs_continuous(J, I)
  opchar_H0 <- opchar_e(delta0, J, e, I, N, Sigma)
  opchar_H1 <- opchar_e(delta1, J, e, I, N, Sigma)
  O         <- penalty*((alpha < opchar_H0$P)*(opchar_H0$P - alpha)/alpha +
                          (beta < 1 - opchar_H1$P)*
                          ((1 - opchar_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + opchar_H0$ESS)
  } else {
    return(O + opchar_H1$ESS)
  }
}

# Function to return objective function in a group sequential design with
# futility stopping boundaries only
obj_fn_gs_continuous_f <- function(parameters, J, delta0, delta1, sigma, alpha,
                                   beta, optimality, equal_n, penalty) {
  if (equal_n) {
    n <- rep(parameters[1], J)
    f <- parameters[2:(J + 1)]
  } else {
    n <- parameters[1:J]
    f <- parameters[(J + 1):(2*J)]
  }
  N         <- cumsum(n)
  I         <- N/sigma^2
  Sigma     <- cov_gs_continuous(J, I)
  opchar_H0 <- opchar_f(delta0, J, f, I, N, Sigma)
  opchar_H1 <- opchar_f(delta1, J, f, I, N, Sigma)
  O         <- penalty*((alpha < opchar_H0$P)*(opchar_H0$P - alpha)/alpha +
                          (beta < 1 - opchar_H1$P)*
                          ((1 - opchar_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + opchar_H0$ESS)
  } else {
    return(O + opchar_H1$ESS)
  }
}

# Function used to find max_delta{ESS(delta)} for a group sequential design with
# efficacy and futility boundaries
max_ess_ef <- function(delta, J, e, f, I, N, Sigma) {
  return(-opchar_ef(delta, J, e, f, I, N, Sigma)[2])
}

# Function for determining operating characteristics of a group sequential
# design with efficacy and futility boundaries
opchar_ef <- function(delta, J, e, f, I, N, Sigma) {
  Fu <- E <- numeric(J)
  E[1]    <- mvtnorm::pmvnorm(e[1], Inf, delta*sqrt(I[1]), sigma = 1)[1]
  Fu[1]   <- mvtnorm::pmvnorm(-Inf, f[1], delta*sqrt(I[1]), sigma = 1)[1]
  for (j in 2:J) {
    E[j]  <- mvtnorm::pmvnorm(c(f[1:(j - 1)], e[j]), c(e[1:(j - 1)], Inf),
                              delta*sqrt(I[1:j]), sigma = Sigma[1:j, 1:j])[1]
    Fu[j] <- mvtnorm::pmvnorm(c(f[1:(j - 1)], -Inf), c(e[1:(j - 1)], f[j]),
                              delta*sqrt(I[1:j]), sigma = Sigma[1:j, 1:j])[1]
  }
  cum_S   <- cumsum(S <- Fu + E)
  P       <- sum(E)
  ESS     <- sum(N*S)
  Med     <- ifelse(any(cum_S == 0.5),
                    0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1]),
                    N[which(cum_S > 0.5)[1]])
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, Fu = Fu,
              E = E, S = S, cum_S = cum_S))
}

# Function for determining operating characteristics of a group sequential
# design with futility boundaries only
opchar_f <- function(delta, J, f, I, N, Sigma) {
  Fu      <- numeric(J)
  Fu[1]   <- mvtnorm::pmvnorm(-Inf, f[1], delta*sqrt(I[1]), sigma = 1)[1]
  for (j in 2:J) {
    Fu[j] <- mvtnorm::pmvnorm(c(f[1:(j - 1)], -Inf), c(rep(Inf, j - 1), f[j]),
                              delta*sqrt(I[1:j]), sigma = Sigma[1:j, 1:j])[1]
  }
  E       <- c(rep(0, J - 1), 1 - sum(Fu))
  cum_S   <- cumsum(S <- Fu + E)
  P       <- sum(E)
  ESS     <- sum(N*S)
  Med     <- ifelse(any(cum_S == 0.5),
                    0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1]),
                    N[which(cum_S > 0.5)[1]])
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, Fu = Fu,
              E = E, S = S, cum_S = cum_S))
}

# Function for determining operating characteristics of a group sequential
# design with efficacy boundaries
opchar_e <- function(delta, J, e, I, N, Sigma) {
  E      <- numeric(J)
  E[1]   <- mvtnorm::pmvnorm(lower = e[1], upper = Inf, mean = delta*sqrt(I[1]),
                             sigma = 1)[1]
  for (j in 2:J) {
    E[j] <- mvtnorm::pmvnorm(c(rep(-Inf, j - 1), e[j]), c(e[1:(j - 1)], Inf),
                             delta*sqrt(I[1:j]), sigma = Sigma[1:j, 1:j])[1]
  }
  Fu     <- c(rep(0, J - 1), 1 - sum(E))
  cum_S  <- cumsum(S <- Fu + E)
  P      <- sum(E)
  ESS    <- sum(N*S)
  Med    <- ifelse(any(cum_S == 0.5),
                   0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1]),
                   N[which(cum_S > 0.5)[1]])
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, Fu = Fu,
              E = E, S = S, cum_S = cum_S))
}

integrand <- function(y, z, k, delta, delta0, I) {
  return(dnorm(y - (delta - delta0)*sqrt(I[1]))*sqrt(I[2]/(I[2] - I[1]))*
           dnorm((z*sqrt(I[2]) - y*sqrt(I[1]) - (I[2] - I[1])*(delta - delta0))/
                   sqrt(I[2] - I[1])))
}

pdf_continuous <- function(z, k, delta, delta0, fb, e, I) {
  pdf <- numeric(length(z))
  for (i in 1:length(z)) {
    if (k == 1) {
      if (all(z[i] >= fb[1], z[i] <= e[1])) {
        pdf[i] <- 0
      } else {
        pdf[i] <- dnorm(z[i] - (delta - delta0)*sqrt(I[1]))
      }
    } else {
      pdf[i] <- integrate(f = integrand, lower = fb[1], upper = e[1], z = z[i],
                          k = k, delta = delta, delta0 = delta0, I = I)$value
    }
  }
  return(pdf)
}

bias_naive_continuous <- function(delta, delta0, f, e, I) {
  return(((I[2] - I[1])/(I[2]*sqrt(I[1]))) *
           (dnorm(e[1] - (delta - delta0)*sqrt(I[1])) -
              dnorm(f[1] - (delta - delta0)*sqrt(I[1]))))
}

est_naive_continuous <- function(z, k, delta0) {
  return(z/I[k]^0.5 + delta0)
}

est_bias_sub_continuous <- function(z, k, delta0, f, e, I) {
  return(z/I[k]^0.5 + delta0 -
           bias_naive_continuous(z/I[k]^0.5 + delta0, delta0, f, e, I))
}

est_bias_adj_continuous <- function(z, k, delta0, f, e, I) {

  int_est_bias_adj_continuous <- function(delta, delta0, fb, e, I, delta_naive) {
    return(delta_naive - bias_naive_continuous(delta, delta0, fb, e, I) -
             delta)
  }

  if (any(z <= -1e307, z >= 1e307)) {
    return(z/I[k]^0.5 + delta0)
  } else {
    return(stats::uniroot(f = int_est_bias_adj_continuous,
                          interval = c(-1e307, 1e307), delta0 = delta0, fb = f,
                          e = e, I = I,
                          delta_naive = z/I[k]^0.5 + delta0)$root)
  }
}

pval_continuous_so <- function(z, k, delta, delta0, f, e, I, Lambda) {
  if (k == 1) {
    return(pnorm(z, mean = (delta - delta0)*sqrt(I[1]), lower.tail = F))
  } else {
    return(pnorm(e[1], mean = (delta - delta0)*sqrt(I[1]), lower.tail = F) +
             mvtnorm::pmvnorm(lower = c(f[1], z), upper = c(e[1], Inf),
                              mean = (delta - delta0)*sqrt(I),
                              sigma = Lambda)[1])
  }
}

est_mue_continuous <- function(z, k, delta0, f, e, I, Lambda) {
  return(stats::uniroot(f = function(delta, z, k, delta0, fb, e, I, Lambda)
                          pval_continuous_so(z, k, delta, delta0, fb, e, I,
                                             Lambda) - 0.5,
                        interval = c(-1e307, 1e307),
                        z = z, k = k, delta0 = delta0, fb = f, e = e, I = I,
                        Lambda = Lambda)$root)
}

#naive <- sub <- adj <- mue <- cond <- numeric(100)
#k <- 2
#delta0 <- 0.1
#z <- seq(from = -10, to = 10, length.out = 100)
#f <- c(0.5, 2)
#e <- c(1, 2)
#I <- c(10, 20)
#Lambda <- matrix(c(1, sqrt(I[1]/I[2]), sqrt(I[1]/I[2]), 1), 2, 2)
#sigma <- 1
#n <- c(10, 20)
#for (i in 1:100) {
#  naive[i] <- est_naive_continuous(z[i], k, delta0)
#  sub[i]   <- est_bias_sub_continuous(z[i], k, delta0, f, e, I)
#  adj[i]   <- est_bias_adj_continuous(z[i], k, delta0, f, e, I)
#  mue[i]   <- est_mue_continuous(z[i], k, delta0, f, e, I, Lambda)
#  cond[i]  <- est_cond_continuous(z[i], k, delta0, f, e, n, sigma, I)
#}
#plot(z, naive, type="l")
#lines(z, sub, col = "green")
#lines(z, adj, col = "red")
#lines(z, mue, col = "blue")
#lines(z, cond, col = "yellow")


est_cond_continuous <- function(z, k, delta0, f, e, n, sigma, I) {

  int_est_cond_continuous <- function(delta, z, delta0, fb, e, n, sigma, I) {
    d_log_Ckmin1 <- numDeriv::grad(func =
                                     function(delta, delta0, fb, e, I)
                                       pnorm(e, mean = (delta - delta0)*sqrt(I)) -
                                         pnorm(fb, mean = (delta -
                                                            delta0)*sqrt(I)),
                                   x = delta, delta0 = delta0, fb = fb, e = e,
                                   I = I[1])
    return(n[2]*(z/I[2]^0.5 + delta0 - delta)/sigma^2 - d_log_Ckmin1)
  }

  if (k == 1) {
    return(z/I[k]^0.5 + delta0)
  } else {
    if (any(z <= -1e307, z >= 1e307)) {
      return(z/I[k]^0.5 + delta0)
    } else {
      return(stats::uniroot(f = int_est_cond_continuous,
                            interval = c(-1e307, 1e307), z = z, delta0 = delta0,
                            fb = f[1], e = e[1], n = n, sigma = sigma,
                            I = I)$root)
    }
  }
}

#s1 <- seq(-10, 10, length.out = 100)
#dens <- numeric(100)
#for (i in 1:100) {
#  dens[i] <- int_est_umvue_continuous(s1[i], z, delta, delta0, fb, e, n, I)
#}

#est_umvue_continuous(z,k,delta0,f,e,n,sigma,I)

est_umvue_continuous <- function(z, k, delta0, f, e, n, sigma, I) {
  # there's a problem here with the definition of Z. Need to go back through these
  # then we need a function which computes Bias of all of these estimators.
  # one which also computes RMSE. Then we just need UMVCUE

  int_est_umvue_continuous <- function(s1, z, delta, delta0, fb, e, n, I) {
    s      <- n[2]*(z/sqrt(I[2]) + delta0)
    z1     <- (s1/n[1] - delta0)*sqrt(I[1])
    Ckstar <- c(s - n[1]*(e[1]/sqrt(I[1]) + delta0),
                s - n[1]*(fb[1]/sqrt(I[1]) + delta0))
    return((s1/n[1])*pdf_continuous(((s - s1)/n[1] - delta0)*sqrt(I[1]),
                                    1, delta, delta0, Ckstar[1], Ckstar[2], I)*
             pdf_continuous(z1, 1, delta, delta0, fb, e, I))
  }

  if (k == 1) {
    return(z/sqrt(I[k]) + delta0)
  } else {
    return(integrate(f = int_est_umvue_continuous,
                     lower = n[1]*(f[1]/sqrt(I[1]) + delta0),
                     upper = n[1]*(e[1]/sqrt(I[1]) + delta0),
                     z = z, delta = 0, delta0 = delta0, fb = f, e = e, n = n[1],
                     I = I[1])/pdf_continuous(z, 2, 0, delta0, f, e, I))
  }
}
