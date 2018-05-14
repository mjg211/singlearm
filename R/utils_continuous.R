# Function to return covariance matrix in a group sequential design with a
# continuous outcome
cov_gs_continuous <- function(J, I){
  Sigma           <- diag(1, J, J)
  for (i in 2:J){
    for (j in 1:(i - 1)){
      Sigma[j, i] <- Sigma[i, j] <- sqrt(I[j]/I[i])
    }
  }
  return(Sigma)
}

# Function to return objective function in a group sequential design with
# efficacy and futility stopping boundaries
obj_fn_gs_continuous_ef <- function(parameters, J, delta, alpha, beta, sigma,
                                    optimality, equal_n, penalty) {
  if (equal_n) {
    n     <- rep(parameters[1], J)
    a     <- parameters[2:(J + 1)]
    r     <- c(parameters[2:J] + parameters[(J + 2):(2*J)], parameters[J + 1])
  } else {
    n     <- parameters[1:J]
    a     <- parameters[(J + 1):(2*J)]
    r     <- c(parameters[(J + 1):(2*J - 1)] + parameters[(2*J + 1):(3*J - 1)],
               parameters[2*J])
  }
  N       <- cumsum(n)
  I       <- N/(sigma^2)
  Sigma   <- cov_gs_continuous(J, I)
  perf_H0 <- opchar_ef(0, J, r, a, I, N, Sigma)
  perf_H1 <- opchar_ef(delta, J, r, a, I, N, Sigma)
  O       <- penalty*((alpha < perf_H0$P)*(perf_H0$P - alpha)/alpha +
                        (beta < 1 - perf_H1$P)*
                        ((1 - perf_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + perf_H0$ESS)
  } else if (optimality == "alt-ess") {
    return(O + perf_H1$ESS)
  } else {
    max_ess <- -suppressWarnings(optim(par = 0.5*delta, fn = max_ess_ef,
                                       J = J, r = r, a = a, I = I, N = N,
                                       Sigma = Sigma)$value)
    return(O + max_ess)
  }
}

# Function to return objective function in a group sequential design with
# efficacy stopping boundaries only
obj_fn_gs_continuous_e <- function(parameters, J, delta, alpha, beta, sigma,
                                   optimality, equal_n, penalty) {
  if (equal_n) {
    n     <- rep(parameters[1], J)
    r     <- parameters[2:(J + 1)]
  } else {
    n     <- parameters[1:J]
    r     <- parameters[(J + 1):(2*J)]
  }
  N       <- cumsum(n)
  I       <- N/(sigma^2)
  Sigma   <- cov_gs_continuous(J, I)
  perf_H0 <- opchar_e(0, J, r, I, N, Sigma)
  perf_H1 <- opchar_e(delta, J, r, I, N, Sigma)
  O       <- penalty*((alpha < perf_H0$P)*(perf_H0$P - alpha)/alpha +
                        (beta < 1 - perf_H1$P)*
                        ((1 - perf_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + perf_H0$ESS)
  } else if (optimality == "alt-ess") {
    return(O + perf_H1$ESS)
  } else {
    return(O + N[J])
  }
}

# Function to return objective function in a group sequential design with
# futility stopping boundaries only
obj_fn_gs_continuous_f <- function(parameters, J, delta, alpha, beta, sigma,
                                   optimality, equal_n, penalty) {
  if (equal_n) {
    n     <- rep(parameters[1], J)
    a     <- parameters[2:(J + 1)]
  } else {
    n     <- parameters[1:J]
    a     <- parameters[(J + 1):(2*J)]
  }
  N       <- cumsum(n)
  I       <- N/(sigma^2)
  Sigma   <- cov_gs_continuous(J, I)
  perf_H0 <- opchar_f(0, J, a, I, N, Sigma)
  perf_H1 <- opchar_f(delta, J, a, I, N, Sigma)
  O       <- penalty*((alpha < perf_H0$P)*(perf_H0$P - alpha)/alpha +
                        (beta < 1 - perf_H1$P)*
                        ((1 - perf_H1$P) - beta)/beta)
  if (optimality == "null-ess") {
    return(O + perf_H0$ESS)
  } else if (optimality == "alt-ess") {
    return(O + perf_H1$ESS)
  } else {
    return(O + N[J])
  }
}

# Function used to find max_delta{ESS(delta)} for a group sequential design with
# efficacy and futility boundaries
max_ess_ef <- function(delta, J, r, a, I, N, Sigma){
  return(-opchar_ef(delta, J, r, a, I, N, Sigma)[2])
}

# Function for determining operating characteristics of a group sequential
# design with efficacy and futility boundaries
opchar_ef <- function(delta, J, r, a, I, N, Sigma){
  R      <- numeric(J)
  A      <- numeric(J)
  R[1]   <- mvtnorm::pmvnorm(lower = r[1], upper = Inf, mean = delta*sqrt(I[1]),
                             sigma = 1)[1]
  A[1]   <- mvtnorm::pmvnorm(lower = -Inf, upper = a[1],
                             mean = delta*sqrt(I[1]), sigma = 1)[1]
  for (j in 2:J){
    R[j] <- mvtnorm::pmvnorm(lower = c(a[1:(j - 1)], r[j]),
                             upper = c(r[1:(j - 1)], Inf),
                             mean = delta*sqrt(I[1:j]),
                             sigma = Sigma[1:j, 1:j])[1]
    A[j] <- mvtnorm::pmvnorm(lower = c(a[1:(j - 1)], -Inf),
                             upper = c(r[1:(j - 1)], a[j]),
                             mean = delta*sqrt(I[1:j]),
                             sigma = Sigma[1:j, 1:j])[1]
  }
  S      <- A + R
  P      <- sum(R)
  ESS    <- sum(N*S)
  cum_S  <- cumsum(S)
  if (any(cum_S == 0.5)){
    Med  <- 0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1])
  } else {
    Med  <- N[which(cum_S > 0.5)[1]]
  }
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, A = A,
              R = R, S = S, cum_S = cum_S))
}

# Function for determining operating characteristics of a group sequential
# design with futility boundaries only
opchar_f <- function(delta, J, a, I, N, Sigma){
  A      <- numeric(J)
  A[1]   <- mvtnorm::pmvnorm(lower = -Inf, upper = a[1],
                             mean = delta*sqrt(I[1]), sigma = 1)[1]
  for (j in 2:J){
    A[j] <- mvtnorm::pmvnorm(lower = c(a[1:(j - 1)], -Inf),
                             upper = c(rep(Inf, j - 1), a[j]),
                             mean = delta*sqrt(I[1:j]),
                             sigma = Sigma[1:j, 1:j])[1]
  }
  R      <- c(rep(0, J - 1), 1 - sum(A))
  S      <- A + R
  P      <- sum(R)
  ESS    <- sum(N*S)
  cum_S  <- cumsum(S)
  if (any(cum_S == 0.5)){
    Med  <- 0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1])
  } else {
    Med  <- N[which(cum_S > 0.5)[1]]
  }
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, A = A,
              R = R, S = S, cum_S = cum_S))
}

# Function for determining operating characteristics of a group sequential
# design with efficacy boundaries
opchar_e <- function(delta, J, r, I, N, Sigma){
  R      <- numeric(J)
  R[1]   <- mvtnorm::pmvnorm(lower = r[1], upper = Inf, mean = delta*sqrt(I[1]),
                             sigma = 1)[1]
  for (j in 2:J){
    R[j] <- mvtnorm::pmvnorm(lower = c(rep(-Inf, j - 1), r[j]),
                             upper = c(r[1:(j - 1)], Inf),
                             mean = delta*sqrt(I[1:j]),
                             sigma = Sigma[1:j, 1:j])[1]
  }
  A      <- c(rep(0, J - 1), 1 - sum(R))
  S      <- A + R
  P      <- sum(R)
  ESS    <- sum(N*S)
  cum_S  <- cumsum(S)
  if (any(cum_S == 0.5)){
    Med  <- 0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1])
  } else {
    Med  <- N[which(cum_S > 0.5)[1]]
  }
  return(list(P = P, ESS = ESS, VSS = sum(N^2*S) - ESS^2, Med = Med, A = A,
              R = R, S = S, cum_S = cum_S))
}
