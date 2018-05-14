# Function for determining pmf of a bivariate design
pmf_bivariate <- function(piR, piT, theta, J, aR, aT, n, dbivariate_piRT) {
  if (missing(dbivariate_piRT)) {
    dbivariate_piRT        <- list()
    for (j in 1:J) {
      dbivariate_piRT[[j]] <- dbivariate(n[j], piR, piT, theta)
    }
  }
  pdf_mat          <- list()
  pdf_mat[[1]]     <- dbivariate_piRT[[1]]
  if (J > 1) {
    for (j in 2:J) {
      pdf_mat[[j]] <- matrix(0, sum(n[1:j]) + 1, sum(n[1:j]) + 1)
      for (sR in (max(aR[1:(j - 1)]) + 1):sum(n[1:(j - 1)])) {
        for (sT in 0:(min(aT[1:2]) - 1)) {
          pdf_mat[[j]][(sR + 1):(sR + 1 + n[j]), (sT + 1):(sT + 1 + n[j])] <-
            pdf_mat[[j]][(sR + 1):(sR + 1 + n[j]), (sT + 1):(sT + 1 + n[j])] +
              dbivariate_piRT[[j]]*pdf_mat[[j - 1]][sR + 1, sT + 1]
        }
      }
    }
  }
  terminal <- terminal_states_bivariate(J, aR, aT, n)
  f        <- numeric(nrow(terminal))
  for (i in 1:nrow(terminal)) {
    f[i]   <- pdf_mat[[terminal$k[i]]][terminal$sR[i] + 1, terminal$sT[i] + 1]
  }
  return(tibble::tibble(piR = piR, piT = piT, theta = theta, sR = terminal$sR,
                        sT = terminal$sT, m = terminal$m, k = terminal$k,
                        f = f))
}

# Function to determine density function of a bivariate distribution
dbivariate <- function(n, piR, piT, theta) {
  if (theta == 0) {
    if (piR + piT < 1) {
      pi   <- c(0, piR, piT, 1 - piR - piT)
    } else if (piR + piT == 1) {
      pi   <- c(0, piR, piT, 0)
    } else {
      pi   <- c(piT - (1 - piR), 1 - piT, 1 - piR, 0)
    }
  } else if (theta == Inf) {
    if (piR == piT) {
      pi   <- c(piR, 0, 0, 1 - piR)
    } else if (piR > piT) {
      pi   <- c(piT, piR - piT, 0, 1- piR)
    } else {
      pi   <- c(piR, 0, piT - piR, 1 - piT)
    }
  } else {
    a      <- theta - 1
    b      <- (piR + piT)*(1 - theta) - 1
    c      <- theta*piR*piT
    if (a == 0) {
      piRT <- -c/b
      pi   <- c(piRT, piR - piRT, piT - piRT, 1 + piRT - piR - piT)
    } else {
      piRT <- (-b +c(-1, 1)*sqrt(b^2 - 4*a*c))/(2*a)
      pi   <- rbind(c(piRT[1], piR - piRT[1], piT - piRT[1],
                      1 + piRT[1] - piR - piT),
                    c(piRT[2], piR - piRT[2], piT - piRT[2],
                      1 + piRT[2] - piR - piT))
      if (all(pi[1, ] >= 0)) {
        pi <- pi[1, ]
      } else if (all(pi[2, ] >= 0)) {
        pi <- pi[2, ]
      } else {
        stop("error: pi")
      }
    }
  }
  dbivariate_piRT                             <- matrix(0, n + 1, n + 1)
  if (piR == 0 & piT == 0) {
    dbivariate_piRT[1, 1]                     <- 1
  } else if (all(pi > 0) | all(pi[1:3] > 0, pi[4] == 0)) {
    for (sR in 0:n) {
      for (sT in 0:n) {
        dbivariate_piRT[sR + 1, sT + 1]       <- stats::dbinom(sT, n, piT)*
                                                   sum(stats::dbinom(0:min(sR, sT), sT,
                                                              pi[1]/piT)*
                                                         stats::dbinom(sR - 0:min(sR, sT),
                                                          n - sT,
                                                          pi[2]/(1 - piT)))
      }
    }
  } else if (piR == 0) {
    dbivariate_piRT[1, 1:(n + 1)]             <- stats::dbinom(0:n, n, piT)
  } else if (piT == 0) {
    dbivariate_piRT[1:(n + 1), 1]             <- stats::dbinom(0:n, n, piT)
  } else if (pi[2] == pi[3] & pi[2] == 0) {
    for (s in 0:n) {
      dbivariate_piRT[s + 1, s + 1]           <- stats::dbinom(s, n, piR)
    }
  } else if (pi[1] == pi[4] & pi[1] == 0) {
    for (s in 0:n) {
      dbivariate_piRT[s + 1, n - s + 1]       <- stats::dbinom(s, n, piR)
    }
  } else if (pi[1] == 0) {
    for (sT in 0:n) {
      dbivariate_piRT[1:(n - sT + 1), sT + 1] <- stats::dbinom(sT, n, piT)*
        stats::dbinom(0:(n - sT), n - sT,
                                                          piR/(1 - piT))
    }
  } else if (pi[2] == 0) {
    for (sT in 0:n) {
      dbivariate_piRT[1:(sT + 1), sT + 1]     <- stats::dbinom(sT, n, piT)*
        stats::dbinom(0:sT, sT, piR/piT)
    }
  } else if (pi[3] == 0) {
    for (sR in 0:n) {
      dbivariate_piRT[sR + 1, 1:(sR + 1)]     <- stats::dbinom(sR, n, piR)*
        stats::dbinom(0:sR, sR, piT/piR)
    }
  }
  return(dbivariate_piRT)
}

# Function for determining terminal states in a bivariate design
terminal_states_bivariate <- function(J, aR, aT, n) {
  if (J == 1) {
    terminal <- expand.grid(0:n[1], 0:n[1], n[1], 1)
  } else if (J == 2) {
    stage_1  <- expand.grid(0:n[1], 0:n[1], n[1], 1)
    stage_2  <- expand.grid((aR[1] + 1):sum(n), 0:(aT[1] - 1 + n[2]),
                            sum(n[1:2]), 2)
    terminal <- rbind(stage_1[which(stage_1[, 1] <= aR[1] |
                                      stage_1[, 2] >= aT[1]), ], stage_2)
  } else {
    stage_1  <- expand.grid(0:n[1], 0:n[1], n[1], 1)
    stage_2  <- expand.grid((aR[1] + 1):(n[1] + n[2]), 0:(aT[1] - 1 + n[2]),
                            sum(n[1:2]), 2)
    stage_3  <- expand.grid((max(aR[1:2]) + 1):sum(n[1:2]),
                            0:(min(aT[1:2]) - 1 + n[3]), sum(n), 3)
    terminal <- rbind(stage_1[which(stage_1[, 1] <= aR[1] |
                                      stage_1[, 2] >= aT[1]), ],
                      stage_2[which(stage_2[, 1] <= aR[2] |
                                      stage_2[, 2] >= aT[2]), ], stage_3)
  }
  return(tibble::tibble(sR = as.integer(terminal[, 1]),
                        sT = as.integer(terminal[, 2]),
                        m = as.integer(terminal[, 3]),
                        k = factor(terminal[, 4], 1:J)))
}
