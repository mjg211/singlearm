# Function for determining pmf of an adaptive design
pmf_adaptive <- function(pi, a1, r1, n1, n2, k) {
  int_pi        <- pi
  len_pi        <- length(int_pi)
  if (2 %in% k) {
    poss_n      <- c(n1, n2)
  } else {
    poss_n      <- n1
  }
  unique_n      <- unique(poss_n)
  unique_n      <- unique_n[which(unique_n > 0)]
  pmf_fixed_pi  <- tibble::tibble(pi = rep(int_pi,
                                           sum(unique_n) + length(unique_n)),
                                  s = rep(unlist(sapply(unique_n,
                                                        function(n) 0:n)),
                                          each = len_pi),
                                  m = rep(unlist(sapply(unique_n,
                                                        function(n)
                                                          rep(n, n + 1))),
                                          each = len_pi),
                                  f = stats::dbinom(s, m, pi))
  terminal      <- terminal_states_adaptive(a1, r1, n1, n2, k)
  rows_terminal <- nrow(terminal)
  pmf           <- tibble::tibble(pi = rep(int_pi, each = rows_terminal),
                                  s1 = rep(terminal$s1, len_pi),
                                  s = rep(terminal$s, len_pi),
                                  m = rep(terminal$m, len_pi),
                                  k = factor(rep(terminal$k, len_pi), k),
                                  `f(s1,s,m|pi)` = NA)
  for (i in 1:len_pi) {
    if (int_pi[i] == 0) {
      f_i <- numeric(rows_terminal)
      if (any(terminal$s == 0)) {
        f_i[which(terminal$s == 0)] <- 1
      }
    } else if (int_pi[i] == 1) {
      f_i <- numeric(rows_terminal)
      if (any(terminal$s == terminal$m)) {
        f_i[which(terminal$s == terminal$m)] <- 1
      }
    } else {
      f_i       <- NULL
      if (1 %in% k) {
        f_i     <- dplyr::filter(pmf_fixed_pi, pi == int_pi[i] & m == n1 &
                                   (s <= a1 | s >= r1))$f
      }
      if (2 %in% k) {
        if (a1 < r1 & a1 < n1) {
          for (s1 in (a1 + 1):min(n1, r1 - 1)) {
            f_i <- c(f_i, dplyr::filter(pmf_fixed_pi, pi == int_pi[i] &
                                          m == n1 & s == s1)$f*
                       dplyr::filter(pmf_fixed_pi, pi == int_pi[i] &
                                       m == n2[s1 + 1])$f)
          }
        }
      }
      f_i       <- f_i/sum(f_i)
    }
    pmf$`f(s1,s,m|pi)`[which(pmf$pi == int_pi[i])] <- f_i
  }
  return(pmf)
}

# Function for determining operating characteristics of an adaptive design
int_opchar_adaptive <- function(pi, a1, r1, a2, r2, n1, n2, k, pmf_pi){
  if (missing(pmf_pi)) {
    pmf_pi <- pmf_adaptive(pi, a1, r1, n1, n2, k)
  }
  int_pi   <- pi
  len_pi   <- length(int_pi)
  opchar   <- matrix(0, nrow = len_pi, ncol = 13)
  poss_n   <- sort(unique(c(n1, n1 + n2)))
  prob_n   <- numeric(len_poss_n <- length(poss_n))
  for (i in 1:len_pi) {
    A <- R <- numeric(2)
    A[1]   <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & s <= a1 &
                                  k == 1)$`f(s1,s,m|pi)`)
    R[1]   <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & s >= r1 &
                                  k == 1)$`f(s1,s,m|pi)`)
    for (int_s1 in which(n2 != 0)) {
      A[2] <- A[2] + sum(dplyr::filter(pmf_pi, (pi == int_pi[i]) &
                                         (s1 == int_s1 - 1) &
                                         (s <= a2[int_s1]) &
                                         (m == n1 + n2[int_s1]) &
                                         (k == 2))$`f(s1,s,m|pi)`)
      R[2] <- R[2] + sum(dplyr::filter(pmf_pi, (pi == int_pi[i]) &
                                         (s1 == int_s1 - 1) &
                                         (s >= r2[int_s1]) &
                                         (m == n1 + n2[int_s1]) &
                                         (k == 2))$`f(s1,s,m|pi)`)
    }
    cum_S       <- cumsum(S <- A + R)
    for (j in 1:len_poss_n) {
      prob_n[j] <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] &
                                       m == poss_n[j])$`f(s1,s,m|pi)`)
    }
    cum_prob_n  <- cumsum(prob_n)
    if (any(cum_prob_n == 0.5)){
      Med       <- 0.5*(poss_n[which(cum_prob_n == 0.5)] +
                          poss_n[which(cum_prob_n == 0.5) + 1])
    } else {
      Med       <- poss_n[which(cum_prob_n > 0.5)[1]]
    }
    ESS         <- sum(poss_n*prob_n)
    opchar[i, ] <- c(int_pi[i], sum(R), ESS, sum(poss_n^2*prob_n) - ESS^2, Med,
                     A, R, S, cum_S)
  }
  colnames(opchar) <- c("pi", "P(pi)", "ESS(pi)", "VSS(pi)", "Med(pi)",
                        paste(rep(c("A", "R", "S"), each = 2), rep(1:2, 3),
                              "(pi)", sep = ""),
                        paste("cum(S", 1:2, "(pi))", sep = ""))
  return(tibble::as_tibble(opchar))
}

# Function to determine terminal states in an adaptive design
terminal_states_adaptive <- function(a1, r1, n1, n2, k) {
  terminal     <- NULL
  if (1 %in% k) {
    if (a1 >= 0) {
      terminal   <- rbind(terminal, cbind(0:a1, 0:a1, n1, 1))
    }
    if (r1 < Inf) {
      terminal   <- rbind(terminal, cbind(r1:n1, r1:n1, n1, 1))
    }
  }
  if (2 %in% k) {
    if (a1 < r1 & a1 < n1) {
      for (s1 in max(0, a1 + 1):min(r1 - 1, n1)) {
        terminal <- rbind(terminal, cbind(s1, s1:(s1 + n2[s1 + 1]),
                                          n1 + n2[s1 + 1], 2))
      }
    }
  }
  return(tibble::tibble(s1 = as.integer(terminal[, 1]),
                        s = as.integer(terminal[, 2]),
                        m = as.integer(terminal[, 3]),
                        k = factor(terminal[, 4], k)))
}

# Function to compute stage 2 sample sizes for a particular first stage sample
# size in an adaptive design
des_adaptive_n2forn1 <- function(pi0, pi1, alpha, beta, n1, n2min, n2max,
                                 monotonic, w) {
  nomalpha      <<- alpha
  nombeta       <<- beta
  n1            <<- n1
  w             <<- w
  dbinomial_pi0 <<- stats::dbinom(0:n1, n1, pi0)
  dbinomial_pi1 <<- stats::dbinom(0:n1, n1, pi1)
  if (monotonic) {
    dc_ef                <<- NULL
    dc_pf                <<- NULL
    dc_ss                <<- NULL
    for (n in n2max:n2min) {
      dc_ef_n            <-  matrix(0, nrow = n1 + 1, ncol = n + 1)
      dc_pf_n            <-  matrix(0, nrow = n1 + 1, ncol = n + 1)
      dc_ss_n            <-  matrix(0, nrow = n1 + 1, ncol = n + 1)
      for (s in 0:n1) {
        dc_ef_n[s + 1, ] <-  1 - stats::pbinom(n:0 - s, n, pi0)
        dc_pf_n[s + 1, ] <-  1 - stats::pbinom(n:0 - s, n, pi1)
        dc_ss_n[s + 1, ] <-  n1 + n
      }
      dc_ef              <<- cbind(dc_ef, dc_ef_n)
      dc_pf              <<- cbind(dc_pf, dc_pf_n)
      dc_ss              <<- cbind(dc_ss, dc_ss_n)
    }
    dc_ef                <<- cbind(0, dc_ef, 1)
    dc_pf                <<- cbind(0, dc_pf, 1)
    dc_ss                <<- cbind(n1, dc_ss, n1)
  } else {
    quantiles <-  unlist(sapply(n2min:n2max, function(x){-1:(x - 1)}))
    sam_sizes <-  rep(n2min:n2max, n2min:n2max + 1)
    dc_ef     <<- stats::pbinom(quantiles, sam_sizes, pi0, lower.tail = F)
    dc_pf     <<- stats::pbinom(quantiles, sam_sizes, pi1, lower.tail = F)
    dc_ss     <<- rep(n1 + n2min:n2max, n2min:n2max + 1)
    dc        <<- cbind(dc_ef, dc_pf, dc_ss)
    dc        <<- dc[order(dc[, 1]), ]
    dc_ef     <<- matrix(c(0, dc[, 1], 1), nrow = n1 + 1, ncol = nrow(dc) + 2,
                         byrow = T)
    dc_pf     <<- matrix(c(0, dc[, 2], 1), nrow = n1 + 1, ncol = nrow(dc) + 2,
                         byrow = T)
    dc_ss     <<- matrix(c(n1, dc[, 3], n1), nrow = n1 + 1, ncol = nrow(dc) + 2,
                         byrow = T)
  }
  card_P2     <<- 0.5*(n2max^2 - n2min^2) + 1.5*n2max - 0.5*n2min + 3
  elem        <<- cbind(1:(n1 + 1), 1)
  store       <<- NULL
  en          <<- n1 + n2max
  combination <<- NULL
  if (monotonic) {
    branch_shan(1, 1)
  } else {
    branch_ek(1, 1)
  }
  if (is.null(combination)){
    return(numeric(14 + 3*(n1 + 1)))
  } else {
    D0      <- dc_ef[combination]
    D1      <- dc_pf[combination]
    n2      <- dc_ss[combination] - n1
    a1      <- suppressWarnings(max(which(D0 == 0)) - 1)
    r1      <- sign(suppressWarnings(min(which(D0 == 1)) - 1))*
                 suppressWarnings(min(which(D0 == 1)) - 1)
    c2      <- rep(NA, n1 + 1)
    c2[which(!(D0 %in% c(0, 1)))] <- stats::qbinom(D0[which(!(D0 %in% c(0, 1)))],
                                            n2[which(!(D0 %in% c(0, 1)))], pi0,
                                            lower.tail = F) + 1
    r2      <- 0:n1 + c2
    a2      <- r2 - 1
    P_pi0   <- sum(D0*dbinomial_pi0)
    P_pi1   <- sum(D1*dbinomial_pi1)
    ESS_pi0 <- sum((n1 + n2)*dbinomial_pi0)
    ESS_pi1 <- sum((n1 + n2)*dbinomial_pi1)
    PET_pi0 <- sum(dbinomial_pi0[which(D0 %in% c(0, 1))])
    PET_pi1 <- sum(dbinomial_pi1[which(D1 %in% c(0, 1))])
    order_N      <- sort(n1 + n2)
    cum_prob_pi0 <- cumsum(dbinomial_pi0[order(n1 + n2)])
    cum_prob_pi1 <- cumsum(dbinomial_pi1[order(n1 + n2)])
    if (any(cum_prob_pi0 == 0.5)){
      Med_pi0    <- 0.5*(order_N[which(cum_prob_pi0 == 0.5)] +
                           order_N[which(cum_prob_pi0 == 0.5) + 1])
    } else {
      Med_pi0    <- order_N[which(cum_prob_pi0 > 0.5)[1]]
    }
    if (any(cum_prob_pi1 == 0.5)){
      Med_pi1    <- 0.5*(order_N[which(cum_prob_pi1 == 0.5)] +
                           order_N[which(cum_prob_pi1 == 0.5) + 1])
    } else {
      Med_pi1    <- order_N[which(cum_prob_pi1 > 0.5)[1]]
    }
    VSS_pi0      <- sum((n1 + n2)^2*dbinomial_pi0) - ESS_pi0^2
    VSS_pi1      <- sum((n1 + n2)^2*dbinomial_pi1) - ESS_pi1^2
    maxN         <- max(dc_ss[combination])
    return(c(n1, n2, a1, r1, a2, r2, P_pi0, P_pi1, ESS_pi0, ESS_pi1, PET_pi0,
             PET_pi1, Med_pi0, Med_pi1, VSS_pi0, VSS_pi1, maxN))
  }
}

# Branching function for Englert-Keiser adaptive design
branch_ek <- function(k, j) {
  if (k <= n1) {
    for (jiter in j:card_P2) {
      elem[k + 1, 2] <<- jiter
      if (bound_ek(k + 1, jiter)){
        branch_ek(k + 1, jiter)
      }
    }
  } else {
    en          <<- sum(dc_ss[elem]*dbinomial_pi0)
    combination <<- elem
  }
}

# Branching function for Shan et al. adaptive design
branch_shan <- function(k, j) {
  if (k <= n1) {
    for (jiter in j:card_P2) {
      elem[k + 1, 2] <<- jiter
      if (bound_shan(k + 1, jiter)){
        branch_shan(k + 1, jiter)
      }
    }
  } else {
    en          <<- sum(dc_ss[elem]*dbinomial_pi0)
    combination <<- elem
  }
}

# Bounding function for Englert-Keiser adaptive design
bound_ek <- function(k, j) {
  if (k <= n1) {
    if (sum(dc_ef[elem[1:k, ]]*dbinomial_pi0[1:k]) +
          dc_ef[j]*sum(dbinomial_pi0[(k + 1):(n1 + 1)]) > nomalpha) {
      return(F)
    } else if (sum(dc_pf[elem[1:k, ]]*dbinomial_pi1[1:k]) +
                 sum(dbinomial_pi1[(k + 1):(n1 + 1)]) < 1 - nombeta) {
      return(F)
    } else if (w[1]*(sum(dc_ss[elem[1:k, ]]*dbinomial_pi0[1:k]) +
                       n1*sum(dbinomial_pi0[(k + 1):(n1 + 1)])) +
                 w[2]*(sum(dc_ss[elem[1:k, ]]*dbinomial_pi1[1:k]) +
                         n1*sum(dbinomial_pi1[(k + 1):(n1 + 1)])) +
                 w[3]*max(dc_ss[elem[1:k, ]]) > en) {
      return(F)
    }
  } else {
    if (sum(dc_ef[elem]*dbinomial_pi0) > nomalpha) {
      return(F)
    } else if (sum(dc_pf[elem]*dbinomial_pi1) < 1 - nombeta) {
      return(F)
    } else if (w[1]*sum(dc_ss[elem]*dbinomial_pi0) +
                 w[2]*sum(dc_ss[elem]*dbinomial_pi1) +
                 w[3]*max(dc_ss[elem]) > en) {
      return(F)
    }
  }
  return(T)
}

# Bounding function for Shan et al. adaptive design
bound_shan <- function(k, j) {
  if (k <= n1) {
    if (sum(dc_ef[elem[1:k, ]]*dbinomial_pi0[1:k]) > nomalpha) {
      return(F)
    } else if (sum(dc_pf[elem[1:k, ]]*dbinomial_pi1[1:k]) +
                 sum(dbinomial_pi1[(k + 1):(n1 + 1)]) < 1 - nombeta) {
      return(F)
    } else if (w[1]*(sum(dc_ss[elem[1:k, ]]*dbinomial_pi0[1:k]) +
                       n1*sum(dbinomial_pi0[(k + 1):(n1 + 1)])) +
                 w[2]*(sum(dc_ss[elem[1:k, ]]*dbinomial_pi1[1:k]) +
                         n1*sum(dbinomial_pi1[(k + 1):(n1 + 1)])) +
                 w[3]*max(dc_ss[elem[1:k, ]]) > en) {
      return(F)
    }
  } else {
    if (sum(dc_ef[elem]*dbinomial_pi0) > nomalpha) {
      return(F)
    } else if (sum(dc_pf[elem]*dbinomial_pi1) < 1 - nombeta) {
      return(F)
    } else if (w[1]*sum(dc_ss[elem]*dbinomial_pi0) +
                 w[2]*sum(dc_ss[elem]*dbinomial_pi1) +
                 w[3]*max(dc_ss[elem])  > en) {
      return(F)
    }
  }
  return(T)
}
