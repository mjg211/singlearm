# Function for determining operating characteristics of an adaptive design
int_opchar_gehan <- function(pi, a1, r1, n1, n2, k, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi      <- pmf_adaptive(pi, a1, r1, n1, n2, k)
  }
  int_pi        <- pi
  len_pi        <- length(int_pi)
  opchar        <- matrix(0, nrow = len_pi, ncol = 8)
  poss_n        <- sort(unique(c(n1, n1 + n2)))
  prob_n        <- numeric(len_poss_n <- length(poss_n))
  for (i in 1:len_pi) {
    S1          <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] &
                                       k == 1)$`f(s1,s,m|pi)`)
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
    opchar[i, ] <- c(int_pi[i], ESS, sum(poss_n^2*prob_n) - ESS^2, Med, S1,
                     1 - S1, S1, 1 - S1)
  }
  colnames(opchar) <- c("pi", "ESS", "VSS", "Med", "S1", "S2", "cum(S1)",
                        "cum(S2)")
  return(tibble::as_tibble(opchar))
}

gehan_dc_ef <- function(pi0, pi1, alpha, n1, n2, dc_ef, dc_pf,
                        length_dc_ef) {
  nomalpha    <<- alpha
  n1          <<- n1
  dbinom_pi0  <<- stats::dbinom(0:n1, n1, pi0)
  dbinom_pi1  <<- stats::dbinom(0:n1, n1, pi1)
  elem        <<- cbind(1:(n1 + 1), 1)
  store       <<- NULL
  maxpower    <<- 0
  combination <<- NULL
  dc_ef        <<- dc_ef
  dc_pf        <<- dc_pf
  length_dc_ef <<- length_dc_ef
  branch_gehan(1, 1)
  if (is.null(combination)) {
    return(list())
  } else {
    D0 <- D1  <- numeric(n1 + 1)
    for (i in 1:(n1 + 1)) {
      D0[i]   <- dc_ef[[i]][combination[i, 2]]
      D1[i]   <- dc_pf[[i]][combination[i, 2]]
    }
    a1        <- suppressWarnings(max(which(D0 == 0)) - 1)
    r1        <- sign(suppressWarnings(min(which(D0 == 1)) - 1))*
                        suppressWarnings(min(which(D0 == 1)) - 1)
    c2        <- rep(NA, n1 + 1)
    c2[which(!(D0 %in% c(0, 1)))] <-
      stats::qbinom(D0[which(!(D0 %in% c(0, 1)))],
                    n2[which(!(D0 %in% c(0, 1)))], pi0, lower.tail = F) + 1
    r2        <- 0:n1 + c2
    a2        <- r2 - 1
    return(list(D = D0, a1 = a1, r1 = r1, a2 = a2, r2 = r2))
  }
}

# Branching function for Shan et al. adaptive design
branch_gehan <- function(k, j) {
  if (k <= n1) {
    for (jiter in 1:length_dc_ef[k + 1]) {
      elem[k + 1, 2] <<- jiter
      if (bound_gehan(k + 1, jiter)){
        branch_gehan(k + 1, jiter)
      }
    }
  } else {
    final_dc_pf      <- numeric(n1 + 1)
    for (i in 1:(n1 + 1)) {
      final_dc_pf[i] <- dc_pf[[i]][elem[i, 2]]
    }
    maxpower         <<- sum(final_dc_pf*dbinom_pi1)
    combination      <<- elem
  }
}

# Bounding function for Shan et al. adaptive design
bound_gehan <- function(k, j) {
  local_dc_ef      <- local_dc_pf <- numeric(n1 + 1)
  for (i in 1:(n1 + 1)) {
    local_dc_ef[i] <- dc_ef[[i]][elem[i, 2]]
    local_dc_pf[i] <- dc_pf[[i]][elem[i, 2]]
  }
  if (k <= n1) {
    if (any(local_dc_ef[2:k] < local_dc_ef[1:(k - 1)])) {
      return(F)
    } else if (sum(local_dc_ef[1:k]*dbinom_pi0[1:k]) > nomalpha) {
      return(F)
    } else if (sum(local_dc_pf[1:k]*dbinom_pi1[1:k]) +
                 sum(dbinom_pi1[(k + 1):(n1 + 1)]) < maxpower) {
      return(F)
    }
  } else {
    if (any(local_dc_ef[2:(n1 + 1)] < local_dc_ef[1:n1])) {
      return(F)
    } else if (sum(local_dc_ef*dbinom_pi0) > nomalpha) {
      return(F)
    } else if (sum(local_dc_pf*dbinom_pi1) < maxpower) {
      return(F)
    }
  }
  return(T)
}
