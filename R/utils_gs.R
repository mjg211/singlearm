# Function for determining pmf of a group sequential design
pmf_gs <- function(pi, J, a, r, n, k) {
  max_k        <- max(k)
  int_pi       <- pi
  len_pi       <- length(int_pi)
  unique_n     <- unique(n[1:max_k])
  pmf_fixed_pi <- tibble::tibble(pi = rep(int_pi,
                                          sum(unique_n) + length(unique_n)),
                                 s = rep(unlist(sapply(unique_n,
                                                       function(n) 0:n)),
                                         each = len_pi),
                                 m = rep(unlist(sapply(unique_n,
                                                       function(n) rep(n,
                                                                       n + 1))),
                                         each = len_pi),
                                 `f(s,m|pi)` = dbinom(s, m, pi))
  terminal      <- terminal_states_gs(J, a, r, n, k)
  rows_terminal <- nrow(terminal)
  pmf           <- tibble::tibble(pi = rep(int_pi, each = rows_terminal),
                                  s = rep(terminal$s, len_pi),
                                  m = rep(terminal$m, len_pi),
                                  k = factor(rep(terminal$k, len_pi), k),
                                  `f(s,m|pi)` = NA)
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
      pmf_mat                  <- matrix(0, nrow = sum(n[1:max_k]) + 1,
                                         ncol = max_k)
      pmf_mat[1:(n[1] + 1), 1] <- dplyr::filter(pmf_fixed_pi, pi == int_pi[i] &
                                                  m == n[1])$`f(s,m|pi)`
      cont                     <- c(max(0, a[1] + 1), min(r[1] - 1, n[1]))
      if (max_k > 1) {
        for (j in 1:(max_k - 1)) {
          for (l in cont[1]:cont[2]) {
            pmf_mat[(l + 1):(l + 1 + n[j + 1]), j + 1] <-
              pmf_mat[(l + 1):(l + 1 + n[j + 1]), j + 1] +
                pmf_mat[l + 1, j]*dplyr::filter(pmf_fixed_pi, pi == int_pi[i] &
                                                  m == n[j + 1])$`f(s,m|pi)`
          }
          pmf_mat[(cont[1] + 1):(cont[2] + 1), j] <- 0
          new                                     <- cont[1]:(cont[2] +
                                                                n[j + 1])
          cont                                    <- c(min(new[new > a[j + 1]]),
                                                       max(new[new < r[j + 1]]))
        }
      }
      if (max_k < J) {
        pmf_mat[(cont[1] + 1):(cont[2] + 1), max_k] <- 0
      }
      if (length(k) < J) {
        pmf_mat[, -k] <- 0
        pmf_mat[, k]  <- pmf_mat[, k]/sum(pmf_mat[, k])
      }
      f_i <- pmf_mat[which(pmf_mat > 0)]
    }
    pmf$`f(s,m|pi)`[which(pmf$pi == int_pi[i])] <- f_i
  }
  return(pmf)
}

# Function for determining operating characteristics of group sequential design
int_opchar_gs <- function(pi, J, a, r, n, N, k, summary, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi         <- pmf_gs(pi, J, a, r, n, k)
  }
  int_pi           <- pi
  len_pi           <- length(int_pi)
  opchar           <- matrix(0, nrow = len_pi, ncol = 6 + 4*J)
  for (i in 1:len_pi) {
    A <- R         <- numeric(J)
    for (j in k) {
      A[j]         <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & s <= a[j] &
                                          k == j)$`f(s,m|pi)`)
      R[j]         <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & s >= r[j] &
                                          k == j)$`f(s,m|pi)`)
    }
    cum_S          <- cumsum(S <- A + R)
    if (any(cum_S == 0.5)) {
      Med          <- 0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1])
    } else {
      Med          <- N[which(cum_S > 0.5)[1]]
    }
    opchar[i, ]    <- c(int_pi[i], sum(R), sum(N*S), sum(N^2*S) - sum(N*S)^2,
                        Med, A, R, S, cum_S, N[J])
  }
  if (all(summary, i%%1000 == 0)) {
    message("...performance for ", i, " elements of pi evaluated...")
  }
  colnames(opchar)  <- c("pi", "P(pi)", "ESS(pi)", "VSS(pi)", "Med(pi)",
                         paste(rep(c("A", "R", "S"), each = J), rep(1:J, 3),
                               "(pi)", sep = ""),
                         paste("cum(S", 1:J, "(pi))", sep = ""), "max(N)")
  opchar[, 6 + 4*J] <- as.integer(opchar[, 6 + 4*J])
  return(tibble::as_tibble(opchar))
}

# Function to determine terminal states in a group sequential design
terminal_states_gs <- function(J, a, r, n, k) {
  max_k        <- max(k)
  terminal     <- NULL
  if (a[1] >= 0) {
    terminal   <- rbind(terminal, cbind(0:a[1], n[1], 1))
  }
  if (r[1] < Inf) {
    terminal   <- rbind(terminal, cbind(r[1]:n[1], n[1], 1))
  }
  cont         <- c(max(0, a[1] + 1), min(r[1] - 1, n[1]))
  if (max_k >= 2) {
    if (max_k >= 3) {
      for (j in 1:(max_k - 2)) {
        new      <- cbind(cont[1]:(cont[2] + n[j + 1]),
                          sum(n[1:(j + 1)]), j + 1)
        terminal <- rbind(terminal, new[which(new[, 1] <= a[j + 1] |
                                                new[, 1] >= r[j + 1]), ])
        cont     <- c(min(new[which(new[, 1] > a[j + 1]), 1]),
                      max(new[which(new[, 1] < r[j + 1]), 1]))
      }
    }
    new          <- cbind(cont[1]:(cont[2] + n[max_k]), sum(n[1:max_k]), max_k)
    if (max_k == J) {
      terminal   <- rbind(terminal, new)
    } else {
      terminal   <- rbind(terminal, new[which(new[, 1] <= a[max_k] |
                                                new[, 1] >= r[max_k]), ])
    }
  }
  terminal       <- terminal[which(terminal[, 3] %in% k), ]
  return(tibble::tibble(s = as.integer(terminal[, 1]),
                        m = as.integer(terminal[, 2]),
                        k = factor(terminal[, 3], k)))
}

# Function for use in determining ESS across a beta prior
ess_gs_beta <- function(pi, J, a, r, n, N, beta_prior) {
  return(int_opchar_gs(pi, J, a, r, n, N, 1:J)$ESS*
           dbeta(pi, beta_prior[1], beta_prior[2]))
}

# Function for finding UMVUE in a group sequential design
est_gs_umvue <- function(s, m, k, a, r, n) {
  if (k == 1) {
    return(s/m)
  } else if (k == 2) {
    s1 <- max(a[1] + 1, s - n[2], 0):min(s, r[1] - 1, n[1])
    return(sum((choose(n[1] - 1, s1 - 1)*choose(n[2], s - s1)))/
             sum((choose(n[1], s1)*choose(n[2], s - s1))))
  } else {
    R_ms        <- iterpc::getall(iterpc::iterpc(s + 1, k, 0:s, T, T))
    R_ms        <- R_ms[which(rowSums(R_ms) == s), ]
    for (j in 1:(k - 1)) {
      cum_R_ms  <- rowSums(R_ms[, 1:j, drop = F])
      R_ms      <- R_ms[which(cum_R_ms > a[j] & cum_R_ms < r[j]), , drop = F]
    }
    sum_num     <- 0
    sum_denom   <- 0
    for (i in 1:nrow(R_ms)) {
      sum_num   <- sum_num + prod(choose(n[1:k] - c(1, numeric(k - 1)),
                                  R_ms[i, ] - c(1, numeric(k - 1))))
      sum_denom <- sum_denom + prod(choose(n[1:k], R_ms[i, ]))
    }
    return(sum_num/sum_denom)
  }
}

# Function for finding UMVCUE in a group sequential design
est_gs_umvcue <- function(s, m, k, a, r, n) {
  if (k == 1) {
    return(s/m)
  } else {
    s1 <- max(a[1] + 1, s - n[2], 0):min(s, r[1] - 1, n[1])
    return(sum((choose(n[1], s1)*choose(n[2] - 1, s - s1 - 1)))/
                sum((choose(n[1], s1)*choose(n[2], s - s1))))
  }
}

# Function for finding bias-subtracted estimate in a group sequential design
est_gs_bias_sub <- function(s, m, J, a, r, n) {
  pmf <- pmf_gs(s/m, J, a, r, n, 1:J)
  return(2*s/m - sum(pmf$s*pmf$`f(s,m|pi)`/pmf$m))
}

# Function for finding bias-adjusted estimate in a group sequential design
est_gs_bias_adj <- function(s, m, J, a, r, n) {

  int_est_gs_bias_adj <- function(pi, J, a, r, n, pi_mle) {
    pmf <- pmf_gs(pi, J, a, r, n, 1:J)
    return(pi_mle - sum(pmf$s*pmf$`f(s,m|pi)`/pmf$m))
  }

  if (s == 0) {
    return(0)
  } else if (s == m) {
    return(1)
  } else {
    return(uniroot(int_est_gs_bias_adj, c(0, 1), J = J, a = a, r = r, n = n,
                   pi_mle = s/m)$root)
  }
}

# Function for finding conditional MLE in a group sequential design
est_gs_cond_mle <- function(s, m, k, a, r, n) {

  int_est_gs_cond_mle <- function(pi, s, m, k, a, r, n) {
    if (k == 2) {
      range <- max(0, a[1] + 1):min(r[1] - 1, n[1])
      return(-(sum(choose(n[1], range)*(pi^(range - 1))*
                     ((1 - pi)^(n[1] - range - 1))*(range - pi*n[1])))/
               sum(dbinom(range, n[1], pi)) + s/pi -
               (m - s)/(1 - pi))
    } else {
      continuation_k   <- iterpc::getall(iterpc::iterpc(max(n[1:(k - 1)]) + 1,
                                                        k - 1,
                                                        0:max(n[1:(k - 1)]),
                                                        T, T))
      for (j in 1:(k - 1)) {
        continuation_k <- continuation_k[which(continuation_k[, j] <= n[j]), ]
        row_sums       <- rowSums(continuation_k[, 1:j, drop = F])
        continuation_k <- continuation_k[which(row_sums > a[j] &
                                                 row_sums < r[j]), ]
      }
      product          <- 1
      for (j in 1:(k - 1)) {
        product        <- product*choose(n[j], continuation_k[, j])
      }
      row_sums         <- rowSums(continuation_k)
      G                <- sum(product*(pi^(row_sums))*
                                ((1 - pi)^(sum(n[1:(j - 1)]) - row_sums)))
      d_G              <- sum(product*(pi^(row_sums - 1))*
                                ((1 - pi)^(sum(n[1:(j - 1)]) - row_sums - 1))*
                                (row_sums - pi*sum(n[1:(j - 1)])))
      return(-d_G/G  + s/pi - (m - s)/(1 - pi))
    }
  }

  if (k == 1) {
    return(s/m)
  } else if (k == 2) {
    if (s <= max(0, a[1] + 1)) {
      return(0)
    } else if (s >= min(r[1] - 1, n[1]) + n[2]) {
      return(1)
    } else {
      return(uniroot(int_est_gs_cond_mle, c(10^-10, 1 - 10^-10), s = s, m = m,
                     k = 2, a = a, r = r, n = n)$root)
    }
  } else {
    if (s <= max(0, a[k] + 1)) {
      return(0)
    } else if (s >= min(r[k - 1] - 1, n[k - 1]) + n[k]) {
      return(1)
    } else {
      return(uniroot(int_est_gs_cond_mle, c(10^-10, 1 - 10^-10), s = s, m = m,
                     k = k, a = a, r = r, n = n)$root)
    }
  }
}

# Function for finding MUE in a group sequential design
est_gs_mue <- function(s, m, k, J, a, r, n) {
  if (s == 0) {
    return(0)
  } else if (s == m) {
    return(1)
  } else {
    return(uniroot(function(pi, s, m, k, J, a, r, n)
                     pval_gs_umvue(pi, s, m, k, J, a, r, n) - 0.5, c(0, 1),
                   s = s, m = m, k = as.numeric(k), J = J, a = a, r = r,
                   n = n)$root)
  }
}

# Function for finding p-value, based on MLE ordering, in a group sequential
# design
pval_gs_mle <- function(pi, s, m, J, a, r, n) {
  pmf <- pmf_gs(pi, J, a, r, n, 1:J)
  return(sum(pmf$`f(s,m|pi)`[pmf$s/pmf$m >= s/m]))
}

# Function for finding p-value, based on UMVUE ordering, in a group sequential
# design
pval_gs_umvue <- function(pi, s, m, k, J, a, r, n) {
  pmf               <- pmf_gs(pi, J, a, r, n, 1:J)
  umvues            <- tibble::tibble(s = pmf$s, m = pmf$m, k = pmf$k,
                                      `f(s,m|pi)` = pmf$`f(s,m|pi)`, umvue = NA)
  for (i in 1:nrow(umvues)) {
    umvues$umvue[i] <- est_gs_umvue(umvues$s[i], umvues$m[i],
                                    as.numeric(umvues$k[i]), a, r, n)
  }
  umvue_sm          <- est_gs_umvue(s, m, k, a, r, n)
  return(sum(dplyr::filter(umvues, umvue >= umvue_sm)$`f(s,m|pi)`))
}

# Function for finding conditional p-value in a group sequential design
pval_gs_cond <- function(pi, s, m, k, J, a, r, n) {
  if (k == 1) {
    return(pbinom(s - 1, n[1], pi, F))
  } else {
    pmf <- pmf_gs(pi, J, a, r, n, k)
    return(sum(pmf$`f(s,m|pi)`[which(pmf$s >= s)]))
  }
}

# Function for finding CI, using the exact method, in a group sequential design
ci_gs_exact <- function(s, m, k, J, a, r, n, alpha) {

  int_ci_gs_exact_clow <- function(pi, s, m, k, J, a, r, n, alpha) {
    return(pval_gs_umvue(pi, s, m, k, J, a, r, n) - alpha)
  }

  int_ci_gs_exact_cupp <- function(pi, s, m, k, J, a, r, n, alpha) {
    pmf                        <- pmf_gs(pi, J, a, r, n, 1:J)
    umvues                     <- dplyr::mutate(pmf, `pi(hat)(s,m)` = NA)
    for (i in 1:nrow(umvues)) {
      umvues$`pi(hat)(s,m)`[i] <-
        est_gs_umvue(umvues$s[i], umvues$m[i],
                     as.numeric(as.character(umvues$k[i])), a, r, n)
    }
    umvue_sm <- est_gs_umvue(s, m, k, a, r, n)
    return(sum(dplyr::filter(umvues, `pi(hat)(s,m)` <= umvue_sm)$`f(s,m|pi)`) -
             alpha)
  }

  if (s == 0) {
    Clow <- 0
  } else {
    Clow <- uniroot(int_ci_gs_exact_clow, c(0, 1), s = s, m = m, k = k, J = J,
                    a = a, r = r, n = n, alpha = alpha/2)$root
  }
  if (s == m) {
    Cupp <- 1
  } else {
    Cupp <- uniroot(int_ci_gs_exact_cupp, c(0, 1), s = s, m = m, k = k, J = J,
                    a = a, r = r, n = n, alpha = alpha/2)$root
  }
  return(c(Clow, Cupp))
}

# Function for finding CI, using the mid-p method, in a group sequential design
ci_gs_mid_p <- function(s, m, k, J, a, r, n, alpha) {

  ci_gs_mid_p_Clow <- function(pi, s, m, k, J, a, r, n, alpha) {
    pmf               <- pmf_gs(pi, J, a, r, n, 1:J)
    umvues            <- tibble::tibble(s = pmf$s, m = pmf$m, k = pmf$k,
                                        `f(s,m|pi)` = pmf$`f(s,m|pi)`,
                                        umvue = NA)
    for (i in 1:nrow(umvues)) {
      umvues$umvue[i] <- est_gs_umvue(umvues$s[i], umvues$m[i],
                                      as.numeric(umvues$k[i]), a, r, n)
    }
    umvue_sm          <- est_gs_umvue(s, m, k, a, r, n)
    prob_g_umvue      <- sum(dplyr::filter(umvues,
                                           umvue > umvue_sm)$`f(s,m|pi)`)
    prob_eq_umvue     <- sum(dplyr::filter(umvues,
                                           umvue == umvue_sm)$`f(s,m|pi)`)
    return(prob_g_umvue + 0.5*prob_eq_umvue - alpha/2)
  }

  ci_gs_mid_p_Cupp <- function(pi, s, m, k, J, a, r, n, alpha) {
    pmf               <- pmf_gs(pi, J, a, r, n, 1:J)
    umvues            <- tibble::tibble(s = pmf$s, m = pmf$m, k = pmf$k,
                                        `f(s,m|pi)` = pmf$`f(s,m|pi)`,
                                        umvue = NA)
    for (i in 1:nrow(umvues)) {
      umvues$umvue[i] <- est_gs_umvue(umvues$s[i], umvues$m[i],
                                      as.numeric(umvues$k[i]), a, r, n)
    }
    umvue_sm          <- est_gs_umvue(s, m, k, a, r, n)
    prob_l_umvue      <- sum(dplyr::filter(umvues,
                                           umvue < umvue_sm)$`f(s,m|pi)`)
    prob_eq_umvue     <- sum(dplyr::filter(umvues,
                                           umvue == umvue_sm)$`f(s,m|pi)`)
    return(prob_l_umvue + 0.5*prob_eq_umvue - alpha/2)
  }

  if (s == 0) {
    Clow <- 0
  } else {
    Clow <- uniroot(ci_gs_mid_p_Clow, c(0, 1), s = s, m = m, k = k, J = J,
                    a = a, r = r, n = n, alpha = alpha)$root
  }
  if (s == m) {
    Cupp <- 1
  } else {
    Cupp <- uniroot(ci_gs_mid_p_Cupp, c(0, 1), s = s, m = m, k = k, J = J,
                    a = a, r = r, n = n, alpha = alpha)$root
  }
  return(c(Clow, Cupp))
}
