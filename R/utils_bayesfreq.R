# Function for determining pmf of bayesian-frequentist design
pmf_bayesfreq <- function(mu, nu, J, a, r, n, k) {
  max_k           <- max(k)
  int_mu          <- mu
  int_nu          <- nu
  len_mu          <- length(int_mu)
  unique_n        <- unique(n[1:max_k])
  pmf_fixed_mu_nu <- tibble::tibble(mu = rep(int_mu,
                                             sum(unique_n) + length(unique_n)),
                                    nu = rep(int_nu,
                                             sum(unique_n) + length(unique_n)),
                                    s = rep(unlist(sapply(unique_n,
                                                          function(n) 0:n)),
                                            each = len_mu),
                                    m = rep(unlist(sapply(unique_n,
                                                          function(n)
                                                            rep(n, n + 1))),
                                            each = len_mu),
                                    f = choose(m, s)*
                                          beta(mu + s, nu + m - s)/
                                          beta(mu, nu))
  if (J == 1) {
    terminal         <- terminal_states_fixed(n)
  } else {
    terminal         <- terminal_states_gs(J, a, r, n, k)
  }
  rows_terminal    <- nrow(terminal)
  pmf              <- tibble::tibble(mu = rep(int_mu, each = rows_terminal),
                                     nu = rep(int_nu, each = rows_terminal),
                                     s = rep(terminal$s, len_mu),
                                     m = rep(terminal$m, len_mu),
                                     k = factor(rep(terminal$k, len_mu), k),
                                     `f(s,m|mu,nu)` = NA)
  if (J == 1) {
    for (i in 1:len_mu) {
      pmf$`f(s,m|mu,nu)`[which(pmf$mu == int_mu[i] & pmf$nu == int_nu[i])] <-
        dplyr::filter(pmf_fixed_mu_nu, mu == int_mu[i] & nu == int_nu[i])$f
    }
  } else {
    for (i in 1:len_mu) {
      pmf_mat                  <- matrix(0, nrow = sum(n[1:max_k]) + 1,
                                         ncol = max_k)
      pmf_mat[1:(n[1] + 1), 1] <- dplyr::filter(pmf_fixed_mu_nu,
                                                mu == int_mu[i] &
                                                  nu == int_nu[i] & m == n[1])$f
      cont                     <- c(max(0, a[1] + 1), min(r[1] - 1, n[1]))
      if (max_k == 2) {
        for (l in cont[1]:cont[2]) {
          pmf_mat[(l + 1):(l + 1 + n[2]), 2] <-
            pmf_mat[(l + 1):(l + 1 + n[2]), 2] +
            pmf_mat[l + 1, 1]*dplyr::filter(pmf_fixed_mu_nu, mu == int_mu[i] &
                                              nu == int_nu[i] & m == n[2])$f
        }
      }
      pmf_mat[(cont[1] + 1):(cont[2] + 1), 1] <- 0
      if (length(k) == 1) {
        pmf_mat[, -k] <- 0
        pmf_mat[, k]  <- pmf_mat[, k]/sum(pmf_mat[, k])
      }
      pmf$`f(s,m|mu,nu)`[which(pmf$mu == int_mu[i] &
                                 pmf$nu == int_nu[i])] <-
        pmf_mat[which(pmf_mat > 0)]
    }
  }
  return(pmf)
}

# Function for determining operating characteristics of bayesian-frequentist design
int_opchar_bayesfreq <- function(mu, nu, J, a, r, n, N, k, pmf_mu_nu) {
  if (missing(pmf_mu_nu)) {
    pmf_mu_nu      <- pmf_bayesfreq(mu, nu, J, a, r, n, k)
  }
  int_mu           <- mu
  int_nu           <- nu
  len_mu           <- length(int_mu)
  opchar           <- matrix(0, nrow = len_mu, ncol = 6 + 4*J)
  if (J == 1) {
    for (i in 1:len_mu) {
      R            <- sum(dplyr::filter(pmf_mu_nu, mu == int_mu[i] &
                                          nu == int_nu[i] &
                                          s >= r)$`f(s,m|mu,nu)`)
      opchar[i, ]  <- c(int_mu[i], int_nu[i], R, n, 0, n, 1 - R, R, 1, 1)
    }
  } else {
    for (i in 1:len_mu) {
      A <- R       <- numeric(J)
      for (j in k) {
        A[j]       <- sum(dplyr::filter(pmf_mu_nu, mu == int_mu[i] &
                                          nu == int_nu[i] & s <= a[j] &
                                          k == j)$`f(s,m|mu,nu)`)
        R[j]       <- sum(dplyr::filter(pmf_mu_nu, mu == int_mu[i] &
                                          nu == int_nu[i] & s >= r[j] &
                                          k == j)$`f(s,m|mu,nu)`)
      }
      cum_S        <- cumsum(S <- A + R)
      if (any(cum_S == 0.5)) {
        Med        <- 0.5*(N[which(cum_S == 0.5)] + N[which(cum_S == 0.5) + 1])
      } else {
        Med        <- N[which(cum_S > 0.5)[1]]
      }
      ESS          <- sum(N*S)
      opchar[i, ]  <- c(int_mu[i], int_nu[i], sum(R), ESS, sum(N^2*S) - ESS^2,
                        Med, A, R, S, cum_S)
    }
  }
  colnames(opchar) <- c("mu", "nu", "P(mu,nu)", "ESS(mu,nu)", "VSS(mu,nu)",
                        "Med(mu,nu)", paste(rep(c("A", "R", "S"), each = J),
                                            rep(1:J, 3), "(mu,nu)", sep = ""),
                        paste("cum(S", 1:J, "(mu,nu))", sep = ""))
  return(tibble::as_tibble(opchar))
}
