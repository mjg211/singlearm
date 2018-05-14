# Function for determining pmf of fixed design
pmf_fixed <- function(pi, n1) {
  len_pi <- length(pi)
  return(tibble::tibble(pi = rep(pi, each = n1 + 1), s = rep(0:n1, len_pi),
                        m = rep(n1, (n1 + 1)*len_pi),
                        `f(s,m|pi)` = stats::dbinom(s, m, pi)))
}

# Function for determining operating characteristics of fixed design
int_opchar_fixed <- function(pi, r1, n1, summary, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- pmf_fixed(pi, n1)
  }
  len_pi   <- length(int_pi <- pi)
  P        <- numeric(len_pi)
  for (i in 1:len_pi) {
    P[i]   <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & s >= r1)$`f(s,m|pi)`)
    if (all(summary, i%%1000 == 0)) {
      message("...performance for ", i, " elements of pi evaluated...")
    }
  }
  return(tibble::tibble(pi = int_pi, `P(pi)` = P, `ESS(pi)` = n1, `VSS(pi)` = 0,
                        `Med(pi)` = n1, `A1(pi)` = 1 - P, `R1(pi)` = P,
                        `S1(pi)` = 1, `cum(S1(pi))` = 1))
}

terminal_states_fixed <- function(n) {
  return(tibble::tibble(s = 0:n, m = rep(n, n + 1),
                        k = factor(rep(1, n + 1), 1)))
}

# Function for finding p-value, using the exact method, in a fixed design
pval_fixed_exact <- function(s, m, pi0) {
  return(stats::pbinom(s - 1, m, pi0, lower.tail = F))
}

# Function for finding p-value, using the normal method, in a fixed design
pval_fixed_normal <- function(s, m, pi0) {
  return(stats::pnorm(s/m, pi0, sqrt(pi0*(1 - pi0)/m), lower.tail = F))
}

# Function for determining Agresti-Coull CI in a fixed design
ci_fixed_agresti_coull <- function(s, m, alpha) {
  req_qnorm <- ifelse(s %in% c(0, m), stats::qnorm(alpha),
                      stats::qnorm(alpha/2))
  m_tilde   <- m + req_qnorm^2
  pi_tilde  <- (s + 0.5*req_qnorm^2)/m_tilde
  factor    <- req_qnorm*sqrt(pi_tilde*(1 - pi_tilde)/m_tilde)
  if (s == 0) {
    Clow    <- 0
  } else {
    Clow    <- pi_tilde + factor
    Clow    <- ifelse(Clow < 0, 0, ifelse(Clow > 1, 1, Clow))
  }
  if (s == m) {
    Cupp    <- 1
  } else {
    Cupp    <- pi_tilde - factor
    Cupp    <- ifelse(Cupp < 0, 0, ifelse(Cupp > 1, 1, Cupp))
  }
  return(c(Clow, Cupp))
}

# Function for determining Clopper-Pearson CI in a fixed design
ci_fixed_clopper_pearson <- function(s, m, alpha) {
  if (s == 0) {
    Clow <- 0
    Cupp <- 1 - (alpha/2)^(1/m)
  } else if (s == m) {
    Clow <- (alpha/2)^(1/m)
    Cupp <- 1
  } else {
    Clow <- stats::qbeta(alpha/2, s, m - s + 1)
    Cupp <- stats::qbeta(1 - alpha/2, s + 1, m - s)
  }
  return(c(Clow, Cupp))
}

# Function for determining Jeffreys CI in a fixed design
ci_fixed_jeffreys <- function(s, m, alpha) {
  if (s == 0) {
    Clow <- 0
    Cupp <- stats::qbeta(1 - alpha, s + 0.5, m - s + 0.5)
  } else if (s == m) {
    Clow <- stats::qbeta(alpha, s + 0.5, m - s + 0.5)
    Cupp <- 1
  } else {
    Clow <- stats::qbeta(alpha/2, s + 0.5, m - s + 0.5)
    Cupp <- stats::qbeta(1 - alpha/2, s + 0.5, m - s + 0.5)
  }
  return(c(Clow, Cupp))
}

# Function for determining mid-p CI in a fixed design
ci_fixed_mid_p <- function(s, m, alpha) {
  if (s == 0) {
    Clow <- 0
    Cupp <- 1 - alpha^(1/m)
  } else if (s == m) {
    Clow <- alpha^(1/m)
    Cupp <- 1
  } else {
    Clow <- stats::uniroot(function(pi, s, m, alpha)
                             0.5*stats::dbinom(s, m, pi) +
                             stats::pbinom(s, m, pi, F) - alpha/2,
                           c(0, s/m), s = s, m = m, alpha = alpha)$root
    Cupp <- stats::uniroot(function(pi, s, m, alpha)
                             0.5*stats::dbinom(s, m, pi) +
                             stats::pbinom(s - 1, m, pi) - alpha/2,
                    c(s/m, 1), s = s, m = m, alpha = alpha)$root
  }
  return(c(Clow, Cupp))
}

# Function for determining Wald CI in a fixed design
ci_fixed_wald <- function(s, m, alpha) {
  if (s == 0) {
    Clow   <- Cupp <- 0
  } else if (s == m) {
    Clow   <- Cupp <- 1
  } else {
    pi_hat <- s/m
    factor <- stats::qnorm(1 - alpha/2)*sqrt(pi_hat*(1 - pi_hat)/m)
    Clow   <- pi_hat - factor
    Cupp   <- pi_hat + factor
    Clow   <- ifelse(Clow < 0, 0, ifelse(Clow > 1, 1, Clow))
    Cupp   <- ifelse(Cupp < 0, 0, ifelse(Cupp > 1, 1, Cupp))
  }
  return(c(Clow, Cupp))
}

# Function for determining Wilson CI in a fixed design
ci_fixed_wilson <- function(s, m, alpha) {
  pi_hat <- s/m
  z      <- ifelse(s %in% c(0, m), stats::qnorm(1 - alpha),
                   -stats::qnorm(alpha/2))
  z2     <- z^2
  factor <- sqrt(pi_hat*(1 - pi_hat)/m + z2/(4*m^2))
  if (s != 0) {
    Clow <- (pi_hat + z2/(2*m) - z*factor)/(1 + z2/m)
    Clow <- ifelse(Clow < 0, 0, ifelse(Clow > 1, 1, Clow))
  } else {
    Clow <- 0
  }
  if (s != m) {
    Cupp <- (pi_hat + z2/(2*m) + z*factor)/(1 + z2/m)
    Cupp <- ifelse(Cupp < 0, 0, ifelse(Cupp > 1, 1, Cupp))
  } else {
    Cupp <- 1
  }
  return(c(Clow, Cupp))
}
