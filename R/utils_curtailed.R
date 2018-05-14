# Function to return conditional rejection probability in a group sequential
# design
cr_gs <- function(pi, s, m, J, a, r, n) {
  if (J == 1) {
    cr         <- stats::pbinom(r - s - 1, n - m, pi, lower.tail = F)
  } else {
    if (m <= n[1]) {
      if (s >= r[1]) {
        cr     <- 1
      } else {
        vals   <- unique(c(s:(s + n[1] - m), 0:n[2]))
        poss_s <- iterpc::getall(iterpc::iterpc(length(vals), 2, vals, T, T))
        poss_s <- poss_s[which(poss_s[, 1] %in% s:(s + n[1] - m) &
                                 poss_s[, 2] %in% 0:n[2]), ]
        poss_s <- poss_s[which(poss_s[, 1] >= r[1] |
                                 (poss_s[, 1] %in% max(0, a[1] + 1):
                                    min(n[1], r[1] - 1) &
                                    poss_s[, 1] + poss_s[, 2] >= r[2])), ]
        cr     <- sum(stats::dbinom(poss_s[, 1] - s, n[1] - m, pi)*
                        stats::dbinom(poss_s[, 2], n[2], pi))
      }
    } else {
      cr       <- stats::pbinom(r[2] - s - 1, sum(n) - m, pi, lower.tail = F)
    }
  }
  return(cr)
}

# Function for determining pmf of curtailed design
pmf_curtailed <- function(pi, J, J_curt, a_curt, r_curt, n, n_curt, k) {
  max_k        <- max(k)
  int_pi       <- pi
  len_pi       <- length(int_pi)
  unique_n     <- unique(n_curt[which(cumsum(n_curt) <= sum(n[1:max_k]))])
  pmf_fixed_pi <- tibble::tibble(pi = rep(int_pi,
                                          sum(unique_n) + length(unique_n)),
                                 s = rep(unlist(sapply(unique_n,
                                                       function(n) 0:n)),
                                         each = len_pi),
                                 m = rep(unlist(sapply(unique_n,
                                                       function(n) rep(n,
                                                                       n + 1))),
                                         each = len_pi),
                                 f = stats::dbinom(s, m, pi))
  N             <- c(0, cumsum(n))
  full_k        <- NULL
  for (i in 1:length(k)) {
    full_k      <- c(full_k, (N[k[i]] + 1):N[k[i] + 1])
  }
  terminal      <- terminal_states_gs(J_curt, a_curt, r_curt, n_curt, full_k)
  rows_terminal <- nrow(terminal)
  pmf           <- tibble::tibble(pi = rep(int_pi, each = rows_terminal),
                                  s = rep(terminal$s, len_pi),
                                  m = rep(terminal$m, len_pi),
                                  k = factor(rep(terminal$k, len_pi),
                                             unique(terminal$k)),
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
      pmf_mat                  <- matrix(0, nrow = sum(n) + 1, ncol = J_curt)
      pmf_mat[1:(n_curt[1] + 1), 1] <- dplyr::filter(pmf_fixed_pi,
                                                     pi == int_pi[i] &
                                                       m == n_curt[1])$f
      cont                     <- c(max(0, a_curt[1] + 1),
                                    min(r_curt[1] - 1, n_curt[1]))
      for (j in 1:(J_curt - 1)) {
        for (l in cont[1]:cont[2]) {
          pmf_mat[(l + 1):(l + 1 + n_curt[j + 1]), j + 1] <-
            pmf_mat[(l + 1):(l + 1 + n_curt[j + 1]), j + 1] +
            pmf_mat[l + 1, j]*dplyr::filter(pmf_fixed_pi, pi == int_pi[i] &
                                              m == n_curt[j + 1])$f
        }
        pmf_mat[(cont[1] + 1):(cont[2] + 1), j] <- 0
        new  <- cont[1]:(cont[2] + n_curt[j + 1])
        cont <- c(min(new[new > a_curt[j + 1]]), max(new[new < r_curt[j + 1]]))
      }
      if (length(full_k) < J_curt) {
        pmf_mat[, -full_k] <- 0
        pmf_mat[, full_k]  <- pmf_mat[, full_k]/sum(pmf_mat[, full_k])
      }
      f_i <- pmf_mat[which(pmf_mat > 0)]
    }
    pmf$`f(s,m|pi)`[which(pmf$pi == int_pi[i])] <- f_i
  }
  true_k      <- numeric(rows_terminal)
  for (i in 1:rows_terminal) {
    true_k[i] <- min(which(as.numeric(terminal$k[i]) <= N)) - 1
  }
  pmf$k       <- factor(rep(true_k, len_pi), unique(true_k))
  return(pmf)
}

# Function for determining operating characteristics of group sequential design
int_opchar_curtailed <- function(pi, J, J_curt, a_curt, r_curt, n, n_curt, N,
                                 N_curt, k, pmf_pi){
  if (missing(pmf_pi)) {
    pmf_pi <- pmf_curtailed(pi, J, J_curt, a_curt, r_curt, n, n_curt, k)
  }
  int_pi           <- pi
  len_pi           <- length(int_pi)
  opchar           <- matrix(0, nrow = len_pi, ncol = 5 + 4*J_curt + 4*J)
  for (i in 1:len_pi) {
    A <- R         <- numeric(J_curt)
    for (j in 1:J_curt) {
      A[j]         <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] &
                                          s <= a_curt[j] &
                                          m == N_curt[j])$`f(s,m|pi)`)
      R[j]         <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] &
                                          s >= r_curt[j] &
                                          m == N_curt[j])$`f(s,m|pi)`)
    }
    cum_S          <- cumsum(S <- A + R)
    if (any(cum_S == 0.5)){
      Med          <- 0.5*(N_curt[which(cum_S == 0.5)] +
                             N_curt[which(cum_S == 0.5) + 1])
    } else {
      Med          <- N_curt[which(cum_S > 0.5)[1]]
    }
    N              <- c(0, N)
    Atilde         <- Rtilde <- numeric(J)
    for (j in 1:J) {
      Atilde[j]    <- sum(A[which(N_curt %in% (N[j] + 1):N[j + 1])])
      Rtilde[j]    <- sum(R[which(N_curt %in% (N[j] + 1):N[j + 1])])
    }
    cum_Stilde     <- cumsum(Stilde <- Atilde + Rtilde)
    ESS            <- sum(N_curt*S)
    opchar[i, ]    <- c(int_pi[i], sum(R), ESS, sum(N_curt^2*S) - ESS^2,
                        Med, A, R, S, cum_S, Atilde, Rtilde, Stilde, cum_Stilde)
  }
  colnames(opchar) <- c("pi", "P(pi)", "ESS(pi)", "VSS(pi)", "Med(pi)",
                        paste(rep(c("A", "R", "S"), each = J_curt),
                              rep(1:J_curt, 3), "(pi)", sep = ""),
                        paste("cum(S", 1:J_curt, "(pi))", sep = ""),
                        paste(rep(c("Atilde", "Rtilde", "Stilde"),
                                  each = J),
                              rep(1:J, 3), "(pi)", sep = ""),
                        paste("cum(Stilde", 1:J, "(pi))", sep = ""))
  return(tibble::as_tibble(opchar))
}
