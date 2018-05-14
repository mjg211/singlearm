# Function for determining the pmf of a variable sample size design
pmf_variable <- function(pi, J, v, a, n, dbinom_pi, dsamsize) {
  if (missing(dbinom_pi)) {
    dbinom_pi                 <- matrix(0, nrow = max(n) + 1, ncol = max(n) + 1)
    unique_n                  <- unique(as.vector(n))
    for (i in unique_n) {
      dbinom_pi[1:(i + 1), i] <- dbinom(0:i, i, pi)
    }
  }
  if (missing(dsamsize)) {
    if (J == 1) {
      dsamsize   <- cbind(n[1, ], 1/v)
    } else {
      dsamsize   <- cbind(expand.grid(n[1, ], n[2, ]), 1/v^2)
    }
  }
  pdf_mat        <- list()
  if (J == 1) {
    pdf_mat[[1]] <- matrix(0, nrow = max(n[1, ]) + 1, ncol = max(n[1, ]) + 1)
    for (i in 1:v) {
      pdf_mat[[1]][1:(n[1, i] + 1), n[1, i]] <- dbinom_pi[1:(n[1, i] + 1),
                                                          n[1, i]]*
                                                  dsamsize[which(dsamsize[, 1] ==
                                                                   n[1, i]), 2]
    }
    terminal <- terminal_states_variable(J, v, a, n)
    f        <- numeric(nrow(terminal))
    for (i in 1:v){
      for (s in 0:n[1, i]) {
        f[which(terminal$s == s & terminal$m == n[1, i])] <-
          pdf_mat[[1]][s + 1, n[1, i] + 1]
      }
    }
  } else {
    test <- 0
    for (i in 1:v) {
      pdf_mat[[i]]        <- list()
      for (j in 1:v) {
        pdf_mat[[i]][[j]] <- matrix(0, nrow = n[1, i] + n[2, j] + 1, ncol = 2)
        pdf_mat[[i]][[j]][1:(n[1, i] + 1), 1] <- dbinom_pi[1:(n[1, i] + 1),
                                                           n[1, i]]
        for (s in (a[1, i] + 1):n[1, i]) {
          pdf_mat[[i]][[j]][(s + 1):(s + 1 + n[2, j]), 2] <-
            pdf_mat[[i]][[j]][(s + 1):(s + 1 + n[2, j]), 2] +
              pdf_mat[[i]][[j]][s + 1, 1]*dbinom_pi[1:(n[2, j] + 1), n[2, j]]
        }
        pdf_mat[[i]][[j]][(a[1, i] + 2):(n[1, i] + 1), 1] <- 0
        pdf_mat[[i]][[j]] <- dsamsize[which(dsamsize[, 1] == n[1, i] &
                                                   dsamsize[, 2] == n[2, j]), 3]*
                                    pdf_mat[[i]][[j]]
        test              <- test + sum(pdf_mat[[i]][[j]][, 1])
      }
    }
    terminal <- terminal_states_variable(J, v, a, n)
    f        <- numeric(nrow(terminal))
    for (i in 1:v) {
      for (j in 1:v) {
        for (s1 in 0:a[1, i]) {
          f[which(terminal$s == s1 & terminal$m == n[1, i])] <-
            f[which(terminal$s == s1 & terminal$m == n[1, i])] +
              pdf_mat[[i]][[j]][s1 + 1, 1]
        }
        for (s in (a[1, i] + 1):(n[1, i] + n[2, j])) {
            f[which(terminal$s == s & terminal$m == n[1, i] + n[2, j])] <-
              f[which(terminal$s == s & terminal$m == n[1, i] + n[2, j])] +
                pdf_mat[[i]][[j]][s + 1, 2]
        }
      }
    }
  }
  return(tibble::tibble(s = terminal$s, m = terminal$m, k = terminal$k,
                        f = f))
}

# Function for determining terminal states in a variable sample size design
terminal_states_variable <- function(J, v, a, n) {
  if (J == 1) {
    terminal <- expand.grid(0:max(n), n, 1)
    terminal <- terminal[which(terminal[, 1] <= terminal[, 2]), ]
  } else {
    terminal <- NULL
    for (i in 1:v) {
      for (j in 1:v) {
        terminal <- rbind(terminal, cbind(0:a[1, i], n[1, i], 1),
                          cbind((a[1, i] + 1):(n[1, i] + n[2, j]),
                                n[1, i] + n[2, j], 2))
      }
    }
    terminal <- terminal[!duplicated(terminal), ]
  }
  return(tibble::tibble(s = terminal[, 1], m = terminal[, 2],
                        k = factor(terminal[, 3], 1:J)))
}
