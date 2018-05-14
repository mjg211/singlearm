opchar_deviation <- function(des, pmf_n1n2, k, pi, summary = F) {

  # need to check J=2
  pi0 <- des$pi0
  a   <- des$des$a
  r   <- des$des$r
  n   <- des$des$n

  D                          <- numeric(n[1] + 1)
  if (r[1] < Inf) {
    D[(r[1] + 1):(n[1] + 1)] <- 1
  }
  for (s in max(0, a[1] + 1):min(r[1] - 1, n[1])) {
    D[s + 1]                 <- pbinom(r[2] - s - 1, n[2], pi0, lower.tail = F)
  }


  poss_n1 <- unique(pmf_n1n2$n1)
  D_mat   <- matrix(0, nrow = max(poss_n1), ncol = max(poss_n1) + 1)
  for (i in 1:length(poss_n1)) {
    if (poss_n1[i] < n[1]) {
      D_i <- numeric(poss_n1[i] + 1)
      for (s1 in 0:poss_n1[i]) {
        range       <- s1:(s1 + n[1] - poss_n1[i])
        D_i[s1 + 1] <- sum(D[range + 1]*dbinom(range - s1,
                                               n[1] - poss_n1[i], pi0))
      }
      D_mat[poss_n1[i], 1:(poss_n1[i] + 1)] <- D_i
    } else if (poss_n1[i] == n[1]) {
      D_mat[poss_n1[i], 1:(poss_n1[i] + 1)] <- D
    } else {
      n1star <- poss_n1[i]
      n1 <- n[1]

      coef_mat         <- matrix(0, n1 + 1, n1star + 1)
      dens             <- dbinom(0:(n1star - n1), n1star - n1, pi0)
      for (s1 in 1:(n1 + 1)) {
        coef_mat[s1, s1 + 0:(n1star - n1)] <- dens
      }
      which_0 <- which(D == 0)
      which_1 <- which(D == 1)
      if (length(which_0) > 0) {
        num0  <- tail(which(coef_mat[which_0[length(which_0)], ] > 0), n = 1)
      } else {
        num0  <- 0
      }
      if (length(which_1) > 0) {
        num1      <- n1star + 1 - tail(which(coef_mat[which_1[1] - 1, ] > 0),
                                       n = 1)
        check     <- T
        row  <- which_1[1] - 1
        col <- tail(which(coef_mat[which_1[1] - 1, ] > 0), n = 1)
        while (check) {
          if (sum(coef_mat[row, col:(n1star + 1)]) <= D[row]) {
            num1 <- num1 + 1
            col  <- col - 1
          } else {
            check <- F
          }
        }
      } else {
        num1  <- 0
      }
      effn1star      <- n1star - num0 - num1
      if (effn1star == 0) {
        diag_effn1star <- diag(effn1star + 1)
        ui <- coef_mat[which(D > 0 & D < 1), (num0 + 1):(n1star + 1 - num1), drop = F]
        ci <- D[which(D > 0 & D < 1)]
        D_mat[poss_n1[i], 1:(poss_n1[i] + 1)] <- c(rep(0, num0),
                                                   min(1, min(ci/ui)),
                                                   rep(1, num1))
      } else if (effn1star < 0) {
        D_mat[poss_n1[i], 1:(poss_n1[i] + 1)] <- c(rep(0, num0),
                                                   rep(1, n1star + 1 - num0))
      } else {
        coef_mat <- coef_mat[which(D > 0 & D < 1), (num0 + 1):(n1star + 1 - num1), drop = F]
        diag_effn1star <- diag(effn1star + 1)
        if (effn1star > 0) {
          ineq           <- -diag_effn1star[-(effn1star + 1), , drop = F]
          for (j in 1:effn1star) {
            ineq[j, j + 1] <- 1
          }
        } else {
          ineq <- NULL
        }
        ui <- -rbind(diag_effn1star, -diag_effn1star, ineq, -coef_mat)
        ci <- -c(rep(0, effn1star + 1), rep(-1, effn1star + 1), rep(0, effn1star),
                 -D[which(D > 0 & D < 1)])

        optimal_i <- CEoptim::CEoptim(f = max_Dstar, maximize = T,
                             f.arg = list(n1star = n1star, pi0 = pi0,
                                          begin = num0, end = n1star - num1),
                             continuous = list(mean = seq(from = 1e-2, to = 1e-1,
                                                          length.out = effn1star + 1),
                                               sd = rep(1, effn1star + 1),
                                               conMat = ui,
                                               conVec = ci,
                                               smoothMean = 0.5,
                                               smoothSd = 0.5),
                             N = as.integer(1000*(effn1star + 1)),
                             rho = 0.01,
                             verbose = T, noImproveThr = 20)
        D_mat[poss_n1[i], 1:(poss_n1[i] + 1)] <- c(rep(0, num0),
                                                   optimal_i$optimizer$continuous,
                                                   rep(1, num1))
      }
    }
  }

  #pi, s1, s2, n1, n2, a1, a2, r1, r2, D, s, m, k, f

  states <- list()
  for (i in 1:nrow(sample_sizes)) {
    n1 <- sample_sizes[i, 1]
    n2 <- sample_sizes[i, 2]
    poss_s1s2 <- getall(iterpc(n = max(n1, n2) + 2, r = 2,
                               labels = c(NA_real_, 0:max(n1, n2)),
                               ordered = T, replace = T))
    s1 <- which(D_mat[n1, 1:(n1 + 1)] > 0)[1]:tail(which(D_mat[n1, 1:(n1 + 1)] < 1), n = 1) - 1
    s2 <- 0:n2
    poss_s1s2 <- poss_s1s2[which((poss_s1s2[, 1] %in% s1 & !is.na(poss_s1s2[, 2])) |
                               (!(poss_s1s2[, 1] %in% s1) & !is.na(poss_s1s2[, 1]) &
                                  is.na(poss_s1s2[, 2]) & poss_s1s2[, 1] <= n1)), ]
    D  <- D_mat[n1, poss_s1s2[, 1] + 1]
    if (any(D_mat[n1, 1:(n1 + 1)] == 0)) {
      a1 <- tail(which(D_mat[n1, 1:(n1 + 1)] == 0), n = 1) - 1
    }
    if (any(D_mat[n1, 1:(n1 + 1)] == 1)) {
      r1 <- which(D_mat[n1, 1:(n1 + 1)] == 1)[1] - 1
    }
    r2 <- qbinom(D, n2, pi0, lower.tail = F) + 1 + poss_s1s2[, 1]
    r2[which(is.na(poss_s1s2[, 2]))] <- NA_real_
    s <- rowSums(poss_s1s2)
    m <- rep(n1 + n2, nrow(poss_s1s2))
    s[which(is.na(s))] <- poss_s1s2[which(is.na(s)), 1]
    m[which(is.na(poss_s1s2[, 2]))] <- n1
    states[[i]] <- tibble::tibble(s1 = as.integer(poss_s1s2[, 1]),
                                  s2 = as.integer(poss_s1s2[, 2]),
                                  n1 = as.integer(n1), n2 = as.integer(n2),
                                  a1 = as.integer(a1), a2 = as.integer(r2 - 1),
                                  r1 = as.integer(r1), r2 = as.integer(r2),
                                  D = D, s = as.integer(s), m = as.integer(m),
                                  k = as.integer(1 + !is.na(poss_s1s2[, 2])))
  }

  pmf <- list()
  counter <- 1
  for (i in 1:nrow(sample_sizes)) {
    for (j in 1:length(pi)) {
      f <- ifelse(states[[i]]$k == 1,
                  dbinom(states[[i]]$s1, states[[i]]$n1, pi0),
                  dbinom(states[[i]]$s1, states[[i]]$n1, pi0)*
                    dbinom(states[[i]]$s2, states[[i]]$n2, pi0))
      pmf[[counter]] <- tibble::as_tibble(cbind(pi[j], states[[i]], f))
      counter        <- counter + 1
    }
  }
  pmf <- plyr::rbind.fill(pmf)


  max_Dstar <- function(Dstar, n1star, pi0, begin, end) {
    return(sum(Dstar*dbinom(begin:end, n1star, pi0)))
  }


}
