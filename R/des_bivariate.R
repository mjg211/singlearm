des_bivariate <- function() {
  #BD
  Nmin <- 1
  Nmax <- 60
  piR0  <- 0.05
  piR1  <- 0.25
  piT0  <- 0.4
  piT1  <- 0.2
  theta <- 1
  dbivar_00_list     <- list()
  dbivar_01_list     <- list()
  dbivar_10_list     <- list()
  dbivar_11_min_list <- list()
  dbivar_11_max_list <- list()
  for (n in Nmin:Nmax) {
    if (is.null(theta)) {
      dbivar_00_list[[n]]     <- dbivariate(n, piR0, piT0, Inf)
      dbivar_01_list[[n]]     <- dbivariate(n, piR0, piT1, Inf)
      dbivar_10_list[[n]]     <- dbivariate(n, piR1, piT0, Inf)
      dbivar_11_min_list[[n]] <- dbivariate(n, piR1, piT1, Inf)
      dbivar_11_max_list[[n]] <- dbivariate(n, piR1, piT1, 0)
    } else {
      dbivar_00_list[[n]]     <- dbivariate(n, piR0, piT0, theta)
      dbivar_01_list[[n]]     <- dbivariate(n, piR0, piT1, theta)
      dbivar_10_list[[n]]     <- dbivariate(n, piR1, piT0, theta)
      dbivar_11_min_list[[n]] <- dbivariate(n, piR1, piT1, theta)
      dbivar_11_max_list[[n]] <- dbivar_11_min_list[[n]]
    }
  }

  #CP
  Nmin <- 1
  Nmax <- 60
  piR0  <- 0.5
  piR1  <- 0.75
  piT0  <- 0.3
  piT1  <- 0.15
  theta <- 2
  dbivar_00_list     <- list()
  dbivar_01_list     <- list()
  dbivar_10_list     <- list()
  dbivar_11_list <- list()
  dbivar_max1_list <- list()
  dbivar_max2_list <- list()
  for (n in Nmin:Nmax) {
    dbivar_00_list[[n]]     <- dbivariate(n, piR0, piT0, theta)
    dbivar_01_list[[n]]     <- dbivariate(n, piR0, piT1, theta)
    dbivar_10_list[[n]]     <- dbivariate(n, piR1, piT0, theta)
    dbivar_11_list[[n]]     <- dbivariate(n, piR1, piT1, theta)
    dbivar_max1_list[[n]]   <- cbind(stats::dbinom(0:n, n, piR0),
                                     matrix(0, nrow = n + 1, ncol = n))
    dbivar_max2_list[[n]]   <- rbind(matrix(0, nrow = n, ncol = n + 1),
                                     stats::dbinom(0:n, n, piT0))
  }

}
