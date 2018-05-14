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
    S1          <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] & k == 1)$f)
    for (j in 1:len_poss_n) {
      prob_n[j] <- sum(dplyr::filter(pmf_pi, pi == int_pi[i] &
                                       m == poss_n[j])$f)
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
