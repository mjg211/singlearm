#' Determine the operating characteristics of Bayesian-frequentist single-arm
#' trial designs for a single binary endpoint
#'
#' \code{opchar_bayesfreq()} supports the simultaneous evaluation of the
#' operating characteristics of multiple Bayesian-frequentist single-arm
#' clinical trial designs for a single binary primary endpoint, determined using
#' \code{des_bayesfreq()}.
#'
#' Note that each of the supplied designs must have been designed for the same
#' values of \ifelse{html}{\out{<i>&mu;</i>}}{\deqn{\mu}},
#' \ifelse{html}{\out{<i>&nu;</i>}}{\deqn{\nu}}, and
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}.
#'
#' For each value of \ifelse{html}{\out{<i>mu</i>}}{\eqn{\mu}},
#' \ifelse{html}{\out{<i>nu</i>}}{\eqn{\nu}}, and
#' \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in
#' the supplied vectors \ifelse{html}{\out{<b><i>mu</i></b>}}{\eqn{\bold{\mu}}},
#' \ifelse{html}{\out{<b><i>nu</i></b>}}{\eqn{\bold{\nu}}}, and
#' \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}},
#' \code{opchar_bayesfreq()} evaluates the Bayesian and frequentist power, ESS,
#' and other key characteristics, of each of the supplied designs.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_bayesfreq"}, as returned by
#' \code{des_bayesfreq()}.
#' @param ... Additional objects of class \code{"sa_des_bayesfreq"}. These will
#' be grouped in to a list named \code{"add_des"}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and the maximal number of possible stages
#' in the supplied designs. If left unspecified, it will internally default to
#' all possible stages.
#' @param mu A vector of the first Beta shape parameters to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&mu;</i>}}{\deqn{\mu}} from the supplied designs if it
#' left unspecified.
#' @param nu A vector of the second Beta shape parameters to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&nu;</i>}}{\deqn{\mu}} from the supplied designs if it
#' left unspecified.
#' @param pi A vector of response probabilities to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
#' supplied designs if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_opchar_bayesfreq"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$opchar_bayes} summarising the Bayesian
#' operating characteristics of the supplied designs.
#' \item A tibble in the slot \code{$opchar_freq} summarising the frequentist
#' operating characteristics of the supplied designs.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal two-stage Bayesian-frequentist design for the default
#' # parameters
#' des    <- des_bayesfreq()
#' # Determine operating characteristics for a range of mu, nu, and pi
#' opchar <- opchar_bayesfreq(des, mu = seq(0.05, 0.2, length.out = 10),
#'                            nu = seq(0.45, 1.8, length.out = 100),
#'                            pi = seq(0, 1, by = 0.01))
#' @seealso \code{\link{des_bayesfreq}}, and their associated \code{plot} family
#' of functions.
#' @export
opchar_bayesfreq <- function(des, ..., k, mu, nu, pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_bayesfreq(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_bayesfreq(eval(add_des[[i]]), paste("add_des", i, sep = ""))
    }
    for (i in 1:num_add_des) {
      if (eval(add_des[[i]])$des$mu != des$des$mu) {
        stop("Each supplied design must have been designed for the same value of mu")
      }
      if (eval(add_des[[i]])$des$nu != des$des$nu) {
        stop("Each supplied design must have been designed for the same value of nu")
      }
      if (eval(add_des[[i]])$des$pi0 != des$des$pi0) {
        stop("Each supplied design must have been designed for the same value of pi0")
      }
    }
  }
  if (!missing(mu)) {
    check_real_range_strict(mu, "mu", c(0, Inf), "any")
  } else {
    mu <- des$des$mu
  }
  if (!missing(nu)) {
    check_real_range_strict(nu, "nu", c(0, Inf), "any")
  } else {
    nu <- des$des$nu
  }
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  if (all(missing(k), num_add_des == 0)) {
    k           <- 1:des$des$J
  } else if (all(missing(k), num_add_des > 0)) {
    Js          <- numeric(num_add_des + 1)
    Js[1]       <- des$des$J
    for (i in 1:num_add_des) {
      Js[i + 1] <- eval(add_des[[i]])$des$J
    }
    k           <- 1:max(Js)
  } else if (!missing(k)){
    check_k(k, des, add_des)
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Operating characteristic determination for Bayesian-frequentist single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k \u2208 {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  mu_nu          <- as.matrix(expand.grid(mu = mu, nu = nu))
  opchar_bayes   <- opchar_freq <- pmf_bayes <- pmf_freq <- list()
  if (num_add_des == 0) {
    pmf_freq     <- pmf_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n, k)
    opchar_freq  <- int_opchar_gs(pi, des$des$J, des$des$a, des$des$r,
                                  des$des$n, cumsum(des$des$n), k, summary,
                                  pmf_freq)
    pmf_bayes    <- pmf_bayesfreq(mu_nu[, 1], mu_nu[, 2], des$des$J, des$des$a,
                                  des$des$r, des$des$n, k)
    opchar_bayes <- int_opchar_bayesfreq(mu_nu[, 1], mu_nu[, 2], des$des$J,
                                         des$des$a, des$des$r, des$des$n,
                                         cumsum(des$des$n), k, pmf_bayes)
    add_des      <- NULL
  } else {
    pmf_freq[[1]]    <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            pmf_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n,
                   k[which(k <= Js[1])]))
    opchar_freq[[1]] <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            int_opchar_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n,
                          cumsum(des$des$n), k[which(k <= Js[1])], F,
                          pmf_freq[[1]]))
    pmf_bayes[[1]]    <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            pmf_bayesfreq(mu_nu[, 1], mu_nu[, 2], des$des$J, des$des$a,
                          des$des$r, des$des$n, k[which(k <= Js[1])]))
    opchar_bayes[[1]] <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            int_opchar_bayesfreq(mu_nu[, 1], mu_nu[, 2], des$des$J, des$des$a,
                                 des$des$r, des$des$n, cumsum(des$des$n),
                                 k[which(k <= Js[1])], pmf_bayes[[1]]))
    if (summary) {
      message("...performance for Design 1 evaluated...")
    }
    for (i in 1:num_add_des) {
      des_i           <- eval(add_des[[i]])
      pmf_freq[[i + 1]]    <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              pmf_gs(pi, des_i$des$J, des_i$des$a, des_i$des$r,
                     des_i$des$n, k[which(k <= Js[i + 1])]))
      opchar_freq[[i + 1]] <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              int_opchar_gs(pi, des_i$des$J, des_i$des$a, des_i$des$r,
                            des_i$des$n, cumsum(des_i$des$n),
                            k[which(k <= Js[i + 1])], F, pmf_freq[[i + 1]]))
      pmf_bayes[[i + 1]]    <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              pmf_bayesfreq(mu_nu[, 1], mu_nu[, 2], des_i$des$J, des_i$des$a,
                            des_i$des$r, des_i$des$n, k[which(k <= Js[i + 1])]))
      opchar_bayes[[i + 1]] <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              int_opchar_bayesfreq(mu_nu[, 1], mu_nu[, 2], des_i$des$J,
                                   des_i$des$a, des_i$des$r, des_i$des$n,
                                   cumsum(des_i$des$n),
                                   k[which(k <= Js[i + 1])],
                                   pmf_bayes[[i + 1]]))
      if (summary) {
        message("...performance for Design ", i + 1, " evaluated...")
      }
    }
    pmf_freq           <- tibble::as_tibble(plyr::rbind.fill(pmf_freq))
    pmf_freq$m         <- as.integer(pmf_freq$m)
    opchar_freq        <- tibble::as_tibble(plyr::rbind.fill(opchar_freq))
    opchar_freq$Design <- as.factor(opchar_freq$Design)
    pmf_bayes           <- tibble::as_tibble(plyr::rbind.fill(pmf_bayes))
    pmf_bayes$m         <- as.integer(pmf_bayes$m)
    opchar_bayes        <- tibble::as_tibble(plyr::rbind.fill(opchar_bayes))
    opchar_bayes$Design <- as.factor(opchar_bayes$Design)
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(opchar_bayes = opchar_bayes, opchar_freq = opchar_freq,
                        pmf_bayes = pmf_bayes, pmf_freq = pmf_freq, des = des,
                        add_des = add_des, mu = mu, nu = nu, pi = pi,
                        summary = summary)
  class(output) <- "sa_opchar_bayesfreq"
  return(output)
}
