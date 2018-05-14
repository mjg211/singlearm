#' Determine the operating characteristics of adaptive two-stage single-arm
#' trial designs for a single binary endpoint
#'
#' \code{opchar_adaptive()} supports the simultaneous evaluation of the
#' operating characteristics of multiple adaptive two-stage single-arm clinical
#' trial designs for a single binary primary endpoint, determined using
#' \code{des_adaptive()}.
#'
#' Note that each of the supplied designs must have been designed for the same
#' value of \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}.
#'
#' For each value of \ifelse{html}{\out{<i>pi</i>)}}{\eqn{\pi}} in
#' the supplied vector \ifelse{html}{\out{<b><i>pi</i></b>)}}{\eqn{\bold{\pi}}},
#' \code{opchar_adaptive()} evaluates the power, ESS, and other key
#' characteristics, of each of the supplied designs.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_adaptive"}, as returned by
#' \code{des_adaptive()}.
#' @param ... Additional objects of class \code{"sa_des_adaptive"}. These will
#' be grouped in to a list named \code{"add_des"}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and two. If left unspecified, it will
#' internally default to all possible stages.
#' @param pi A vector of response probabilities to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
#' supplied designs if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_opchar_adaptive"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$opchar} summarising the operating
#' characteristics of the supplied designs.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal adaptive design for the default parameters
#' des    <- des_adaptive()
#' # Find its operating characteristics for a range of possible response
#' # probabilities
#' opchar <- opchar_adaptive(des)
#' @seealso \code{\link{des_adaptive}}, and their associated \code{plot} family
#' of functions.
#' @export
opchar_adaptive <- function(des, ..., k, pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_adaptive(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_adaptive(eval(add_des[[i]]), paste("add_des", i, sep = ""))
    }
    for (i in 1:num_add_des) {
      if (eval(add_des[[i]])$des$pi0 != des$des$pi0) {
        stop("Each supplied design must have been designed for the same value of pi0")
      }
    }
  }
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  if (missing(k)) {
    k <- 1:2
  } else if (!missing(k)){
    check_k(k, des, add_des)
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Operating characteristic determination for adaptive single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k \u2208 {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  pmf       <- list()
  if (num_add_des == 0) {
    pmf     <- pmf_adaptive(pi, des$des$a1, des$des$r1, des$des$n1, des$des$n2,
                            k)
    opchar  <- int_opchar_adaptive(pi, des$des$a1, des$des$r1, des$des$a2,
                                   des$des$r2, des$des$n1, des$des$n2, k, pmf)
    add_des <- NULL
  } else {
    opchar      <- list()
    pmf[[1]]    <- cbind("Design" = "Design 1",
                         pmf_adaptive(pi, des$des$a1, des$des$r1, des$des$n1,
                                      des$des$n2, k))
    opchar[[1]] <- cbind("Design" = "Design 1",
                         int_opchar_adaptive(pi, des$des$a1, des$des$r1,
                                             des$des$a2, des$des$r2, des$des$n1,
                                             des$des$n2, k, pmf[[1]]))
    if (summary) {
      message("...performance for Design 1 evaluated...")
    }
    for (i in 1:num_add_des) {
      des_i           <- eval(add_des[[i]])
      pmf[[i + 1]]    <- cbind("Design" = paste("Design", i + 1),
                               pmf_adaptive(pi, des_i$des$a1, des_i$des$r1,
                                            des_i$des$n1, des_i$des$n2, k))
      opchar[[i + 1]] <- cbind("Design" = paste("Design ", i + 1),
                               int_opchar_adaptive(pi, des_i$des$a1,
                                                   des_i$des$r1, des_i$des$a2,
                                                   des_i$des$r2, des_i$des$n1,
                                                   des_i$des$n2, k,
                                                   pmf[[i + 1]]))
      if (summary) {
        message("...performance for Design ", i + 1, " evaluated...")
      }
    }
    pmf                <- tibble::as_tibble(plyr::rbind.fill(pmf))
    pmf$m              <- as.integer(pmf$m)
    opchar             <- tibble::as_tibble(plyr::rbind.fill(opchar))
    opchar$Design      <- as.factor(opchar$Design)
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(opchar = opchar, pmf = pmf, des = des,
                        add_des = add_des, pi = pi, summary = summary)
  class(output) <- "sa_opchar_adaptive"
  return(output)
}
