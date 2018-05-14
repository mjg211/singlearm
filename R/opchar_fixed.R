#' Determine the operating characteristics of single-stage single-arm trial
#' designs for a single binary endpoint
#'
#' \code{opchar_fixed()} supports the simultaneous evaluation of the operating
#' characteristics of multiple single-stage single-arm clinical trial designs
#' for a single binary primary endpoint, determined using \code{des_fixed()}.
#'
#' Note that each of the supplied designs must have been designed for the same
#' value of \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}.
#'
#' For each value of \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} in
#' the supplied vector
#' \ifelse{html}{\out{<b><i>&pi;</i></b>}}{\eqn{\bold{\pi}}} (\code{pi}),
#' \code{opchar_fixed()} evaluates the power of each of the supplied designs.
#' Additional operating characteristics are also returned for convenience,
#' specifically to simplify combining the output with that from other
#' \code{opchar_###()} functions.
#'
#' @param des An object of class \code{"sa_des_fixed"}, as returned by
#' \code{des_fixed()}.
#' @param ... Additional objects of class \code{"sa_des_fixed"}. These will be
#' grouped in to a list named \code{"add_des"}.
#' @param pi A vector of response probabilities to evaluate operating
#' characteristics at. This will internally default to be
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and the
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the supplied
#' designs if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_opchar_fixed"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$opchar} summarising the operating
#' characteristics of the supplied designs. The power for each
#' \ifelse{html}{\out{<i>&pi;</i>}}{\deqn{\pi}} will be of principle use, but
#' other characteristics are also provided for convenience.
#' \item A tibble in the slot \code{$pmf} containing details of the PMF of each
#' of the supplied designs.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal single-stage design for the default parameters
#' des         <- des_fixed()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_fixed(des, pi = seq(0, 1, 0.01))
#' # Find the optimal single-stage design for a 10% type-I error rate
#' des_10      <- des_fixed(alpha = 0.1)
#' # Determine operating characteristics for both designs at the same time
#' opchar_both <- opchar_fixed(des, des_10,
#'                             pi = seq(0, 1, 0.01))
#' @references A'Hern RP (2001) Sample size tables for exact single-stage phase
#' II designs.  \emph{Statistics in Medicine} \strong{20:}859-66.
#' @references Fleming TR (1982) One-sample multiple testing procedure for phase
#' II clinical trials. \emph{Biometrics} \strong{38:}143-51.
#' @seealso \code{\link{des_fixed}}, \code{\link{est_fixed}},
#' \code{\link{pval_fixed}}, \code{\link{ci_fixed}}, and their associated
#' \code{plot} family of functions.
#' @export
opchar_fixed <- function(des, ..., pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_fixed(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_fixed(eval(add_des[[i]]),
                         paste("add_des[[", i, "]]", sep = ""))
    }
    for (i in 1:num_add_des) {
      if (eval(add_des[[i]])$des$pi0 != des$des$pi0) {
        stop("Each supplied design must have been designed for the same value of \u03c0\u2080")
      }
    }
  }
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi     <- c(des$des$pi0, des$des$pi1)
    if (num_add_des > 0) {
      for (i in 1:num_add_des) {
        pi <- c(pi, eval(add_des[[i]])$des$pi1)
      }
      pi   <- unique(pi)
    }
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Operating characteristic determination for single-stage single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  if (num_add_des == 0) {
    pmf     <- pmf_fixed(pi, des$des$n)
    opchar  <- int_opchar_fixed(pi, des$des$r, des$des$n, summary, pmf)
    add_des <- NULL
  } else {
    opchar             <- pmf <- list()
    pmf[[1]]           <- cbind("Design" = paste("Design 1: (", des$des$a, ", ",
                                                 des$des$r, ")/", des$des$n,
                                                 sep = ""),
                                pmf_fixed(pi, des$des$n))
    opchar[[1]]        <- cbind("Design" = paste("Design 1: (", des$des$a, ", ",
                                                 des$des$r, ")/", des$des$n,
                                                 sep = ""),
                                int_opchar_fixed(pi, des$des$r, des$des$n, F,
                                                 pmf[[1]]))
    if (summary) {
      message("...performance for Design 1 evaluated...")
    }
    for (i in 1:num_add_des) {
      des_i            <- eval(add_des[[i]])
      pmf[[i + 1]]     <- cbind("Design" = paste("Design ", i + 1, ": (",
                                                 des_i$des$a, ", ", des_i$des$r,
                                                 ")/", des_i$des$n, sep = ""),
                                pmf_fixed(pi, des_i$des$n))
      opchar[[i + 1]]  <- cbind("Design" = paste("Design ", i + 1, ": (",
                                                 des_i$des$a, ", ", des_i$des$r,
                                                 ")/", des_i$des$n, sep = ""),
                                int_opchar_fixed(pi, des_i$des$r, des_i$des$n,
                                                 F, pmf[[i + 1]]))
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
  class(output) <- "sa_opchar_fixed"
  return(output)
}
