#' Determine point estimates in a single-stage single-arm trial design for a
#' single binary endpoint
#'
#' Determines all possible point estimates (the MLEs) at the end of a
#' single-stage single-arm trial for a single binary endpoint, as determined
#' using \code{des_fixed()}.
#'
#' In addition, the performance of the point estimation procedure (including
#' its expected value, bias, and variance) for each value of
#' \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in the supplied vector
#' \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}}, will also be
#' evaluated.
#'
#' @param des An object of class \code{"sa_des_fixed"}, as returned by
#' \code{des_fixed()}.
#' @param pi A vector of response probabilities to evaluate the expected
#' performance of the point estimation procedure at. This will internally
#' default to be the \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}
#' and \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from
#' \code{des} if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_est_fixed"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$est} summarising the possible point
#' estimates at the end of the trial for the supplied design.
#' \item A tibble in the slot \code{$perf} summarising the performance of the
#' point estimation procedure for each specified value of
#' \ifelse{html}{\out{<i>&pi;</i>}}{\deqn{\pi}}.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal single-stage design for the default parameters
#' des <- des_fixed()
#' # Determine the performance of the point estimation procedure for a range of
#' # possible response probabilities
#' est <- opchar_fixed(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{pval_fixed}}, \code{\link{ci_fixed}}, and their associated
#' \code{plot} family of functions.
#' @export
est_fixed <- function(des, pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_fixed(des, "des")
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des$pi0, des$des$pi1)
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Point estimation for single-stage single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  pmf  <- pmf_fixed(pi, des$des$n)
  est  <- tibble::tibble(s = 0:des$des$n, m = rep(des$des$n, des$des$n + 1),
                         `hat(pi)(s,m)` = s/m)
  perf <- tibble::tibble(pi = pi, `E(hat(pi)|pi)` = pi, `Var(hat(pi)|pi)` = NA,
                         `Bias(hat(pi)|pi)` = 0)
  for (i in 1:nrow(perf)) {
    perf$`Var(hat(pi)|pi)`[i] <-
      sum(dplyr::filter(pmf,pi == perf$pi[i])$`f(s,m|pi)`*
            est$`hat(pi)(s,m)`^2) - perf$`E(hat(pi)|pi)`[i]^2
    if (all(summary, i%%100 == 0)) {
      message("...performance for ", i, " elements of pi evaluated...")
    }
  }
  perf <- dplyr::mutate(perf, `RMSE(hat(pi)|pi)`= sqrt(perf$`Var(hat(pi)|pi)`))

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(est = est, perf = perf, des = des, pi = pi,
                        summary = summary)
  class(output) <- "sa_est_fixed"
  return(output)
}
