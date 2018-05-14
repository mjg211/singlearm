#' Determine the operating characteristics of admissable group sequential
#' single-arm trial designs for a single binary endpoint
#'
#' \code{opchar_admissable()} supports the evaluation of the operating
#' characteristics of admissable group sequential single-arm clinical trial
#' designs for a single binary primary endpoint, determined using
#' \code{des_admissable()}.
#'
#' For each value of \ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in
#' the supplied vector \ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}},
#' \code{opchar_admissable()} evaluates the power, ESS, and other key
#' characteristics, of each of the supplied designs.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_admissable"}, as returned by
#' \code{des_admissable()}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and the maximal number of possible stages
#' in the supplied designs. If left unspecified, it will internally default to
#' all possible stages.
#' @param pi A vector of response probabilities to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
#' supplied designs if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' unction's progress should be printed to the console.
#' @return A list of class \code{"sa_opchar_admissable"} containing the
#' following elements
#' \itemize{
#' \item A tibble in the slot \code{$opchar} summarising the operating
#' characteristics of the supplied designs.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the admissable two-stage design for the default parameters
#' des         <- des_admissable()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_admissable(des, pi = seq(0, 1, 0.01))
#' @seealso \code{\link{des_admissable}}, and their associated \code{plot}
#' family of functions.
#' @export
opchar_admissable <- function(des, k, pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_admissable(des, "des")
  if (!missing(pi)) {
    check_pi(pi, "any")
  } else {
    pi <- c(des$des[[1]]$pi0, des$des[[1]]$pi1)
  }
  if (missing(k)) {
    k <- 1:des$des[[1]]$J
  }
  check_logical(summary, "summary")

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Operating characteristic determination for admissable group-sequential single-arm trials for a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k ∈ {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  num_des       <- length(des$des)
  opchar        <- pmf <- list()
  for (i in 1:num_des) {
    pmf[[i]]    <-
      cbind("Design" = paste("Design ", i, ": ",
                             paste("(", des$des[[i]]$a,
                                   ",", des$des[[i]]$r, ")/",
                                   cumsum(des$des[[i]]$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            pmf_gs(pi, des$des[[i]]$J, des$des[[i]]$a, des$des[[i]]$r,
                   des$des[[i]]$n, k))
    opchar[[i]] <-
      cbind("Design" = paste("Design ", i, ": ",
                             paste("(", des$des[[i]]$a,
                                   ",", des$des[[i]]$r, ")/",
                                   cumsum(des$des[[i]]$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            int_opchar_gs(pi, des$des[[i]]$J, des$des[[i]]$a, des$des[[i]]$r,
                          des$des[[i]]$n, cumsum(des$des[[i]]$n), k, F,
                          pmf[[i]]))
    if (summary) {
      message("...performance for Design ", i, " evaluated...")
    }
  }
  pmf           <- tibble::as_tibble(plyr::rbind.fill(pmf))
  pmf$m         <- as.integer(pmf$m)
  opchar        <- tibble::as_tibble(plyr::rbind.fill(opchar))
  opchar$Design <- as.factor(opchar$Design)

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(opchar = opchar, pmf = pmf, des = des, pi = pi,
                        summary = summary)
  class(output) <- "sa_opchar_admissable"
  return(output)
}
