#' Determine the operating characteristics of group sequential single-arm trial
#' designs for a single binary endpoint
#'
#' \code{opchar_gs()} supports the simultaneous evaluation of the operating
#' characteristics of multiple group sequential single-arm clinical trial
#' designs for a single binary primary endpoint, determined using
#' \code{des_gs()}.
#'
#' Note that each of the supplied designs must have been designed for the same
#' value of \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}.
#'
#' For each value of \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} in
#' the supplied vector \ifelse{html}{\out{<b><i>&pi;</i></b>}}{\eqn{\bold{\pi}}},
#' \code{opchar_gs()} evaluates the power, ESS, and other key characteristics,
#' of each of the supplied designs.
#'
#' Calculations are performed conditional on the trial stopping in one of the
#' stages specified using the input (vector) \code{k}.
#'
#' @param des An object of class \code{"sa_des_gs"}, as returned by
#' \code{des_gs()}.
#' @param ... Additional objects of class \code{"sa_gs_gs"}. These will be
#' grouped in to a list named \code{"add_des"}.
#' @param k Calculations are performed conditional on the trial stopping in one
#' of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
#' integers, with elements between one and the maximum number of possible stages
#' across the supplied designs. If left unspecified, it will internally default
#' to all possible stages.
#' @param pi A vector of response probabilities to evaluate operating
#' characteristics at. This will internally default to be the
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
#' supplied designs if it is left unspecified.
#' @param summary A logical variable indicating whether a summary of the
#' function's progress should be printed to the console.
#' @return A list of class \code{"sa_opchar_gs"} containing the following
#' elements
#' \itemize{
#' \item A tibble in the slot \code{$opchar} summarising the operating
#' characteristics of the supplied designs.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal two-stage design for the default parameters
#' des         <- des_gs()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_gs(des, pi = seq(0, 1, 0.01))
#' # Find the optimal two-stage design for a 10% type-I error rate
#' des_10      <- des_gs(alpha = 0.1)
#' # Determine operating characteristics for both designs at the same time
#' opchar_both <- opchar_gs(des, des_10, pi = seq(0, 1, 0.01))
#' @references Simon R (1989) Optimal two-stage designs for phase II clinical
#' trials. \emph{Controlled Clinical Trials} \strong{10:}1-10.
#' @seealso \code{\link{des_gs}}, \code{\link{est_gs}}, \code{\link{pval_gs}},
#' \code{\link{ci_gs}}, and their associated \code{plot} family of functions.
#' @export
opchar_gs <- function(des, ..., k, pi, summary = F) {

  ##### Input Checking #########################################################

  check_sa_des_gs(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_gs(eval(add_des[[i]]), paste("add_des", i, sep = ""))
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
    message("Operating Characteristic Determination for Group-Sequential Single-arm Trials with a Single Binary Endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("You have chosen to make your calculations conditional on k ∈ {", k[1], ",...,", k[length(k)], "}.\n")
    Sys.sleep(2)
    message("Beginning the required calculations...")
  }

  ##### Main Computations ######################################################

  if (num_add_des == 0) {
    pmf     <- pmf_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n, k)
    opchar  <- int_opchar_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n,
                             cumsum(des$des$n), k, summary, pmf)
    add_des <- NULL
  } else {
    opchar      <- pmf <- list()
    pmf[[1]]    <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            pmf_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n,
                   k[which(k <= Js[1])]))
    opchar[[1]] <-
      cbind("Design" = paste("Design 1: ",
                             paste("(", des$des$a, ",", des$des$r, ")/",
                                   cumsum(des$des$n), sep = "",
                                   collapse = ", "), sep = "", collapse = ""),
            int_opchar_gs(pi, des$des$J, des$des$a, des$des$r, des$des$n,
                          cumsum(des$des$n), k[which(k <= Js[1])], F, pmf[[1]]))
    if (summary) {
      message("...performance for Design 1 evaluated...")
    }
    for (i in 1:num_add_des) {
      des_i           <- eval(add_des[[i]])
      pmf[[i + 1]]    <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              pmf_gs(pi, des_i$des$J, des_i$des$a, des_i$des$r,
                     des_i$des$n, k[which(k <= Js[i + 1])]))
      opchar[[i + 1]] <-
        cbind("Design" = paste("Design ", i + 1, ": ",
                               paste("(", des_i$des$a,
                                     ",", des_i$des$r, ")/",
                                     cumsum(des_i$des$n), sep = "",
                                     collapse = ", "), sep = "", collapse = ""),
              int_opchar_gs(pi, des_i$des$J, des_i$des$a, des_i$des$r,
                            des_i$des$n, cumsum(des_i$des$n),
                            k[which(k <= Js[i + 1])], F, pmf[[i + 1]]))
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
  class(output) <- "sa_opchar_gs"
  return(output)

}
