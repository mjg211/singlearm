#' Design a group sequential single-arm trial for a single binary endpoint
#'
#' Determines group sequential single-arm clinical trial designs for a single
#' binary primary endpoint. In particular, this allows Simon's two-stage designs
#' (Simon, 1989) to be identified.
#'
#' \code{des_gs()} supports the determination of a variety of (optimised)
#' group sequential single-arm clinical trial designs for a single binary
#' primary endpoint. For all supported designs, the following hypotheses are
#' tested for the response probability
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}
#'
#' \ifelse{html}{\out{<center><i>H</i><sub>0</sub> : <i>&pi;</i> &le; <i>&pi;
#' </i><sub>0</sub>, <i>H</i><sub>1</sub> : <i>&pi;</i> > <i>&pi;</i><sub>
#' 0</sub>,</center>}}{\deqn{H_0 : \pi \le \pi_0,\qquad H_1 : \pi > \pi_0,}}
#'
#' for \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\eqn{\pi_0}},  specified
#' using the argument \code{pi0}.
#'
#' In each instance, the optimal design is required to meet the following
#' operating characteristics
#'
#' \ifelse{html}{\out{<center><i>P</i>(<i>&pi;</i><sub>0</sub>) &le;
#' <i>&alpha;</i>, <i>P</i>(<i>&pi;</i><sub>1</sub>) &ge; 1 - <i>&beta;</i>,
#' </center>}}{\deqn{P(\pi_0) \le \alpha,\qquad P(\pi_1) \ge 1 - \beta,}}
#'
#' where \ifelse{html}{\out{<i>P</i>(<i>&pi;</i>)}}{\eqn{P(\pi)}} is the
#' probability of rejecting \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}
#' when the true response probability is
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}, and the values of
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}} and
#' \ifelse{html}{\out{<i>&beta;</i>}}{\eqn{\beta}} are specified using the
#' arguments \code{alpha} and \code{beta} respectively. Moreover,
#' \ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\eqn{\pi_1}}, satisfying
#' \ifelse{html}{\out{<i>&pi;</i><sub>0</sub> &lt;
#' <i>&pi;</i><sub>1</sub>}}{\eqn{\pi_0 < \pi_1}}, is specified using the
#' argument \code{pi1}.
#'
#' A group sequential single-arm design for a single binary endpoint, with a
#' maximum of \ifelse{html}{\out{<i>J</i>}}{\eqn{J}} allowed stages (specifying
#' \ifelse{html}{\out{<i>J</i>}}{\eqn{J}} through the argument \code{J}) is
#' then indexed by three vectors: \ifelse{html}{\out{<b><i>a</i></b> =
#' (<i>a</i><sub>1 </sub>,&hellip;,<i>a</i><sub><i>J</i></sub>)}}{\eqn{\bold{a}=
#' (a_1,\dots,a_J)}}, \ifelse{html}{\out{<b><i>r</i></b> = (<i>r</i><sub>1
#' </sub>,&hellip;,<i>r</i><sub><i>J</i></sub>)}}{\eqn{\bold{r}=
#' (r_1,\dots,r_J)}}, and \ifelse{html}{\out{<b><i>n</i></b> = (<i>n</i><sub>1
#' </sub>,&hellip;,<i>n</i><sub><i>J</i></sub>)}}{\eqn{\bold{n}=
#' (n_1,\dots,n_J)}}.
#'
#' With these vectors, and denoting the number of responses after
#' \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} patients have been observed by
#' \ifelse{html}{\out{<i>s</i><sub><i>m</i></sub>}}{\eqn{s_m}}, the stopping
#' rules for the trial are then as follows
#'
#' \itemize{
#' \item For \ifelse{html}{\out{<i>j</i> = 1,&hellip;,<i>J - 1</i>}}{\eqn{j=
#' 1\dots,J-1}}
#' \itemize{
#' \item If \ifelse{html}{\out{<i>s</i><sub><i>N</i><sub><i>j</i></sub></sub>
#' &le; <i>a</i><sub><i>j</i></sub>}}{\eqn{s_{N_j} \le a_j}}, then stop the
#' trial and do not reject \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' \item Else if \ifelse{html}{\out{<i>s</i><sub><i>N</i><sub><i>j</i></sub>
#' </sub> &ge; <i>r</i><sub><i>j</i></sub>}}{\eqn{s_{N_j} \ge r_j}}, then stop
#' the trial and reject \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' \item Else if \ifelse{html}{\out{<i>a</i><sub><i>j</i></sub> &lt;
#' <i>s</i><sub><i>N</i><sub><i>j</i></sub></sub> &lt; <i>r</i><sub><i>j</i>
#' </sub>}}{\eqn{a_j < s_{N_j} < r_j}}, then continute to stage
#' \ifelse{html}{\out{<i>j</i> + 1}}{\eqn{j+1}}.
#' }
#' \item For \ifelse{html}{\out{<i>j</i> = <i>J</i>}}{\eqn{j=J}}
#' \itemize{
#' \item If \ifelse{html}{\out{<i>s</i><sub><i>N</i><sub><i>j</i></sub></sub>
#' &le; <i>a</i><sub><i>j</i></sub>}}{\eqn{s_{N_j} \le a_j}}, then
#' do not reject \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' \item Else if \ifelse{html}{\out{<i>s</i><sub><i>N</i><sub><i>j</i></sub>
#' </sub> &ge; <i>r</i><sub><i>j</i></sub>}}{\eqn{s_{N_j} \ge r_j}}, then reject
#' \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}.
#' }
#' }
#' Here, \ifelse{html}{\out{<i>N</i><sub><i>j</i></sub> = <i>n</i><sub><i>1</i>
#' </sub> + &ctdot; + <i>n</i><sub><i>j</i></sub>}}{\eqn{N_j = n_1 + \dots
#' n_j}}.
#'
#' The purpose of this function is then to optimise \ifelse{html}{\out{<b><i>a
#' </i></b>}}{\eqn{\bold{a}}}, \ifelse{html}{\out{<b><i>r</i>
#' </b>}}{\eqn{\bold{r}}}, and \ifelse{html}{\out{<b><i>n</i>
#' </b>}}{\eqn{\bold{n}}}, accounting for the chosen restrictions placed on
#' these vectors, and the chosen optimality criteria.
#'
#' The arguments \code{Nmin}, \code{Nmax}, and \code{equal_n} allow restrictions
#' to be placed on \ifelse{html}{\out{<b><i>n</i></b>}}{\eqn{\bold{n}}}.
#' Precisely, \code{Nmin} and \code{Nmax} set an inclusive range of allowed
#' values for \ifelse{html}{\out{<i>N</i><sub><i>J</i></sub>}}{\eqn{N_J}}.
#' While, if set to \code{TRUE}, \code{equal_n} enforces
#' \ifelse{html}{\out{<i>n</i><sub>1</sub> = &ctdot; = <i>n</i><sub><i>J</i>
#' </sub>}}{\eqn{n_1 = \dots = n_J}}.
#'
#' The arguments \code{futility}, \code{efficacy}, and \code{ensign}
#' allow restrictions to be placed on
#' \ifelse{html}{\out{<b><i>a</i></b>}}{\eqn{\bold{a}}} and
#' \ifelse{html}{\out{<b><i>r</i></b>}}{\eqn{\bold{r}}}. If \code{futility} is
#' set to \code{FALSE}, early stopping for futility (to not reject
#' \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}) is prevented by
#' enforcing \ifelse{html}{\out{<i>a</i><sub>1</sub> = &ctdot; = <i>a</i><sub>
#' <i>J</i> - 1</sub> = -&infin;}}{\eqn{a_1 = \dots = a_{J-1} = -\infty}}.
#' Similarly, if \code{efficacy} is set to \code{FALSE}, early stopping for
#' efficacy (to reject \ifelse{html}{\out{<i>H</i><sub>0</sub>}}{\eqn{H_0}}) is
#' prevented by enforcing \ifelse{html}{\out{<i>r</i><sub>1</sub> = &ctdot; =
#' <i>r</i><sub><i>J</i> - 1</sub> = &infin;}}{\eqn{r_1 = \dots = r_{J-1} =
#' \infty}}. Finally, if set to \code{TRUE}, \code{ensign} enforces the
#' restriction that \ifelse{html}{\out{<i>a</i><sub>1</sub> = 0}}{\eqn{a_1 =
#' 0}}, as suggested in Ensign \emph{et al} (1994) for 3-stage designs.
#'
#' Note that to ensure a decision is made about \ifelse{html}{\out{<i>H</i>
#' <sub>0</sub>}}{\eqn{H_0}}, this function enforces \ifelse{html}{\out{<i>a</i>
#' <sub><i>J</i></sub> + 1 = <i>r</i><sub><i>J</i>}}{\eqn{a_J + 1 = r_J}}.
#'
#' To describe the supported optimality criteria, denote the expected sample
#' size and median required sample size when the true response probability is
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} by
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i>)}}{\eqn{ESS(\pi)}} and
#' \ifelse{html}{\out{<i>Med</i>(<i>&pi;</i>)}}{\eqn{Med(\pi)}} respectively.
#' Then, the following optimality criteria are currently supported:
#' \itemize{
#' \item \code{"minimax"}: The design which minimises
#' \ifelse{html}{\out{<i>N</i><sub><i>J</i></sub>}}{\eqn{N_J}}.
#' \item \code{"null_ess"}: The design which minimises
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i><sub>0</sub>)}}{\eqn{ESS(\pi_0)}}.
#' \item \code{"alt_ess"}: The design which minimises
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i><sub>1</sub>)}}{\eqn{ESS(\pi_1)}}.
#' \item \code{"null_med"}: The design which minimises
#' \ifelse{html}{\out{<i>Med</i>(<i>&pi;</i><sub>0</sub>)}}{\eqn{Med(\pi_0)}}.
#' \item \code{"alt_med"}: The design which minimises
#' \ifelse{html}{\out{<i>Med</i>(<i>&pi;</i><sub>1</sub>)}}{\eqn{Med(\pi_1)}}.
#' \item \code{"prior"}: Either the design which minimises
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i>)}}{\eqn{ESS(\pi)}} for the value of
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} specified using
#' \code{"point_prior"}. Or, the design which minimises
#' \ifelse{html}{\out{&#x222b;<i>ESS</i>(<i>&pi;</i>)<i>Beta</i>(<i>&pi;</i>,
#' <i>a</i>,<i>b</i>)d<i>&pi;</i>}}{\eqn{\int_0^1ESS(\pi)Beta(\pi,a,b)d\pi}}
#' over [0,1], where \ifelse{html}{\out{<i>Beta</i>(<i>&pi;</i>,<i>x</i>,
#' <i>y</i>)}}{\eqn{Beta(\pi,x,y)}} is the PDF of a beta distribution with shape
#' parameters \ifelse{html}{\out{<i>x</i>}}{\eqn{x}} and \ifelse{html}{\out{<i>
#' y</i>}}{\eqn{y}}, specified through \code{"beta_prior"} as a vector,
#' evaluated at point \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}.
#' }
#'
#' Note that when \ifelse{html}{\out{<i>J</i> &le; 3}}{\eqn{J \le 3}}, the
#' optimal design is determined by an exhaustive search. This means that vast
#' speed improvements can be made by carefully choosing the values of
#' \code{Nmin} and \code{Nmax}. In contrast, if \ifelse{html}{\out{<i>J</i> &gt;
#' 3}}{\eqn{J > 3}}, simulated annealing is employed to stochastically search
#' for the optimal design, as proposed by Chen and Lee (2013).
#'
#' @param J The maximal number of stages to allow.
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desirable) response probability at which the trial is
#' powered.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param Nmin The minimal total sample size to allow in considered
#' designs.
#' @param Nmax The maximal total sample size to allow in considered designs.
#' @param futility A logical variable indicating whether early stopping for
#' futility should be allowed.
#' @param efficacy A logical variable indicating whether early stopping for
#' efficacy should be allowed.
#' @param optimality Choice of optimal design criteria. Must be one of
#' \code{"null_ess"}, \code{"alt_ess"}, \code{"null_med"}, \code{"alt_med"},
#' \code{"minimax"} or \code{"prior"}.
#' @param point_prior Value of the response probability to minimise the expected
#' sample size at. Only (potentially) required if \code{optimality == "prior"}.
#' @param beta_prior Shape parameters of the beta distribution to optimise the
#' expected sample over. Only (potentially) required if
#' \code{optimality == "prior"}.
#' @param equal_n A logical variable indicating that the sample size of each
#' stage should be equal.
#' @param ensign A logical variable indicating that the design of Ensign
#' \emph{et al.} (1994) should be mimicked, and the first stage futility
#' boundary forced to be 0.
#' @param summary A logical variable indicating a summary of the function's
#' progress should be printed to the console.
#' @return A list of class \code{"sa_des_gs"} containing the following elements
#' \itemize{
#' \item A list in the slot \code{$des} containing details of the identified
#' optimal design.
#' \item A tibble in the slot \code{$feasible}, consisting of the
#' identified designs which met the required operating characteristics.
#' \item Each of the input variables as specified.
#' }
#' @examples
#' # The minimax design for the default parameters
#' minimax   <- des_gs(optimality = "minimax")
#' # The corresponding design minimising the expected
#' # sample size under the null hypothesis
#' null_ess  <- des_gs(optimality = "null_ess")
#' # The corresponding 3-stage minimax design
#' minimax_3 <- des_gs(J = 3)
#' @references Chen K, Shan M (2008) Optimal and minimax three-stage designs for
#' phase II oncology trials. \emph{Contemporary Clinical Trials}
#' \strong{29:}32-41.
#' @references Chen N, Lee JJ (2013) Optimal continuous-monitoring design of
#' single-arm phase II trial based on the simulated annealing method
#' \emph{Contemporary Clinical Trials} \strong{35:}170-8.
#' @references Chen TT (1997) Optimal three-stage designs for phase II cancer
#' clinical trials. \emph{Statistics in Medicine} \strong{16:}2701-11.
#' @references Ensign LG \emph{et al.} (1994) An optimal three-stage design for
#' phase II clinical trials. \emph{Statistics in Medicine} \strong{13:}1727-36.
#' @references Hanfelt JJ \emph{et al.} (1999) A modification of Simon's optimal
#' design for phase II trials when the criterion is median sample size.
#' \emph{Controlled Clinical Trials} \strong{20:}555-66.
#' @references Mander AP, Thompson SG (2010) Two-stage designs optimal under the
#' alternative hypothesis for phase II cancer clinical trials.
#' \emph{Contemporary Clinical Trials} \strong{31:}572-8.
#' @references Simon R (1989) Optimal two-stage designs for phase II clinical
#' trials. \emph{Controlled Clinical Trials} \strong{10:}1-10.
#' @references Shuster J (2002) Optimal two-stage designs for single arm phase
#' II cancer trials. \emph{Journal of Biopharmaceutical Statistics}
#' \strong{12:} 39-51.
#' @seealso \code{\link{opchar_gs}}, \code{\link{est_gs}},
#' \code{\link{pval_gs}}, \code{\link{ci_gs}}, and their associated \code{plot}
#' family of functions. Note that similar functionality is available through
#' \code{\link[clinfun]{ph2simon}} and
#' \code{\link[OneArmPhaseTwoStudy]{getSolutions}}.
#' @export
des_gs <- function(J = 2, pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                   Nmin = 1, Nmax = 30, futility = T, efficacy = F,
                   optimality = "null_ess", point_prior, beta_prior,
                   equal_n = F, ensign = F, summary = F){

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(1, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), 1)
  check_real_range_strict(beta, "beta", c(0, 1), 1)
  check_integer_pair_range(Nmin, Nmax, "Nmin", "Nmax", c(0, Inf))
  check_stopping(futility, efficacy)
  if (missing(point_prior)) {
    point_prior <- NULL
  }
  if (missing(beta_prior)) {
    beta_prior  <- NULL
  }
  check_optimality(optimality, point_prior, beta_prior, futility, efficacy)
  check_logical(equal_n, "equal_n")
  check_ensign(ensign, futility)
  check_logical(summary, "summary")
  if (J > 3) {
    stop("J > 3 is not currently supported.")
  }

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Design of a ", J, "-Stage Group Sequential Single-arm Trial with a Single Binary Endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("\nYou have chosen to test the following hypotheses\n")
    message("     H\u2080: \u03c0 \u2264 \u03c0\u2080 = ", pi0, ", H\u2081: \u03c0 > \u03c0\u2080 = ", pi0, ".\n")
    message("with the following error constraints\n")
    message("     P(\u03c0\u2080) = P(", pi0, ") \u2264 \u03b1 = ", alpha, ", P(\u03c0\u2081) = P(", pi1, ") \u2265 1 - \u03b2 = ", 1 - beta, ".\n")
    Sys.sleep(2)
    message("You have chosen to restrict the allowed possible sample size N = n such that\n")
    message("  \u2022 N \u2265 ", Nmin, ".")
    message("  \u2022 N \u2264 ", Nmax, ".\n")
    if (equal_n) {
      Sys.sleep(2)
      if (J == 2) {
        message("You have chosen to restrict the allowed values of the n\u2c7c, j = 1,2, such that\n")
        message("  \u2022 n\u2081 = n\u2082.\n")
      } else {
        message("You have chosen to restrict the allowed values of the n\u2c7c, j = 1,\u2026,", J, " such that\n")
        message("  \u2022 n\u2081 = \u22ef = n", sub_num(J), ".\n")
      }
    }
    Sys.sleep(2)
    message("You have chosen to restrict the allowed values in a and r such that\n")
    message("  \u2022 a", sub_num(J), " + 1 = r", sub_num(J), ".")
    if (!futility) {
      if (J == 2) {
        message("  \u2022 a\u2081 = -\u221e.\n")
      } else if (J == 3) {
        message("  \u2022 a\u2081 = a\u2082 = -\u221e.\n")
      } else {
        message("  \u2022 a\u2081 = \u22ef = a", sub_num(J - 1), "= -\u221e.")
      }
    }
    if (ensign) {
      message("  \u2022 a\u2081 = 0.")
    }
    if (!efficacy) {
      if (J == 2) {
        message("  \u2022 r\u2081 = \u221e.\n")
      } else if (J == 3) {
        message("  \u2022 r\u2081 = r\u2082 = \u221e.\n")
      } else {
        message("  \u2022 r\u2081 = \u22ef = r", sub_num(J - 1), "= \u221e.")
      }
    }
    Sys.sleep(2)
    if (J <= 3){
      message("\nNote: For J = ", J, ", the optimal design will be determined using an exhaustive search.\n")
    } else {
      message("\nNote: For J = ", J, ", the optimal design will be determined using simulated annealing.\n")
    }
    Sys.sleep(2)
    message("\nBeginning the required calculations...")
  }

  ##### Main Computations ######################################################

  if (J <= 3) {
    feasible <- saGS(J, pi0, pi1, alpha, beta, Nmin, Nmax, as.numeric(futility),
                     as.numeric(efficacy), as.numeric(equal_n),
                     as.numeric(ensign), as.numeric(summary))
    if (feasible[1, 1] > 0) {
      if (summary) {
        message("...feasible designs in range of considered maximal allowed sample size identified....")
        Sys.sleep(2)
        message("...now identifying optimal design(s)...")
      }
      if (nrow(feasible) == 1){
        feasible         <- matrix(c(feasible[, 1:(2*J + 1)],
                                     feasible[, 2*J] + 1,
                                     feasible[, (2*J + 2):ncol(feasible)]),
                                   nrow = 1, ncol = 16 + 5*(J == 3))
      } else {
        feasible         <- cbind(feasible[, 1:(2*J + J - 1)], feasible[, 2*J] + 1,
                                  feasible[, (2*J + J):ncol(feasible)])
      }
      feasible           <- tibble::as.tibble(feasible)
      colnames(feasible) <- c(paste(rep(c("n", "a", "r"), each = J),
                                    rep(1:J, 3), sep = ""), "P(pi0)", "P(pi1)",
                              paste(rep("PET", 2*(J - 1)),
                                    rep(1:(J - 1), each = 2),
                                    rep("(pi", 2*(J - 1)), rep(0:1, J - 1), ")",
                                    sep = ""), "ESS(pi0)", "ESS(pi1)",
                              "Med(pi0)", "Med(pi1)", "VSS(pi0)", "VSS(pi1)")
      feasible            <- dplyr::mutate(feasible,
                                           `max(N)` = rowSums(feasible[, 1:J]))
      feasible[, 1:(3*J)] <- dplyr::mutate_if(feasible[, 1:(3*J)], is.double,
                                              as.integer)
      if (!futility) {
        feasible[, (J + 1):(2*J - 1)]   <- -Inf
      }
      if (!efficacy) {
        feasible[, (2*J + 1):(3*J - 1)] <- Inf
      }
      if (optimality == "null_ess") {
        feasible <- dplyr::arrange(feasible, `ESS(pi0)`, `max(N)`)
      } else if (optimality == "alt_ess") {
        feasible <- dplyr::arrange(feasible, `ESS(pi1)`, `max(N)`)
      } else if (optimality == "null_med") {
        feasible <- dplyr::arrange(feasible, `Med(pi0)`, `max(N)`)
      } else if (optimality == "alt_med") {
        feasible <- dplyr::arrange(feasible, `Med(pi1)`, `max(N)`)
      } else if (optimality == "minimax") {
        feasible <- dplyr::arrange(feasible, `max(N)`, `ESS(pi0)`)
      } else if (optimality == "prior") {
        if (!is.null(point_prior)) {
          ESS_point      <- numeric(nrow(feasible))
          for (i in 1:nrow(feasible)) {
            ESS_point[i] <- int_opchar_gs(point_prior, J,
                                          as.numeric(feasible[i, (J + 1):(2*J)]),
                                          as.numeric(feasible[i, (2*J + 1):(3*J)]),
                                          as.numeric(feasible[i, 1:J]),
                                          cumsum(as.numeric(feasible[i, 1:J])),
                                          1:J)$ESS
            if (all(summary, i%%1000 == 0)) {
              message("...", i, " feasible designs evaluated for optimality...")
            }
          }
          feasible       <- dplyr::mutate(feasible, `ESS(point)` = ESS_point)
          feasible       <- dplyr::arrange(feasible, `ESS(point)`, `ESS(pi0)`)
        } else {
          ESS_beta      <- numeric(nrow(feasible))
          for (i in 1:nrow(feasible)) {
            ESS_beta[i] <- stats::integrate(ess_gs_beta, lower = 0, upper = 1, J = J,
                                     a = as.numeric(feasible[i, (J + 1):(2*J)]),
                                     r = as.numeric(feasible[i,
                                                             (2*J + 1):(3*J)]),
                                     n = as.numeric(feasible[i, 1:J]),
                                     N = cumsum(as.numeric(feasible[i, 1:J])),
                                     beta_prior = beta_prior)$value
            if (all(summary, i%%100 == 0)) {
              message("...", i, " feasible designs evaluated for optimality...")
            }
          }
          feasible      <- dplyr::mutate(feasible, `ESS(beta)` = ESS_beta)
          feasible      <- dplyr::arrange(feasible, `ESS(beta)`, `ESS(pi0)`)
        }
      }
      des <- list(J = J, n = as.numeric(feasible[1, 1:J]),
                  a = as.numeric(feasible[1, (J + 1):(2*J)]),
                  r = as.numeric(feasible[1, (2*J + 1):(3*J)]),
                  pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                  opchar = feasible[1, ])
    } else {
      feasible <- des <- NULL
      if (summary) {
        message("...no feasible designs found in range of considered maximal allowed sample size. Consider decreasing Nmin and increasing Nmax...")
      }
    }
  } else {
    # Required simulated annealing code
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = feasible, J = J, pi0 = pi0,
                        pi1 = pi1, alpha = alpha, beta = beta, Nmin = Nmin,
                        Nmax = Nmax, futility = futility, efficacy = efficacy,
                        optimality = optimality, point_prior = point_prior,
                        beta_prior = beta_prior, equal_n = equal_n,
                        ensign = ensign, summary = summary)
  class(output) <- "sa_des_gs"
  return(output)

}
