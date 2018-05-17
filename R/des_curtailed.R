#' Design a curtailed group sequential single-arm trial for a single binary
#' endpoint
#'
#' Determines curtailed group sequential single-arm clinical trial designs for a
#' single binary primary endpoint.
#'
#' \code{des_curtailed()} supports the determination of a variety of (optimised)
#' curtailed group sequential single-arm clinical trial designs for a single
#' binary primary endpoint. For all supported designs, the following hypotheses
#' are tested for the response probability
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}}
#'
#' \ifelse{html}{\out{<center><i>H</i><sub>0</sub> : <i>&pi;</i> &le; <i>&pi;
#' </i><sub>0</sub>, <i>H</i><sub>1</sub> : <i>&pi;</i> > <i>&pi;</i><sub>
#' 0</sub>,</center>}}{\deqn{H_0 : \pi \le \pi_0,\qquad H_1 : \pi > \pi_0,}}
#'
#' for \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\eqn{\pi_0}}, specified
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
#' these vectors, the chosen optimality criteria, and the chosen curtailment
#' rule.
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
#' In addition, two vectors \ifelse{html}{\out{<i>&theta;</i><sub>F</sub>
#' )}}{\eqn{\theta_F}} and \ifelse{html}{\out{<i>&theta;</i><sub>E</sub>
#' )}}{\eqn{\theta_E}}, each of length \ifelse{html}{\out{<i>J</i>}}{\eqn{J}},
#' must be specified which determine the chosen curtailment rule.
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
#' Note that the optimal design is determined by an exhaustive search. This
#' means that vast speed improvements can be made by carefully choosing the
#' values of \code{Nmin} and \code{Nmax}.
#'
#' @param J The maximal number of stages to allow.
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desirable) response probability at which the trial is
#' powered.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param thetaF The vector of futility curtailment probabilities.
#' @param thetaE The vector of efficacy curtailment probabilities.
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
#' # The minimax design for the default parameters with NSC
#' minimax <- des_curtailed(optimality = "minimax")
#' @seealso \code{\link{opchar_curtailed}}, and their associated \code{plot}
#' family of functions.
#' @export
des_curtailed <- function(J = 2, pi0 = 0.1, pi1 = 0.3, alpha = 0.05, beta = 0.2,
                          thetaF = rep(0, J), thetaE = rep(1, J), Nmin = 1,
                          Nmax = 50, futility = T, efficacy = F,
                          optimality = "null_ess", point_prior, beta_prior,
                          equal_n = F, ensign = F, summary = F) {

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(0, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), "1")
  check_real_range_strict(beta, "beta", c(0, 1), "1")
  check_real_range_strict(thetaF, "thetaF", c(-Inf, Inf), J)
  check_real_range_strict(thetaE, "thetaE", c(-Inf, Inf), J)
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
    message("Design of ", J, "-stage curtailed group-sequential single-arm trials with a single binary endpoint")
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

  if (J == 1) {
    des <- des_fixed(pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                     Nmin = Nmin, Nmax = Nmax, exact = T, summary = summary)
  } else {
    des <- des_gs(J = J, pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                  Nmin = Nmin, Nmax = Nmax, futility = futility,
                  efficacy = efficacy, optimality = optimality,
                  point_prior = point_prior, beta_prior = beta_prior,
                  equal_n = equal_n, ensign = ensign, summary = summary)
  }
  if (!is.null(des$des)) {
    a      <- des$des$a
    r      <- des$des$r
    n      <- des$des$n
    cp_mat <- matrix(0, nrow = sum(n) + 1, ncol = sum(n))
    for (m in 1:sum(n)) {
      for (s in 0:m) {
        cp_mat[s + 1, m] <- cr_gs(pi1, s, m, J, a, r, n)
      }
    }
    thetaEog <- thetaE
    thetaFog <- thetaF
    thetaE[which(thetaE == 1)] <- 1 - 1e-10
    thetaF[which(thetaF == 0)] <- 1e-10
    a_curt <- rep(-Inf, sum(n))
    r_curt <- rep(Inf, sum(n))
    n_curt <- rep(1, sum(n))
    for (i in 1:sum(n)) {
      j <- min(which(i <= cumsum(n)))
      if (any(cp_mat[1:(i + 1), i] <= thetaF[j])) {
        a_curt[i] <- max(which(cp_mat[1:(i + 1), i] <= thetaF[j])) - 1
      }
      if (any(cp_mat[, i] >= thetaE[j])) {
        r_curt[i] <- min(which(cp_mat[, i] >= thetaE[j])) - 1
      }
    }
    opchar <- int_opchar_curtailed(c(pi0, pi1), J, sum(n), a_curt, r_curt, n,
                                   n_curt, cumsum(n), cumsum(n_curt), 1:J)
    des <- list(J = J, a = des$des$a, r = des$des$r, n = des$des$n, pi0 = pi0,
                pi1 = pi1, alpha = alpha, beta = beta, J_curt = sum(n),
                a_curt = a_curt, r_curt = r_curt, n_curt = n_curt,
                thetaF = thetaFog, thetaE = thetaEog, opchar = opchar)
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output        <- list(des = des, feasible = des$feasible, J = J, pi0 = pi0,
                        pi1 = pi1, alpha = alpha, beta = beta,
                        thetaF = thetaFog, thetaE = thetaEog, Nmin = Nmin,
                        Nmax = Nmax, futility = futility, efficacy = efficacy,
                        optimality = optimality, point_prior = point_prior,
                        beta_prior = beta_prior, equal_n = equal_n,
                        ensign = ensign, summary = summary, thetaF = thetaF,
                        thetaE = thetaE)
  class(output) <- "sa_des_curtailed"
  return(output)
}
