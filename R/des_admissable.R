#' Design an admissable group sequential single-arm trial for a single binary
#' endpoint
#'
#' Determines admissable group sequential single-arm clinical trial designs for
#' a single binary primary endpoint, using exact binomial calculations.
#'
#' \code{des_admissable()} supports the determination of admissable
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
#' (<i>a</i><sub>1</sub>,&hellip;,<i>a</i><sub><i>J</i></sub>)}}{\eqn{\bold{a}=
#' (a_1,\dots,a_J)}}, \ifelse{html}{\out{<b><i>r</i></b> = (<i>r</i><sub>1
#' </sub>,&hellip;,<i>r</i><sub><i>J</i></sub>)}}{\eqn{\bold{r}=
#' (r_1,\dots,r_J)}}, and \ifelse{html}{\out{<b><i>n</i></b> = (<i>n</i><sub>1
#' </sub>,&hellip;,<i>n</i><sub><i>J</i></sub>)}}{\eqn{\bold{n}=
#' (n_1,\dots,n_J)}}.
#'
#' With these vectors, and denoting the number of responses after
#' \ifelse{html}{\out{<i>m</i>}}{\eqn{m}} patients have been observed by
#' \ifelse{html}{\out{<i>s</i><sub><i>m</i></sub>}}{\eqn{s_m}}, the stopping
#' rules for a trial are then as follows
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
#' values for \ifelse{html}{\out{<i>N</i><sub><i>J</i></sub>}}{\eqn{N_J}}, while
#' if set to \code{TRUE}, \code{equal_n} enforces that \ifelse{html}{\out{<i>n
#' </i><sub>1</sub> = &ctdot; = <i>n</i><sub><i>J</i></sub>}}{\eqn{n_1 = \dots
#' = n_J}}.
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
#' 0}}, as suggested for 3-stage designs by Ensign \emph{et al} (1994).
#'
#' Note that to ensure a decision is made about \ifelse{html}{\out{<i>H</i><sub>
#' 0</sub>}}{\eqn{H_0}}, this function enforces that
#' \ifelse{html}{\out{<i>a</i><sub><i>J</i></sub> + 1 = <i>r</i><sub>
#' <i>J</i>}}{\eqn{a_J + 1 = r_J}}.
#'
#' To describe the admissability criteria, denote the expected sample size
#' when the true response probability is
#' \ifelse{html}{\out{<i>&pi;</i>}}{\eqn{\pi}} by
#' \ifelse{html}{\out{<i>ESS</i>(<i>&pi;</i>)}}{\eqn{ESS(\pi)}}.
#' Then, the following admissability criteria are currently supported:
#'
#' The designs which minimise
#' \ifelse{html}{\out{<i>w</i><sub>0</sub><i>ESS</i>(<i>&pi;</i><sub>0</sub>) +
#' <i>w</i><sub>1</sub><i>ESS</i>(<i>&pi;</i><sub>1</sub>) +
#' (1 - <i>w</i><sub>0</sub> - <i>w</i><sub>1</sub>)<i>N</i><sub><i>J</i></sub>
#' }}{\eqn{w_0ESS(\pi_0)+w_1ESS(\pi_1)+(1-w_0-w_1)N_J}} for
#' \ifelse{html}{\out{0 &le; <i>w</i><sub>0</sub> + <i>w</i><sub>1</sub> &le;
#' 1}}{\eqn{0 \le w_0+w_1 \le 1}}.
#'
#' @param J The maximal number of stages to allow.
#' @param pi0 The (undesirable) response probability used in the definition of
#' the null hypothesis.
#' @param pi1 The (desiable) response probability at which the trial is powered.
#' @param alpha The desired maximal type-I error-rate.
#' @param beta The desired maximal type-II error-rate.
#' @param Nmin The minimal total sample size to allow in considered
#' designs.
#' @param Nmax The maximal total sample size to allow in considered designs.
#' @param futility A logical variable indicating whether early stopping for
#' futility should be allowed.
#' @param efficacy A logical variable indicating whether early stopping for
#' efficacy should be allowed.
#' @param equal_n A logical variable indicating that the sample size of each
#' stage should be equal.
#' @param ensign A logical variable indicating that the design of Ensign
#' \emph{et al} (1994) should be mimicked, and the first stage futility
#' boundary forced to be 0.
#' @param summary A logical variable indicating a summary of the function's
#' progress should be printed to the console.
#' @return A list of class \code{"sa_des_admissable"} containing the following
#' elements
#' \itemize{
#' \item A list of lists in the slot \code{$des} containing details of the
#' identified admissable design(s). Each element will contain
#' details of each admissable design.
#' \item A tibble in the slot \code{$feasible}, consisting of the
#' identified designs which met the required operating characteristics.
#' \item A tibble in the slot \code{$admissable},
#' summarising the performance of the admissable designs.
#' \item A tibble in the slot \code{$weights} containing details of which design
#' is admissable for each considered combination of \ifelse{html}{\out{<i>w</i>
#' <sub>0</sub>}}{\eqn{w_0}} and \ifelse{html}{\out{<i>w</i><sub>1
#' </sub>}}{\eqn{w_1}}.
#' \item Each of the input variables as specified.
#' }
#' @examples
#' # The admissable designs for the default parameters
#' admissable <- des_admissable()
#' @seealso \code{\link{opchar_admissable}}, and their associated \code{plot}
#' family of functions.
#' @export
des_admissable <- function(J = 2, pi0 = 0.1, pi1 = 0.3, alpha = 0.05,
                           beta = 0.2, Nmin = 1, Nmax = 50, futility = T,
                           efficacy = F, equal_n = F, ensign = F, summary = F) {

  ##### Input Checking #########################################################

  check_integer_range(J, "J", c(1, Inf))
  check_real_pair_range_strict(pi0, pi1, "pi0", "pi1", c(0, 1))
  check_real_range_strict(alpha, "alpha", c(0, 1), "1")
  check_real_range_strict(beta, "beta", c(0, 1), "1")
  check_integer_pair_range(Nmin, Nmax, "Nmin", "Nmax", c(0, Inf))
  check_stopping(futility, efficacy)
  check_logical(equal_n, "equal_n")
  check_ensign(ensign, futility)
  check_logical(summary, "summary")
  if (J > 3) {
    stop("J > 3 is not currently supported.")
  }

  ##### Print Summary ##########################################################

  if (summary){
    message(rep("-", 10))
    message("Design of ", J, "-stage admissable group sequential single-arm trials with a single binary endpoint")
    message(rep("-", 10))
    Sys.sleep(2)
    message("\nYou have chosen to test the following hypotheses\n")
    message("     H₀: π ≤ π₀ = ", pi0, ", H₁: π > π₀ = ", pi0, ".\n")
    message("with the following error constraints\n")
    message("     P(π₀) = P(", pi0, ") ≤ α = ", alpha, ", P(π₁) = P(", pi1, ") ≥ 1 - β = ", 1 - beta, ".\n")
    Sys.sleep(2)
    message("You have chosen to restrict the allowed maximal possible sample size N such that\n")
    message("  • N ≥ ", Nmin, ".")
    message("  • N ≤ ", Nmax, ".\n")
    if (equal_n) {
      Sys.sleep(2)
      if (J == 2) {
        message("You have chosen to restrict the allowed values of the nⱼ, j = 1,2, such that\n")
        message("  • n₁ = n₂.\n")
      } else {
        message("You have chosen to restrict the allowed values of the nⱼ, j = 1,…,", J, " such that\n")
        message("  • n₁ = ⋯ = n", sub_num(J), ".\n")
      }
    }
    Sys.sleep(2)
    message("You have chosen to restrict the allowed values in a and r such that\n")
    message("  • a", sub_num(J), " + 1 = r", sub_num(J), ".")
    if (!futility) {
      if (J == 2) {
        message("  • a₁ = -∞.\n")
      } else if (J == 3) {
        message("  • a₁ = a₂ = -∞.\n")
      } else {
        message("  • a₁ = ⋯ = a", sub_num(J - 1), "= -∞.")
      }
    }
    if (ensign) {
      message("  • a", sub_num(1), " = 0.")
    }
    if (!efficacy) {
      if (J == 2) {
        message("  • r₁ = ∞.\n")
      } else if (J == 3) {
        message("  • r₁ = r₂ = ∞.\n")
      } else {
        message("  • r₁ = ⋯ = r", sub_num(J - 1), "= ∞.")
      }
    }
    Sys.sleep(2)
    if (J <= 3){
      message("\nNote: For J = ", J, ", the optimal design will be determined using an exhaustive search.\n")
    } else {
      message("\nNote: For J = ", J, ", the optimal design will be determined using simulated annealing.\n")
    }
    Sys.sleep(2)
    message("\nNow beginning the required calculations...")
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
        message("...now identifying admissable design(s)...")
      }
      if (nrow(feasible) == 1){
        feasible         <- matrix(c(feasible[, 1:(2*J + 1)],
                                     feasible[, 2*J] + 1,
                                     feasible[, (2*J + 2):ncol(feasible)]),
                                   nrow = 1, ncol = 16 + 5*(J == 3))
      } else {
        feasible         <- cbind(feasible[, 1:(2*J + 1)], feasible[, 2*J] + 1,
                                  feasible[, (2*J + 2):ncol(feasible)])
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
      feasible          <- dplyr::arrange(feasible, `ESS(pi0)`, `max(N)`)
      weights           <- expand.grid(w0 = seq(from = 0, to = 1, by = 0.005),
                                       w1 = seq(from = 0, to = 1,
                                                by = 0.005))[-1, ]
      weights           <- weights[which(weights[, 1] + weights[, 2] <= 1), ]
      opt_indices       <- numeric(nrow(weights))
      for (i in 1:nrow(weights)){
        opt_indices[i]  <- which.min(weights[i, 1]*feasible$`ESS(pi0)` +
                                       weights[i, 2]*feasible$`ESS(pi1)` +
                                       (1 - weights[i, 1] - weights[i, 2])*
                                       feasible$`max(N)`)
      }
      uniq_opt_indices  <- unique(opt_indices)
      admissable        <- feasible[uniq_opt_indices, ]
      des               <- list()
      optimal_des       <- numeric(nrow(weights))
      for (i in 1:length(uniq_opt_indices)) {
        des[[i]]        <- list(J = J,
                                n = as.numeric(feasible[uniq_opt_indices[i],
                                                        1:J]),
                                a = as.numeric(feasible[uniq_opt_indices[i],
                                                        (J + 1):(2*J)]),
                                r = as.numeric(feasible[uniq_opt_indices[i],
                                                        (2*J + 1):(3*J)]),
                                pi0 = pi0, pi1 = pi1, alpha = alpha,
                                beta = beta,
                                opchar = feasible[uniq_opt_indices[i], ])
        class(des[[i]]) <- "sa_des_gs"
        optimal_des[which(opt_indices == uniq_opt_indices[i])] <-
          paste("Design", i)
      }
      weights <- tibble::tibble(w0 = weights[, 1], w1 = weights[, 2],
                                Design = factor(optimal_des,
                                                levels = unique(optimal_des)))
      new_levels      <- levels(weights$Design)
      if (efficacy == F) {
        for (i in 1:length(new_levels)) {
          new_levels[i] <- paste(paste(new_levels[i], ": ",
                                       paste(des[[i]]$a[1:(J - 1)], "/",
                                             cumsum(des[[i]]$n[1:(J - 1)]),
                                             sep = "",
                                             collapse = ", "), sep = "",
                                       collapse = ""), ", ",
                                 des[[i]]$a[J], "/",
                                 sum(des[[i]]$n), sep = "", collapse = "")
        }
      } else {
        for (i in 1:length(new_levels)) {
          new_levels[i] <- paste(paste(new_levels[i], ": ",
                                       paste("(", des[[i]]$a[1:(J - 1)], ",",
                                             des[[i]]$r[1:(J - 1)], ")/",
                                             cumsum(des[[i]]$n[1:(J - 1)]),
                                             sep = "",
                                             collapse = ", "), sep = "",
                                       collapse = ""), ", ",
                                 des[[i]]$a[J], "/",
                                 sum(des[[i]]$n), sep = "", collapse = "")
        }
      }
      weights$Design   <- plyr::mapvalues(weights$Design,
                                          from = levels(weights$Design),
                                          to = new_levels)
      admissable <- tibble::as_tibble(cbind(new_levels, admissable))
      colnames(admissable)[1] <- "Design"
      admissable$Design <- as.factor(admissable$Design)
    } else {
      feasible <- des <- admissable <- weights <- NULL
      if (summary) {
        message("...no feasible designs found in range of considered maximal allowed sample size. Consider decreasing Nmin and increasing Nmax...")
      }
    }
  }

  ##### Outputting #############################################################

  if (summary) {
    message("...outputting.")
  }
  output          <- list(des = des, feasible = feasible,
                          admissable = admissable, weights = weights, J = J,
                          pi0 = pi0, pi1 = pi1, alpha = alpha, beta = beta,
                          Nmin = Nmin, Nmax = Nmax, futility = futility,
                          efficacy = efficacy, equal_n = equal_n,
                          ensign = ensign, summary = summary)
  class(output) <- "sa_des_admissable"
  return(output)
}
