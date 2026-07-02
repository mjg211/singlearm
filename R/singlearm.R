#' singlearm: A package for designing and analysing single-arm clinical trials
#'
#' singlearm provides five categories of important function: \code{des},
#' \code{opchar}, \code{est}, \code{pval}, and \code{ci}. Each of the
#' functions listed below also has its own S3 plotting function.
#'
#' @section \code{des} functions:
#' \itemize{
#' \item \code{\link{des_adaptive}}: Determines a two-stage adaptive single-arm
#' trial design for a single binary primary outcome variable.
#' \item \code{\link{des_admissable}}: Determines two- or three-stage admissable
#' group sequential trial designs for a single binary primary outcome variable.
#' \item \code{\link{des_bayesfreq}}: Determines a single- or two-stage
#' Bayesian-frequentist single-arm trial design for a single binary primary
#' outcome variable.
#' \item \code{\link{des_curtailed}}: Determines curtailed single- or two-stage
#' group sequential single-arm trial designs for a single binary primary outcome
#' variable.
#' \item \code{\link{des_fixed}}: Determines a single-stage single-arm trial
#' design for a single binary primary outcome variable.
#' \item \code{\link{des_gehan}}: Determines two-stage Gehan designs for a
#' single primary outcome variable.
#' \item \code{\link{des_gs}}: Determines two- or three-stage group sequential
#' trial designs for a single binary primary outcome variable.
#' }
#'
#' @section \code{opchar} functions:
#' \itemize{
#' \item \code{\link{opchar_adaptive}}: Determines the operating characteristics
#' of designs determined using \code{\link{des_adaptive}}.
#' \item \code{\link{opchar_admissable}}: Determines the operating
#' characteristics of designs determined using \code{\link{des_admissable}}.
#' \item \code{\link{opchar_bayesfreq}}: Determines the operating
#' characteristics of designs determined using \code{\link{des_bayesfreq}}.
#' \item \code{\link{opchar_curtailed}}: Determines the operating
#' characteristics of designs determined using \code{\link{des_curtailed}}.
#' \item \code{\link{opchar_fixed}}: Determines the operating characteristics of
#' designs determined using \code{\link{des_fixed}}.
#' \item \code{\link{opchar_gs}}: Determines the operating characteristics of
#' designs determined using \code{\link{des_gs}}.
#' }
#'
#' @section \code{est} functions:
#' \itemize{
#' \item \code{\link{est_fixed}}: Determine and evaluate the performance of
#' point estimates for a design determined using \code{\link{des_fixed}}.
#' \item \code{\link{est_gs}}: Determine and evaluate the performance of point
#' estimates for a design determined using \code{\link{des_gs}}.
#' }
#'
#' @section \code{pval} functions:
#' \itemize{
#' \item \code{\link{pval_fixed}}: Determine and evaluate the performance of
#' p-values for a design determined using \code{\link{des_fixed}}.
#' \item \code{\link{pval_gs}}: Determine and evaluate the performance of
#' p-values for a design determined using \code{\link{des_gs}}.
#' }
#'
#' @section \code{ci} functions:
#' \itemize{
#' \item \code{\link{ci_fixed}}: Determine and evaluate the performance of
#' confidence intervals for a design determined using \code{\link{des_fixed}}.
#' \item \code{\link{ci_gs}}: Determine and evaluate the performance of
#' confidence intervals for a design determined using \code{\link{des_gs}}.
#' }
#'
#' @docType package
#' @name singlearm
#' @importFrom stats dnorm integrate pbinom pnorm qnorm
NULL

utils::globalVariables(c(
  # Column names used in dplyr/ggplot2 pipelines
  "Bias(hat(pi)|pi)", "Cover(C|pi)", "Design", "E(L|pi)", "E(hat(pi)|pi)",
  "E(p|pi)", "ESS(beta)", "ESS(mu,nu)", "ESS(pi)", "ESS(pi0)", "ESS(pi1)",
  "ESS(point)", "Med(mu,nu)", "Med(pi)", "Med(pi0)", "Med(pi1)", "O",
  "P(mu,nu)", "P(pi)", "P(pi0)", "P(pi1)", "PB(pi1)",
  "RMSE(hat(pi)|pi)", "VSS(mu,nu)", "VSS(pi)", "Var(L|pi)",
  "Var(hat(pi)|pi)", "Var(p|pi)",
  "a", "bar(L)", "card_P2", "clow(s,m)", "combination", "cupp(s,m)",
  "dbinom_pi0", "dbinom_pi1", "dbinomial_pi0", "dbinomial_pi1",
  "dc", "dc_ef", "dc_pf", "dc_ss", "delta",
  "elem", "en", "end", "f(s,m|pi)", "hat(pi)(s,m)",
  "key", "l(s,m)", "length_dc_ef", "limit", "m", "margPB(pi1)",
  "max(L)", "max(N)", "maxpower", "method", "mu",
  "n1", "n2", "nomalpha", "nombeta", "nu",
  "p(s,m)", "pi(hat)(s,m)",
  "s", "s1", "sample_sizes", "start", "status",
  "store", "umvue", "value", "w", "w0", "w1", "where"
))
