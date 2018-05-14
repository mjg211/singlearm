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
NULL
