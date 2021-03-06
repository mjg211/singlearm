% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opchar_gehan.R
\name{opchar_gehan}
\alias{opchar_gehan}
\title{Determine the operating characteristics of Gehan two-stage single-arm
trial designs for a single binary endpoint}
\usage{
opchar_gehan(des, ..., k, pi, summary = F)
}
\arguments{
\item{des}{An object of class \code{"sa_des_gehan"}, as returned by
\code{des_gehan()}.}

\item{...}{Additional objects of class \code{"sa_des_gehan"}. These will
be grouped in to a list named \code{"add_des"}.}

\item{k}{Calculations are performed conditional on the trial stopping in one
of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
integers, with elements between one and two. If left unspecified, it will
internally default to all possible stages.}

\item{pi}{A vector of response probabilities to evaluate operating
characteristics at. This will internally default to be the
\ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
\ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
supplied designs if it is left unspecified.}

\item{summary}{A logical variable indicating whether a summary of the
function's progress should be printed to the console.}
}
\value{
A list of class \code{"sa_opchar_gehan"} containing the following
elements
\itemize{
\item A tibble in the slot \code{$opchar} summarising the operating
characteristics of the supplied designs.
\item Each of the input variables as specified, subject to internal
modification.
}
}
\description{
\code{opchar_gehan()} supports the simultaneous evaluation of the
operating characteristics of multiple Gehan two-stage single-arm clinical
trial designs for a single binary primary endpoint, determined using
\code{des_gehan()}.
}
\details{
Note that each of the supplied designs must have been designed for the same
values of \ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
\code{find_D}.

For each value of \ifelse{html}{\out{<i>pi</i>)}}{\eqn{\pi}} in
the supplied vector \ifelse{html}{\out{<b><i>pi</i></b>)}}{\eqn{\bold{\pi}}},
\code{opchar_gehan()} evaluates the power (if find_D = \code{T}), ESS, and
other key characteristics, of each of the supplied designs.

Calculations are performed conditional on the trial stopping in one of the
stages specified using the input (vector) \code{k}.
}
\examples{
# Find the optimal adaptive design for the default parameters
des    <- des_gehan()
# Find its operating characteristics for a range of possible response
# probabilities
opchar <- opchar_gehan(des)
}
\seealso{
\code{\link{des_gehan}}, and their associated \code{plot} family
of functions.
}
