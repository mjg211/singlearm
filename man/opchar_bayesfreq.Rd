% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opchar_bayesfreq.R
\name{opchar_bayesfreq}
\alias{opchar_bayesfreq}
\title{Determine the operating characteristics of Bayesian-frequentist single-arm
trial designs for a single binary endpoint}
\usage{
opchar_bayesfreq(des, ..., k, mu, nu, pi, summary = F)
}
\arguments{
\item{des}{An object of class \code{"sa_des_bayesfreq"}, as returned by
\code{des_bayesfreq()}.}

\item{...}{Additional objects of class \code{"sa_des_bayesfreq"}. These will
be grouped in to a list named \code{"add_des"}.}

\item{k}{Calculations are performed conditional on the trial stopping in one
of the stages listed in vector \code{k}. Thus, \code{k} should be a vector of
integers, with elements between one and the maximal number of possible stages
in the supplied designs. If left unspecified, it will internally default to
all possible stages.}

\item{mu}{A vector of the first Beta shape parameters to evaluate operating
characteristics at. This will internally default to be the
\ifelse{html}{\out{<i>&mu;</i>}}{\deqn{\mu}} from the supplied designs if it
left unspecified.}

\item{nu}{A vector of the second Beta shape parameters to evaluate operating
characteristics at. This will internally default to be the
\ifelse{html}{\out{<i>&nu;</i>}}{\deqn{\mu}} from the supplied designs if it
left unspecified.}

\item{pi}{A vector of response probabilities to evaluate operating
characteristics at. This will internally default to be the
\ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}} and
\ifelse{html}{\out{<i>&pi;</i><sub>1</sub>}}{\deqn{\pi_1}} from the
supplied designs if it is left unspecified.}

\item{summary}{A logical variable indicating whether a summary of the
function's progress should be printed to the console.}
}
\value{
A list of class \code{"sa_opchar_bayesfreq"} containing the following
elements
\itemize{
\item A tibble in the slot \code{$opchar_bayes} summarising the Bayesian
operating characteristics of the supplied designs.
\item A tibble in the slot \code{$opchar_freq} summarising the frequentist
operating characteristics of the supplied designs.
\item Each of the input variables as specified, subject to internal
modification.
}
}
\description{
\code{opchar_bayesfreq()} supports the simultaneous evaluation of the
operating characteristics of multiple Bayesian-frequentist single-arm
clinical trial designs for a single binary primary endpoint, determined using
\code{des_bayesfreq()}.
}
\details{
Note that each of the supplied designs must have been designed for the same
values of \ifelse{html}{\out{<i>&mu;</i>}}{\deqn{\mu}},
\ifelse{html}{\out{<i>&nu;</i>}}{\deqn{\nu}}, and
\ifelse{html}{\out{<i>&pi;</i><sub>0</sub>}}{\deqn{\pi_0}}.

For each value of \ifelse{html}{\out{<i>mu</i>}}{\eqn{\mu}},
\ifelse{html}{\out{<i>nu</i>}}{\eqn{\nu}}, and
\ifelse{html}{\out{<i>pi</i>}}{\eqn{\pi}} in
the supplied vectors \ifelse{html}{\out{<b><i>mu</i></b>}}{\eqn{\bold{\mu}}},
\ifelse{html}{\out{<b><i>nu</i></b>}}{\eqn{\bold{\nu}}}, and
\ifelse{html}{\out{<b><i>pi</i></b>}}{\eqn{\bold{\pi}}},
\code{opchar_bayesfreq()} evaluates the Bayesian and frequentist power, ESS,
and other key characteristics, of each of the supplied designs.

Calculations are performed conditional on the trial stopping in one of the
stages specified using the input (vector) \code{k}.
}
\examples{
# Find the optimal two-stage Bayesian-frequentist design for the default
# parameters
des    <- des_bayesfreq()
# Determine operating characteristics for a range of mu, nu, and pi
opchar <- opchar_bayesfreq(des, mu = seq(0.05, 0.2, length.out = 10),
                           nu = seq(0.45, 1.8, length.out = 100),
                           pi = seq(0, 1, by = 0.01))
}
\seealso{
\code{\link{des_bayesfreq}}, and their associated \code{plot} family
of functions.
}
