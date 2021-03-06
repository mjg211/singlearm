% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.sa_des_bayesfreq.R
\name{plot.sa_des_bayesfreq}
\alias{plot.sa_des_bayesfreq}
\title{Plot the stopping boundaries of Bayesian-frequentist single-arm trial designs
for a single binary endpoint}
\usage{
\method{plot}{sa_des_bayesfreq}(x, ..., output = F)
}
\arguments{
\item{x}{An object of class \code{"sa_des_bayesfreq"}, as returned by
\code{des_bayesfreq()}.}

\item{...}{Additional objects of class \code{"sa_des_bayesfreq"}. These will
be grouped in to a list named \code{"add_des"}.}

\item{output}{A logical variable indicating whether the outputs described
below should be returned.}
}
\value{
If \code{output = TRUE}, a list containing the following elements is
returned
\itemize{
\item A list in the slot \code{$plot_des} containing the available plots.
\item Each of the input variables as specified, subject to internal
modification.
}
}
\description{
Plots the stopping boundaries of Bayesian-frequentist single-arm trial
designs determined using \code{des_bayesfreq()}. The possible
\ifelse{html}{\out{(<i>s</i>,<i>m</i>)}}{\eqn{(s,m)}} states during the trial
are plotted in a colour coded manner to indicate their associated decision
rules.
}
\details{
Support is available to simultaneously plot the stopping boundaries of
multiple Bayesian-frequentist single-arm clinical trial designs for a single
binary primary endpoint, using faceting.
}
\examples{
# Find the optimal two-stage Bayesian-frequentist design for the default
# parameters
des <- des_bayesfreq()
# Plot the stopping boundaries
plot(des)
}
\seealso{
\code{\link{des_bayesfreq}}, \code{\link{opchar_bayesfreq}}, and
 their associated \code{plot} family of functions.
}
