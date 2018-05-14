#' Plot the operating characteristics of curtailed group sequential single-arm trial designs for a
#' single binary endpoint
#'
#' Plots the operating characteristics of curtailed group sequential single-arm trial designs
#' determined using \code{opchar_curtailed()}. A range of plots are available, of which
#' the expected sample size and power curves will be printed by default.
#'
#' @param opchar An object of class \code{"sa_opchar_curtailed"}, as returned by \code{opchar_curtailed()}.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the optimal group sequential design for the default parameters with NSC
#' des         <- des_curtailed()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_curtailed(des, pi = seq(0, 1, 0.01))
#' # Plot the expected sample size and power curve
#' plot(opchar)
#' @seealso \code{\link{des_curtailed}}, \code{\link{opchar_curtailed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_opchar_curtailed <- function(opchar, output = F) {

  ##### Input Checking #########################################################

  check_sa_opchar_curtailed(opchar)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_opchar <- list()

  if (min(opchar$pi) < opchar$des$des$pi0) {
    red   <- tibble::tibble(start = min(opchar$pi),
                            end = min(opchar$des$des$pi0, max(opchar$pi)))
  }
  if (all(min(opchar$pi) <= opchar$des$des$pi1,
          max(opchar$pi) >= opchar$des$des$pi0)) {
    amber <- tibble::tibble(start = max(opchar$des$des$pi0, min(opchar$pi)),
                            end = min(opchar$des$des$pi1, max(opchar$pi)))
  }
  if (max(opchar$pi) > opchar$des$des$pi1) {
    green <- tibble::tibble(start = max(opchar$des$des$pi1,
                                        min(opchar$pi)), end = max(opchar$pi))
  }

  if (is.null(opchar$add_des)) {

    if (length(unique(opchar$pmf$pi)) > 1) {
      plot_opchar$`ESS(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `ESS(pi)`)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(ESS), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))
      print(plot_opchar$`ESS(pi)`)
      plot_opchar$`VSS(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `VSS(pi)`)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(VSS), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))

      plot_opchar$`Med(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `Med(pi)`)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(Med), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))

      plot_opchar$`P(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
        ggplot2::geom_hline(yintercept =
                              c(opchar$des$des$alpha,
                                1 - opchar$des$des$beta),
                            linetype = 2) +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `P(pi)`)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(P), "(", pi,
                                       ")", sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))
      print(plot_opchar$`P(pi)`)
    } else {
      plot_opchar$`ESS(pi)` <- NULL
      plot_opchar$`VSS(pi)` <- NULL
      plot_opchar$`Med(pi)` <- NULL
      plot_opchar$`P(pi)`   <- NULL
    }
  } else {
    if (length(unique(opchar$pmf$pi)) > 1) {
      plot_opchar$`ESS(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `ESS(pi)`,
                                        colour = Design)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(ESS), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))
      print(plot_opchar$`ESS(pi)`)
      plot_opchar$`VSS(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `VSS(pi)`,
                                        colour = Design)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(VSS), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))

      plot_opchar$`Med(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `Med(pi)`,
                                        colour = Design)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(Med), "(",
                                       pi, ")",
                                       sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))

      plot_opchar$`P(pi)` <- ggplot2::ggplot()
      if (min(opchar$pi) < opchar$des$des$pi0) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
      if (all(min(opchar$pi) <= opchar$des$des$pi1,
              max(opchar$pi) >= opchar$des$des$pi0)) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = amber,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "orange")
      }
      if (max(opchar$pi) > opchar$des$des$pi1) {
        plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
          ggplot2::geom_rect(data = green,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "green4")
      }
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
        ggplot2::geom_hline(yintercept =
                              c(opchar$des$des$alpha,
                                1 - opchar$des$des$beta),
                            linetype = 2) +
        ggplot2::geom_line(data = opchar$opchar,
                           ggplot2::aes(x = pi,
                                        y = `P(pi)`,
                                        colour = Design)) +
        ggplot2::xlab(expression(pi)) +
        ggplot2::ylab(expression(paste(italic(P), "(", pi,
                                       ")", sep = ""))) +
        theme_singlearm() +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                                    limits = c(min(opchar$opchar$pi),
                                               max(opchar$opchar$pi)))
      print(plot_opchar$`P(pi)`)
    } else {
      plot_opchar$`ESS(pi)` <- NULL
      plot_opchar$`VSS(pi)` <- NULL
      plot_opchar$`Med(pi)` <- NULL
      plot_opchar$`P(pi)`   <- NULL
    }
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_opchar = plot_opchar, opchar = opchar))
  }
}

