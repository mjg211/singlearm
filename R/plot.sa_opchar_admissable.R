#' Plot the operating characteristics of admissable group sequential single-arm
#' trial designs for a single binary endpoint
#'
#' Plots the operating characteristics of admissable group sequential single-arm trial designs
#' determined using \code{opchar_admissable()}. A range of plots are available, of which
#' the expected sample size and power curves will be printed by default.
#'
#' @param opchar An object of class \code{"sa_opchar_admissable"}, as returned by \code{opchar_admissable()}.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the admissable group sequential design for the default parameters
#' des         <- des_admissable()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_admissable(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the expected sample size and power curve
#' plot(opchar)
#' @seealso \code{\link{des_admissable}}, \code{\link{opchar_admissable}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_opchar_admissable <- function(opchar, output = F) {

  ##### Input Checking #########################################################

  check_sa_opchar_admissable(opchar)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_opchar <- list()
  num_des     <- length(levels(opchar$opchar$Design))

  pi          <- opchar$pi
  des         <- opchar$des$des[[1]]
  if (min(pi) < des$pi0) {
    red   <- tibble::tibble(start = min(pi),
                            end = min(des$pi0, max(pi)))
  }
  if (all(min(pi) <= des$pi1,
          max(pi) >= des$pi0)) {
    amber <- tibble::tibble(start = max(des$pi0, min(pi)),
                            end = min(des$pi1, max(pi)))
  }
  if (max(pi) > des$pi1) {
    green <- tibble::tibble(start = max(des$pi1,
                                        min(pi)), end = max(pi))
  }

  if (length(unique(opchar$pmf$pi)) > 1) {

    plot_opchar$`ESS(pi)` <- ggplot2::ggplot()
    if (min(opchar$pi) < opchar$des$des[[1]]$pi0) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(opchar$pi) <= opchar$des$des[[1]]$pi1,
            max(opchar$pi) >= opchar$des$des[[1]]$pi0)) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(opchar$pi) > opchar$des$des[[1]]$pi1) {
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
                                      color = Design)) +
      ggthemes::scale_color_ptol("Design") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(ESS), "(",
                                     pi, ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(opchar$opchar$pi),
                                                               max(opchar$opchar$pi))) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = ceiling(num_des/3), byrow = T))
    print(plot_opchar$`ESS(pi)`)
    plot_opchar$`VSS(pi)` <- ggplot2::ggplot()
    if (min(opchar$pi) < opchar$des$des[[1]]$pi0) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(opchar$pi) <= opchar$des$des[[1]]$pi1,
            max(opchar$pi) >= opchar$des$des[[1]]$pi0)) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(opchar$pi) > opchar$des$des[[1]]$pi1) {
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
                                      color = Design)) +
      ggthemes::scale_color_ptol("Design") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(VSS), "(",
                                     pi, ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(opchar$opchar$pi),
                                                               max(opchar$opchar$pi))) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = ceiling(num_des/3), byrow = T))

    plot_opchar$`Med(pi)` <- ggplot2::ggplot()
    if (min(opchar$pi) < opchar$des$des[[1]]$pi0) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(opchar$pi) <= opchar$des$des[[1]]$pi1,
            max(opchar$pi) >= opchar$des$des[[1]]$pi0)) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(opchar$pi) > opchar$des$des[[1]]$pi1) {
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
                                      color = Design)) +
      ggthemes::scale_color_ptol("Design") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Med), "(",
                                     pi, ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(opchar$opchar$pi),
                                                               max(opchar$opchar$pi))) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = ceiling(num_des/3), byrow = T))
    plot_opchar$`P(pi)` <- ggplot2::ggplot()
    if (min(opchar$pi) < opchar$des$des[[1]]$pi0) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(opchar$pi) <= opchar$des$des[[1]]$pi1,
            max(opchar$pi) >= opchar$des$des[[1]]$pi0)) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(opchar$pi) > opchar$des$des[[1]]$pi1) {
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
      ggplot2::geom_line(data = opchar$opchar,
                         ggplot2::aes(x = pi,
                                      y = `P(pi)`,
                                      color = Design)) +
      ggthemes::scale_color_ptol("Design") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(P), "(", pi,
                                     ")", sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(opchar$opchar$pi),
                                                               max(opchar$opchar$pi))) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = ceiling(num_des/3), byrow = T))
    print(plot_opchar$`P(pi)`)

  } else {
    plot_opchar$`ESS(pi)` <- NULL
    plot_opchar$`VSS(pi)` <- NULL
    plot_opchar$`Med(pi)` <- NULL
    plot_opchar$`P(pi)`   <- NULL
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_opchar = plot_opchar, opchar = opchar))
  }
}
