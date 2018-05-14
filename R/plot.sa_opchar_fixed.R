#' Plot the operating characteristics of single-stage single-arm trial designs
#' for a single binary endpoint
#'
#' Plots the operating characteristics of single-stage single-arm trial designs
#' determined using \code{opchar_fixed()}. A range of plots are available, of
#' which the power curve will be printed by default.
#'
#' @param x An object of class \code{"sa_opchar_fixed"}, as returned by
#' \code{opchar_fixed()}.
#' @param ... Included for compatibility with the generic. Not currently used.
#' @param output A logical variable indicating whether the outputs described
#' below should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is
#' returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal single-stage design for the default parameters
#' des    <- des_fixed()
#' # Determine operating characteristics for a range of response probabilities
#' opchar <- opchar_fixed(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the power curve
#' plot(opchar)
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, \code{\link{ci_fixed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_opchar_fixed <- function(x, ..., output = F) {

  opchar <- x

  ##### Input Checking #########################################################

  check_sa_opchar_fixed(opchar)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_opchar <- list()
  pi          <- opchar$pi
  des         <- opchar$des$des
  if (min(pi) < des$pi0) {
    red   <- tibble::tibble(start = min(pi), end = min(des$pi0, max(pi)))
  }
  if (all(min(pi) <= des$pi1,
          max(pi) >= des$pi0)) {
    amber <- tibble::tibble(start = max(des$pi0, min(pi)),
                            end = min(des$pi1, max(pi)))
  }
  if (max(pi) > des$pi1) {
    green <- tibble::tibble(start = max(des$pi1, min(pi)), end = max(pi))
  }
  if (is.null(opchar$add_des)) {
    if (all(des$pi0 %in% pi, !(des$pi1 %in% pi))) {
      pmf_pi0        <- dplyr::filter(opchar$pmf, pi == des$pi0)
      outcome        <- numeric(nrow(pmf_pi0))
      for (i in 1:nrow(pmf_pi0)) {
        if (pmf_pi0$s[i] >= des$r) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi0        <- dplyr::mutate(pmf_pi0,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      plot_opchar$pmf_pi0 <- ggplot2::ggplot(data = pmf_pi0,
                                             ggplot2::aes(x = s,
                                                          y = `f(s,m|pi)`,
                                                          fill = outcome)) +
                               ggplot2::geom_bar(stat = "identity") +
                               ggplot2::xlab(expression(italic(s))) +
                               ggplot2::ylab(expression(paste(italic(f), "(",
                                                              italic(s), ",",
                                                              italic(m), "|",
                                                              pi[0], ")",
                                                              sep = ""))) +
                               ggplot2::scale_fill_manual(values =
                                                            c("firebrick2",
                                                              "green4")) +
                               theme_singlearm()
    } else if (all(!(des$pi0 %in% pi), des$pi1 %in% pi)) {
      pmf_pi1        <- dplyr::filter(opchar$pmf, pi == des$pi1)
      outcome        <- numeric(nrow(pmf_pi1))
      for (i in 1:nrow(pmf_pi1)) {
        if (pmf_pi1$s[i] >= des$r) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi1        <- dplyr::mutate(pmf_pi1,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      plot_opchar$pmf_pi1 <- ggplot2::ggplot(data = pmf_pi1,
                                             ggplot2::aes(x = s,
                                                          y = `f(s,m|pi)`,
                                                          fill = outcome)) +
                               ggplot2::geom_bar(stat = "identity") +
                               ggplot2::xlab(expression(italic(s))) +
                               ggplot2::ylab(expression(paste(italic(f), "(",
                                                              italic(s), ",",
                                                              italic(m), "|",
                                                              pi[1], ")",
                                                              sep = ""))) +
                               ggplot2::scale_fill_manual(values =
                                                            c("firebrick2",
                                                              "green4")) +
                               theme_singlearm()
    } else if (all(des$pi0 %in% pi, des$pi1 %in% pi)) {
      pmf_pi0pi1     <- dplyr::filter(opchar$pmf, pi %in% c(des$pi0, des$pi1))
      outcome        <- numeric(nrow(pmf_pi0pi1))
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$s[i] >= des$r) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi0pi1     <- dplyr::mutate(pmf_pi0pi1,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      facet_pi        <- numeric(nrow(pmf_pi0pi1))
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$pi[i] == des$pi0) {
          facet_pi[i] <- "pi == pi[0]"
        } else {
          facet_pi[i] <- "pi == pi[1]"
        }
      }
      pmf_pi0pi1$pi   <- factor(facet_pi, levels = unique(facet_pi))
      plot_opchar$pmf_pi0pi1 <- ggplot2::ggplot(data = pmf_pi0pi1,
                                                ggplot2::aes(x = s,
                                                             y = `f(s,m|pi)`,
                                                             fill = outcome)) +
                                  ggplot2::geom_bar(stat = "identity") +
                                  ggplot2::xlab(expression(italic(s))) +
                                  ggplot2::ylab(expression(paste(italic(f), "(",
                                                                 italic(s), ",",
                                                                 italic(m), "|",
                                                                 pi, ")",
                                                                 sep = ""))) +
                                  ggplot2::scale_fill_manual(values =
                                                               c("firebrick2",
                                                                 "green4")) +
                                  ggplot2::facet_grid(.~pi,
                                                      labeller =
                                                        ggplot2::label_parsed) +
                                  theme_singlearm()
    }
    plot_opchar$`ESS(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = red,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`ESS(pi)`   <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y =
                                                                   `ESS(pi)`)) +
                                 ggplot2::xlab(expression(pi)) +
                                 ggplot2::ylab(expression(paste(italic(ESS),
                                                                "(", pi, ")",
                                                                sep = ""))) +
                                 theme_singlearm()
    plot_opchar$`VSS(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_rect(data = red,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`VSS(pi)`   <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y =
                                                                   `VSS(pi)`)) +
                                 ggplot2::xlab(expression(pi)) +
                                 ggplot2::ylab(expression(paste(italic(VSS),
                                                                "(", pi, ")",
                                                                sep = ""))) +
                                 theme_singlearm()
    plot_opchar$`Med(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = red,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`Med(pi)`   <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y =
                                                                   `Med(pi)`)) +
                                 ggplot2::xlab(expression(pi)) +
                                 ggplot2::ylab(expression(paste(italic(Med),
                                                                "(", pi, ")",
                                                                sep = ""))) +
                                 theme_singlearm()

    plot_opchar$`P(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = red,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1,
                                                  fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = amber,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1, fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = green,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1, fill = "green4")
    }
    plot_opchar$`P(pi)`   <- plot_opchar$`P(pi)` +
                               ggplot2::geom_hline(yintercept =
                                                     c(opchar$des$des$alpha,
                                                       1 - opchar$des$des$beta),
                                                   linetype = 2) +
                               ggplot2::geom_line(data = opchar$opchar,
                               ggplot2::aes(x = pi, y = `P(pi)`)) +
                               ggplot2::xlab(expression(pi)) +
                               ggplot2::ylab(expression(paste(italic(P), "(",
                                                              pi, ")",
                                                              sep = ""))) +
                               theme_singlearm()
    print(plot_opchar$`P(pi)`)
  } else {
    plot_opchar$`ESS(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = red,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`ESS(pi)` <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`ESS(pi)`   <- plot_opchar$`ESS(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y = `ESS(pi)`,
                                                                 color =
                                                                   Design)) +
                                ggthemes::scale_color_ptol("Design") +
                                ggplot2::xlab(expression(pi)) +
                                ggplot2::ylab(expression(paste(italic(ESS), "(",
                                                               pi, ")",
                                                               sep = ""))) +
                                theme_singlearm()
    plot_opchar$`VSS(pi)`   <- ggplot2::ggplot()
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
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`VSS(pi)` <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`VSS(pi)`   <- plot_opchar$`VSS(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y = `VSS(pi)`,
                                                                 color =
                                                                   Design)) +
                                 ggthemes::scale_color_ptol("Design") +
                                 ggplot2::xlab(expression(pi)) +
                                 ggplot2::ylab(expression(paste(italic(VSS),
                                                                "(", pi, ")",
                                                                sep = ""))) +
                                 theme_singlearm()
    plot_opchar$`Med(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = red,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = amber,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "orange")
    }
    if (max(pi) > des$pi1) {
      plot_opchar$`Med(pi)` <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_rect(data = green,
                                                    ggplot2::aes(xmin = start,
                                                                 xmax = end,
                                                                 ymin = -Inf,
                                                                 ymax = Inf),
                                                    alpha = 0.1,
                                                    fill = "green4")
    }
    plot_opchar$`Med(pi)`   <- plot_opchar$`Med(pi)` +
                                 ggplot2::geom_line(data = opchar$opchar,
                                                    ggplot2::aes(x = pi,
                                                                 y = `Med(pi)`,
                                                                 color =
                                                                   Design)) +
                                 ggthemes::scale_color_ptol("Design") +
                                 ggplot2::xlab(expression(pi)) +
                                 ggplot2::ylab(expression(paste(italic(Med),
                                                                "(", pi, ")",
                                                                sep = ""))) +
                                 theme_singlearm()
    plot_opchar$`P(pi)`   <- ggplot2::ggplot()
    if (min(pi) < des$pi0) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = red,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1,
                                                  fill = "firebrick2")
    }
    if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = amber,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1, fill = "orange")
    }
    if (max(opchar$pi) > opchar$des$des$pi1) {
      plot_opchar$`P(pi)` <- plot_opchar$`P(pi)` +
                               ggplot2::geom_rect(data = green,
                                                  ggplot2::aes(xmin = start,
                                                               xmax = end,
                                                               ymin = -Inf,
                                                               ymax = Inf),
                                                  alpha = 0.1, fill = "green4")
    }
    plot_opchar$`P(pi)`   <- plot_opchar$`P(pi)` +
                               ggplot2::geom_line(data = opchar$opchar,
                                                  ggplot2::aes(x = pi,
                                                               y = `P(pi)`,
                                                               color =
                                                                 Design)) +
                               ggthemes::scale_color_ptol("Design") +
                               ggplot2::xlab(expression(pi)) +
                               ggplot2::ylab(expression(paste(italic(P), "(",
                                                              pi, ")",
                                                              sep = ""))) +
                               theme_singlearm()
    print(plot_opchar$`P(pi)`)
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_opchar = plot_opchar, opchar = opchar))
  }
}
