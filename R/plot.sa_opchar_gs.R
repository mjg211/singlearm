#' Plot the operating characteristics of group sequential single-arm trial designs for a
#' single binary endpoint
#'
#' Plots the operating characteristics of group sequential single-arm trial designs
#' determined using \code{opchar_gs()}. A range of plots are available, of which
#' the expected sample size and power curves will be printed by default.
#'
#' @param x An object of class \code{"sa_opchar_gs"}, as returned by \code{opchar_gs()}.
#' @param ... Included for compatibility with the generic. Not currently used.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the optimal group sequential design for the default parameters
#' des         <- des_gs()
#' # Determine operating characteristics for a range of response probabilities
#' opchar      <- opchar_gs(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the expected sample size and power curve
#' plot(opchar)
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}}, \code{\link{est_gs}}, \code{\link{pval_gs}},
#' and \code{\link{ci_gs}}, and their associated \code{plot} family of functions.
#' @export
plot.sa_opchar_gs <- function(x, ..., output = F) {

  opchar <- x

  ##### Input Checking #########################################################

  check_sa_opchar_gs(opchar)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_opchar <- list()
  if (is.null(opchar$add_des)) {
    if (all(opchar$des$des$pi0 %in% opchar$pmf$pi,
            !opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi0        <- dplyr::filter(opchar$pmf, pi == opchar$des$des$pi0)
      outcome        <- numeric(nrow(pmf_pi0))
      for (i in 1:nrow(pmf_pi0)) {
        if (pmf_pi0$s[i] >= opchar$des$des$r[as.numeric(pmf_pi0$k[i])]) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi0$k      <- plyr::mapvalues(pmf_pi0$k, from = levels(pmf_pi0$k),
                                           to = paste("k =",
                                                      levels(pmf_pi0$k)))
      pmf_pi0        <- dplyr::mutate(pmf_pi0,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      plot_opchar$`f(s,m|pi)` <- ggplot2::ggplot(data = pmf_pi0,
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
        ggplot2::scale_fill_manual(values = c("firebrick2",
                                              "green4")) +
        theme_singlearm() + ggplot2::facet_grid(.~k)
      print(plot_opchar$`f(s,m|pi)`)
    } else if (all(!opchar$des$des$pi0 %in% opchar$pmf$pi,
                   opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi1        <- dplyr::filter(opchar$pmf, pi == opchar$des$des$pi1)
      outcome        <- numeric(nrow(pmf_pi1))
      for (i in 1:nrow(pmf_pi1)) {
        if (pmf_pi1$s[i] >= opchar$des$des$r[as.numeric(pmf_pi1$k[i])]) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi1        <- dplyr::mutate(pmf_pi1,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      pmf_pi1$k      <- plyr::mapvalues(pmf_pi1$k, from = levels(pmf_pi1$k),
                                        to = paste("k =",
                                                   levels(pmf_pi1$k)))
      plot_opchar$`f(s,m|pi)` <- ggplot2::ggplot(data = pmf_pi1,
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
        ggplot2::scale_fill_manual(values = c("firebrick2",
                                              "green4")) +
        theme_singlearm() + ggplot2::facet_grid(.~k)
      print(plot_opchar$`f(s,m|pi)`)
    } else if (all(opchar$des$des$pi0 %in% opchar$pmf$pi,
                   opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi0pi1     <- dplyr::filter(opchar$pmf, pi %in% c(opchar$des$des$pi0,
                                                            opchar$des$des$pi1))
      outcome        <- numeric(nrow(pmf_pi0pi1))
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$s[i] >=
              opchar$des$des$r[as.numeric(pmf_pi0pi1$k[i])]) {
          outcome[i] <- "Reject"
        } else {
          outcome[i] <- "Accept"
        }
      }
      pmf_pi0pi1     <- dplyr::mutate(pmf_pi0pi1,
                                      outcome = factor(outcome,
                                                       levels = c("Accept",
                                                                  "Reject")))
      pmf_pi0pi1$k   <- plyr::mapvalues(pmf_pi0pi1$k,
                                        from = levels(pmf_pi0pi1$k),
                                        to = paste("k =",
                                                   levels(pmf_pi0pi1$k)))
      facet_pi        <- numeric(nrow(pmf_pi0pi1))
      facet_k         <- numeric(nrow(pmf_pi0pi1))
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$pi[i] == opchar$des$des$pi0) {
          facet_pi[i] <- "pi == pi[0]"
        } else {
          facet_pi[i] <- "pi == pi[1]"
        }
        local    <- strsplit(as.character(pmf_pi0pi1$k[i]), " ")[[1]]
        local[1] <- "italic(k)"
        local[2] <- "=="
        facet_k[i] <- paste(local, collapse = " ")
      }
      pmf_pi0pi1$pi   <- factor(facet_pi, levels = unique(facet_pi))
      pmf_pi0pi1$k    <- factor(facet_k, levels = unique(facet_k))

      plot_opchar$`f(s,m|pi)` <- ggplot2::ggplot(data = pmf_pi0pi1,
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
        theme_singlearm() +
        ggplot2::facet_grid(pi~k,
                            labeller = ggplot2::label_parsed)
      print(plot_opchar$`f(s,m|pi)`)
    }
  }

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
    int_tibble            <- tidyr::gather(opchar$opchar, "key", "value",
                                           6:(5 + 2*opchar$des$des$J))
    plot_opchar$rejection <- ggplot2::ggplot(int_tibble,
                                             ggplot2::aes(pi, value,
                                                          fill = key)) +
                               ggplot2::geom_area() +
                               ggplot2::xlab(expression(pi)) +
                               ggplot2::ylab("Stopping probability") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(opchar$opchar$pi),
                                             max(opchar$opchar$pi))) +
      ggplot2::scale_y_continuous(expand = c(0, 0))

    int_tibble            <- tidyr::gather(opchar$opchar, "key", "value",
                                           (6 + 2*opchar$des$des$J):
                                             (5 + 3*opchar$des$des$J))
    plot_opchar$stopping  <- ggplot2::ggplot(int_tibble,
                                             ggplot2::aes(pi, value,
                                                          fill = key)) +
                               ggplot2::geom_area() +
                               ggplot2::xlab(expression(pi)) +
                               ggplot2::ylab("Stopping probability") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(opchar$opchar$pi),
                                             max(opchar$opchar$pi))) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
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
                                               max(opchar$opchar$pi))) +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
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
    output <- list(plot_opchar = plot_opchar, opchar = opchar)
    return(output)
  }

}
