#' @export
plot.sa_opchar_curtailed <- function(opchar, output = F) {

  ##### Input Checking #########################################################

  #check_sa_opchar_curtailed(opchar)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_opchar <- list()
  if (is.null(opchar$add_des)) {
    if (all(opchar$des$des$pi0 %in% opchar$pmf$pi,
            !opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi0        <- dplyr::filter(opchar$pmf, pi == opchar$des$des$pi0)
      outcome        <- numeric(nrow(pmf_pi0))
      for (i in 1:nrow(pmf_pi0)) {
        if (pmf_pi0$s[i] >= opchar$des$des$r_curt[as.numeric(pmf_pi0$k[i])]) {
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
        ggplot2::theme(legend.title =
                         ggplot2::element_blank(),
                       legend.position = "bottom") + ggplot2::facet_grid(.~k)
      print(plot_opchar$`f(s,m|pi)`)
    } else if (all(!opchar$des$des$pi0 %in% opchar$pmf$pi,
                   opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi1        <- dplyr::filter(opchar$pmf, pi == opchar$des$des$pi1)
      outcome        <- numeric(nrow(pmf_pi1))
      for (i in 1:nrow(pmf_pi1)) {
        if (pmf_pi1$s[i] >= opchar$des$des$r_curt[as.numeric(pmf_pi1$k[i])]) {
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
        ggplot2::theme(legend.title =
                         ggplot2::element_blank(),
                       legend.position = "bottom") + ggplot2::facet_grid(.~k)
      print(plot_opchar$`f(s,m|pi)`)
    } else if (all(opchar$des$des$pi0 %in% opchar$pmf$pi,
                   opchar$des$des$pi1 %in% opchar$pmf$pi)) {
      pmf_pi0pi1     <- dplyr::filter(opchar$pmf, pi %in% c(opchar$des$des$pi0,
                                                            opchar$des$des$pi1))
      outcome        <- numeric(nrow(pmf_pi0pi1))
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$s[i] >=
            opchar$des$des$r_curt[as.numeric(pmf_pi0pi1$k[i])]) {
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
      for (i in 1:nrow(pmf_pi0pi1)) {
        if (pmf_pi0pi1$pi[i] == opchar$des$des$pi0) {
          facet_pi[i] <- "pi == pi[0]"
        } else {
          facet_pi[i] <- "pi == pi[1]"
        }
      }
      pmf_pi0pi1$pi   <- factor(facet_pi, levels = unique(facet_pi))

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
        ggplot2::theme(legend.title =
                         ggplot2::element_blank(),
                       legend.position = "bottom") +
        ggplot2::facet_grid(.~pi,
                            labeller = ggplot2::label_parsed) + ggplot2::facet_grid(.~k)
      print(plot_opchar$`f(s,m|pi)`)
    }
  }

  pi          <- opchar$pi
  des         <- opchar$des$des
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

  if (is.null(opchar$add_des)) {
    if (length(unique(opchar$pmf$pi)) > 1) {
      plot_opchar$`ESS(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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
      suppressWarnings(print(plot_opchar$`ESS(pi)`))

      plot_opchar$`VSS(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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

      plot_opchar$`Med(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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

      plot_opchar$`P(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1, fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1, fill = "green4") +
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
      suppressWarnings(print(plot_opchar$`P(pi)`))
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
      ggplot2::geom_bar(stat="identity") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab("Stopping probability") +
      theme_singlearm()

    int_tibble            <- tidyr::gather(opchar$opchar, "key", "value",
                                           (6 + 2*opchar$des$des$J):
                                             (5 + 3*opchar$des$des$J))
    plot_opchar$stopping  <- ggplot2::ggplot(int_tibble,
                                             ggplot2::aes(pi, value,
                                                          fill = key)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab("Stopping probability") +
      theme_singlearm()
  } else {
    if (length(unique(opchar$pmf$pi)) > 1) {
      plot_opchar$`ESS(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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
        ggplot2::guides(colour =
                          ggplot2::guide_legend(nrow =
                                                  ceiling(length(unique(opchar$opchar$Design))/2)))
      suppressWarnings(print(plot_opchar$`ESS(pi)`))

      plot_opchar$`VSS(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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
        ggplot2::guides(colour =
                          ggplot2::guide_legend(nrow =
                                                  ceiling(length(unique(opchar$opchar$Design))/2)))

      plot_opchar$`Med(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4") +
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
        ggplot2::guides(colour =
                          ggplot2::guide_legend(nrow =
                                                  ceiling(length(unique(opchar$opchar$Design))/2)))

      plot_opchar$`P(pi)` <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2") +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1, fill = "orange") +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1, fill = "green4") +
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
        ggplot2::guides(colour =
                          ggplot2::guide_legend(nrow =
                                                  ceiling(length(unique(opchar$opchar$Design))/2)))
      suppressWarnings(print(plot_opchar$`P(pi)`))
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
