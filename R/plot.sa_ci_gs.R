#' Plot the confidence intervals, and the confidence interval determination procedures
#' performance, in a group sequential single-arm trial design for a single binary endpoint
#'
#' Plots the confidence intervals, and the performance of the confidence interval determination procedures,
#' in a group sequential single-arm trial design for a single binary endpoint
#' determined using \code{ci_gs()}. A range of plots are available, of which
#' the confidence intervals and the coverage probability curve will be printed by default.
#'
#' @param x An object of class \code{"sa_ci_gs"}, as returned by \code{ci_gs()}.
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
#' des <- des_gs()
#' # Determine the performance of the confidence interval determination procedures for a range of
#' # possible response probabilities
#' ci  <- ci_gs(des, pi = seq(0, 1, 0.01))
#' # Plot the confidence intervals and the determination procedures performance
#' plot(ci)
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}}, \code{\link{est_gs}},
#' \code{\link{pval_gs}}, \code{\link{ci_gs}}, and their associated \code{plot}
#' family of functions.
#' @export
plot.sa_ci_gs <- function(x, ..., output = F) {

  ci <- x

  ##### Input Checking #########################################################

  check_sa_ci_gs(ci)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_ci    <- list()
  new_levels <- levels(ci$ci$method)
  for (i in 1:length(new_levels)) {
    if (new_levels[i] == "exact") {
      new_levels[i] <- "Exact"
    } else if (new_levels[i] == "mid_p") {
      new_levels[i] <- "Mid-p"
    } else {
      new_levels[i] <- "Naive"
    }
  }
  ci$ci$method <- plyr::mapvalues(ci$ci$method,
                                      from = levels(ci$ci$method),
                                      to = new_levels)

  local_ci <- tidyr::gather(ci$ci, key = "limit", value = "pi",
                            `clow(s,m)`:`cupp(s,m)`)

  local_ci$k <- plyr::mapvalues(ci$ci$k, from = levels(ci$ci$k),
                             to = paste("k =", levels(ci$ci$k)))

  plot_ci$ci <- ggplot2::ggplot() +
    ggplot2::xlab(expression(italic(s))) +
    ggplot2::ylab(expression(paste(italic(c)[low], "(", italic(s),
                                   ",", italic(n), "),  ",
                                   italic(c)[upp], "(", italic(s),
                                   ",", italic(n), ")",
                                   sep = ""))) +
    ggplot2::geom_hline(yintercept = c(ci$des$pi0, ci$des$pi1),
                        linetype = 2) +
    ggplot2::geom_point(data = dplyr::filter(local_ci,
                                             limit == "clow(s,m)"),
                        ggplot2::aes(x = s, y = pi,
                                     color = method)) +
    ggplot2::geom_point(data = dplyr::filter(local_ci,
                                             limit == "cupp(s,m)"),
                        ggplot2::aes(x = s, y = pi,
                                     color = method)) +
    ggplot2::geom_line(data = dplyr::filter(local_ci,
                                            limit == "clow(s,m)"),
                       ggplot2::aes(x = s, y = pi,
                                    color = method)) +
    ggplot2::geom_line(data = dplyr::filter(local_ci,
                                            limit == "cupp(s,m)"),
                       ggplot2::aes(x = s, y = pi,
                                    color = method)) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm() + ggplot2::facet_grid(.~k)
  print(plot_ci$ci)

  plot_ci$l <- ggplot2::ggplot() +
    ggplot2::geom_point(data = ci$ci,
                        ggplot2::aes(x = s, y = `l(s,m)`,
                                     color = method)) +
    ggplot2::geom_line(data = ci$ci,
                       ggplot2::aes(x = s, y = `l(s,m)`,
                                    color = method)) +
    ggplot2::xlab(expression(italic(s))) +
    ggplot2::ylab(expression(italic(l), "(", italic(s), ",", italic(m), ")")) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm() + ggplot2::facet_grid(.~k)

  if (is.null(ci$perf)) {
    perf <- tibble::tibble(method = factor(unique(ci$ci$method),
                                           levels = unique(ci$ci$method)),
                           `bar(L)` = NA, `max(L)` = NA)
    for (i in 1:nrow(perf)) {
      ci_i             <- dplyr::filter(ci$ci, method == ci$ci$method[i])
      perf$`bar(L)`[i] <- mean(ci_i$`l(s,m)`)
      perf$`max(L)`[i] <- max(ci_i$`l(s,m)`)
    }

    plot_ci$`bar(L)` <- ggplot2::ggplot(data = perf,
                                        ggplot2::aes(x = method, y = `bar(L)`,
                                                     fill = method)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Method") +
      ggplot2::ylab(expression(bar(italic(L)))) +
      ggthemes::scale_fill_ptol("method") +
      theme_singlearm() +
      ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = 45,
                                             hjust = 1),
                     legend.position = "none")

    plot_ci$`max(L)` <- ggplot2::ggplot(data = perf,
                                        ggplot2::aes(x = method, y = `max(L)`,
                                                     fill = method)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Method") +
      ggplot2::ylab(expression(paste(italic(max), "(",
                                     italic(L), ")",
                                     sep = ""))) +
      ggthemes::scale_fill_ptol("method") +
      theme_singlearm() +
      ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = 45,
                                             hjust = 1),
                     legend.position = "none")
  } else {

    new_levels <- levels(ci$perf$method)
    for (i in 1:length(new_levels)) {
      if (new_levels[i] == "exact") {
        new_levels[i] <- "Exact"
      } else if (new_levels[i] == "mid_p") {
        new_levels[i] <- "Mid-p"
      } else {
        new_levels[i] <- "Naive"
      }
    }
    ci$perf$method <- plyr::mapvalues(ci$perf$method,
                                    from = levels(ci$perf$method),
                                    to = new_levels)

    plot_ci$`bar(L)` <- ggplot2::ggplot(data = dplyr::filter(ci$perf,
                                                             pi ==
                                                               ci$perf$pi[1]),
                                        ggplot2::aes(x = method,
                                                     y = `bar(L)`,
                                                     fill = method)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Method") +
      ggplot2::ylab(expression(bar(italic(L)))) +
      ggthemes::scale_fill_ptol("method") +
      theme_singlearm() +
      ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = 45,
                                             hjust = 1),
                     legend.position = "none")

    plot_ci$`max(L)` <- ggplot2::ggplot(data = dplyr::filter(ci$perf,
                                                             pi ==
                                                               ci$perf$pi[1]),
                                        ggplot2::aes(x = method,
                                                     y = `max(L)`,
                                                     fill = method)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Method") +
      ggplot2::ylab(expression(paste(italic(max), "(",
                                     italic(L), ")",
                                     sep = ""))) +
      ggthemes::scale_fill_ptol("method") +
      theme_singlearm() +
      ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = 45,
                                             hjust = 1),
                     legend.position = "none")

    if (min(ci$pi) < ci$des$des$pi0) {
      red   <- tibble::tibble(start = min(ci$pi),
                              end = min(ci$des$des$pi0, max(ci$pi)))
    }
    if (all(min(ci$pi) <= ci$des$des$pi1,
            max(ci$pi) >= ci$des$des$pi0)) {
      amber <- tibble::tibble(start = max(ci$des$des$pi0, min(ci$pi)),
                              end = min(ci$des$des$pi1, max(ci$pi)))
    }
    if (max(ci$pi) > ci$des$des$pi1) {
      green <- tibble::tibble(start = max(ci$des$des$pi1,
                                          min(ci$pi)), end = max(ci$pi))
    }

    plot_ci$`E(L|pi)` <- ggplot2::ggplot()
      if (min(ci$pi) < ci$des$des$pi0) {
        plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
    if (all(min(ci$pi) <= ci$des$des$pi1,
            max(ci$pi) >= ci$des$des$pi0)) {
      plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(ci$pi) > ci$des$des$pi1) {
      plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
      ggplot2::geom_line(data = ci$perf,
                         ggplot2::aes(x = pi,
                                      y = `E(L|pi)`,
                                      color = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste("E(", italic(L), "|",
                                     pi, ")", sep = ""))) +
      ggthemes::scale_color_ptol("method") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(ci$perf$pi),
                                             max(ci$perf$pi)))

    plot_ci$`Var(L|pi)` <- ggplot2::ggplot()
      if (min(ci$pi) < ci$des$des$pi0) {
        plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
    if (all(min(ci$pi) <= ci$des$des$pi1,
            max(ci$pi) >= ci$des$des$pi0)) {
      plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(ci$pi) > ci$des$des$pi1) {
      plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
      ggplot2::geom_line(data = ci$perf,
                         ggplot2::aes(x = pi,
                                      y = `Var(L|pi)`,
                                      color = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Var), "(",
                                     italic(L), "|", pi,
                                     ")", sep = ""))) +
      ggthemes::scale_color_ptol("method") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(ci$perf$pi),
                                             max(ci$perf$pi)))

    plot_ci$`Cover(C|pi)` <- ggplot2::ggplot()
      if (min(ci$pi) < ci$des$des$pi0) {
        plot_ci$`Cover(C|pi)` <- plot_ci$`Cover(C|pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
    if (all(min(ci$pi) <= ci$des$des$pi1,
            max(ci$pi) >= ci$des$des$pi0)) {
      plot_ci$`Cover(C|pi)` <- plot_ci$`Cover(C|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(ci$pi) > ci$des$des$pi1) {
      plot_ci$`Cover(C|pi)` <- plot_ci$`Cover(C|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_ci$`Cover(C|pi)` <- plot_ci$`Cover(C|pi)` +
      ggplot2::geom_hline(yintercept =
                            1 - ci$des$alpha,
                          linetype = 2) +
      ggplot2::geom_line(data = ci$perf,
                         ggplot2::aes(x = pi,
                                      y = `Cover(C|pi)`,
                                      color = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Cover), "(",
                                     italic(C), "|", pi,
                                     ")", sep = ""))) +
      ggthemes::scale_color_ptol("method") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(ci$perf$pi),
                                             max(ci$perf$pi)))
    print(plot_ci$`Cover(C|pi)`)
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_ci = plot_ci, ci = ci))
  }
}
