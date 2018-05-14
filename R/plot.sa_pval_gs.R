#' Plot the p-values, and the p-value calculation procedures performance, in a group sequential
#' single-arm trial design for a single binary endpoint
#'
#' Plots the p-values, and the performance of the p-value calculation procedures,
#' in a group sequential single-arm trial design for a single binary endpoint
#' determined using \code{pval_gs()}. A range of plots are available, of which
#' the p-values and the expected p-value curve will be printed by default.
#'
#' @param est An object of class \code{"sa_pval_gs"}, as returned by \code{pval_gs()}.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the optimal group sequential design for the default parameters
#' des  <- des_gs()
#' # Determine the performance of the p-value calculation procedures for a range of
#' # possible response probabilities
#' pval <- pval_gs(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the p-values and the calculation procedures performance
#' plot(pval)
#' @seealso \code{\link{des_gs}}, \code{\link{opchar_gs}}, \code{\link{est_gs}}, \code{\link{pval_gs}},
#' and \code{\link{ci_gs}}, and their associated \code{plot} family of functions.
#' @export
plot.sa_pval_gs <- function(pval, output = F) {

  ##### Input Checking #########################################################

  check_sa_pval_gs(pval)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_pval <- list()

  new_levels <- levels(pval$pval$method)
  for (i in 1:length(new_levels)) {
    if (new_levels[i] == "conditional") {
      new_levels[i] <- "Conditional"
    } else if (new_levels[i] == "mle") {
      new_levels[i] <- "MLE-ordering"
    } else if (new_levels[i] == "naive") {
      new_levels[i] <- "Naive"
    } else {
      new_levels[i] <- "UMVUE-ordering"
    }
  }
  pval$pval$method <- plyr::mapvalues(pval$pval$method,
                                      from = levels(pval$pval$method),
                                      to = new_levels)

  pval$pval$k <- plyr::mapvalues(pval$pval$k, from = levels(pval$pval$k),
                               to = paste("k =", levels(pval$pval$k)))
  plot_pval$pval <- ggplot2::ggplot(data = pval$pval,
                                    ggplot2::aes(x = s, y = `p(s,m)`,
                                                 color = method)) +
    ggplot2::xlab(expression(italic(s))) +
    ggplot2::ylab(expression(paste(italic(p), "(", italic(s),
                                   ",", italic(n), "|", pi[0],
                                   ")", sep = ""))) +
    ggplot2::geom_hline(yintercept = pval$des$alpha,
                        linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_line() + ggplot2::facet_grid(.~k) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm()
  print(plot_pval$pval)
  if (!is.null(pval$perf)) {

    new_levels <- levels(pval$perf$method)
    for (i in 1:length(new_levels)) {
      if (new_levels[i] == "conditional") {
        new_levels[i] <- "Conditional"
      } else if (new_levels[i] == "mle") {
        new_levels[i] <- "MLE-ordering"
      } else if (new_levels[i] == "naive") {
        new_levels[i] <- "Naive"
      } else {
        new_levels[i] <- "UMVUE-ordering"
      }
    }
    pval$perf$method <- plyr::mapvalues(pval$perf$method,
                                       from = levels(pval$perf$method),
                                       to = new_levels)

    if (min(pval$pi) < pval$des$des$pi0) {
      red   <- tibble::tibble(start = min(pval$pi),
                              end = min(pval$des$des$pi0, max(pval$pi)))
    }
    if (all(min(pval$pi) <= pval$des$des$pi1,
            max(pval$pi) >= pval$des$des$pi0)) {
      amber <- tibble::tibble(start = max(pval$des$des$pi0, min(pval$pi)),
                              end = min(pval$des$des$pi1, max(pval$pi)))
    }
    if (max(pval$pi) > pval$des$des$pi1) {
      green <- tibble::tibble(start = max(pval$des$des$pi1,
                                          min(pval$pi)), end = max(pval$pi))
    }

    plot_pval$`E(p|pi)` <- ggplot2::ggplot()
      if (min(pval$pi) < pval$des$des$pi0) {
        plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
    if (all(min(pval$pi) <= pval$des$des$pi1,
            max(pval$pi) >= pval$des$des$pi0)) {
      plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(pval$pi) > pval$des$des$pi1) {
      plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
      ggplot2::geom_hline(yintercept = pval$des$alpha,
                          linetype = 2) +
      ggplot2::geom_line(data = pval$perf,
                         ggplot2::aes(x = pi,
                                      y = `E(p|pi)`,
                                      color = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste("E(", italic(p),
                                     "|", pi, ")",
                                     sep = ""))) +
      ggthemes::scale_color_ptol("method") +
      ggplot2::theme(legend.title =
                       ggplot2::element_blank(),
                     legend.position = "bottom") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(pval$perf$pi),
                                             max(pval$perf$pi)))
    print(plot_pval$`E(p|pi)`)

    plot_pval$`Var(p|pi)` <- ggplot2::ggplot()
      if (min(pval$pi) < pval$des$des$pi0) {
        plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
          ggplot2::geom_rect(data = red,
                             ggplot2::aes(xmin = start,
                                          xmax = end,
                                          ymin = -Inf,
                                          ymax = Inf),
                             alpha = 0.1,
                             fill = "firebrick2")
      }
    if (all(min(pval$pi) <= pval$des$des$pi1,
            max(pval$pi) >= pval$des$des$pi0)) {
      plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(pval$pi) > pval$des$des$pi1) {
      plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
      ggplot2::geom_line(data = pval$perf,
                         ggplot2::aes(x = pi,
                                      y = `Var(p|pi)`,
                                      color = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Var), "(",
                                     italic(p), "|",
                                     pi, ")",
                                     sep = ""))) +
      ggthemes::scale_color_ptol("method") +
      ggplot2::theme(legend.title =
                       ggplot2::element_blank(),
                     legend.position = "bottom") +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(pval$perf$pi),
                                             max(pval$perf$pi)))
  }

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_pval = plot_pval, pval = pval)
    return(output)
  }
}
