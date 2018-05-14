#' Plot the confidence intervals, and the confidence interval determination
#' procedures performance, in a single-stage single-arm trial design for a
#' single binary endpoint
#'
#' Plots the confidence intervals, and the performance of the confidence
#' interval determination procedures, in a single-stage single-arm trial design
#' for a single binary endpoint determined using \code{ci_fixed()}. A range of
#' plots are available, of which the confidence intervals and the coverage
#' probability curve will be printed by default.
#'
#' @param ci An object of class \code{"sa_ci_fixed"}, as returned by
#' \code{ci_fixed()}.
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
#' des <- des_fixed()
#' # Determine the performance of the confidence interval determination
#' # procedures for a range of possible response probabilities
#' ci  <- ci_fixed(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the confidence intervals and the determination procedures performance
#' plot(ci)
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, \code{\link{ci_fixed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_ci_fixed <- function(ci, output = F) {

  ##### Input Checking #########################################################

  check_sa_ci_fixed(ci)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_ci    <- list()
  new_levels <- levels(ci$ci$method)
  for (i in 1:length(new_levels)) {
    if (new_levels[i] == "agresti_coull") {
      new_levels[i] <- "Agresti-Coull"
    } else if (new_levels[i] == "clopper_pearson") {
      new_levels[i] <- "Clopper-Pearson"
    } else if (new_levels[i] == "jeffreys") {
      new_levels[i] <- "Jeffreys"
    } else if (new_levels[i] == "mid_p") {
      new_levels[i] <- "Mid-p"
    } else if (new_levels[i] == "wald") {
      new_levels[i] <- "Wald"
    } else {
      new_levels[i] <- "Wilson"
    }
  }
  ci$ci$method <- plyr::mapvalues(ci$ci$method, from = levels(ci$ci$method),
                                  to = new_levels)
  local_ci   <- tidyr::gather(ci$ci, key = "limit", value = "c",
                              `clow(s,m)`:`cupp(s,m)`)
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
                                      ggplot2::aes(x = s, y = c,
                                                   color = method)) +
                  ggplot2::geom_point(data = dplyr::filter(local_ci,
                                                           limit == "cupp(s,m)"),
                                      ggplot2::aes(x = s, y = c,
                                                   color = method)) +
                  ggplot2::geom_line(data = dplyr::filter(local_ci,
                                                          limit == "clow(s,m)"),
                                     ggplot2::aes(x = s, y = c,
                                                  color = method)) +
                  ggplot2::geom_line(data = dplyr::filter(local_ci,
                                                          limit == "cupp(s,m)"),
                                     ggplot2::aes(x = s, y = c,
                                                  color = method)) +
                  ggthemes::scale_color_ptol("method") +
    theme_singlearm()
  print(plot_ci$ci)
  plot_ci$l <- ggplot2::ggplot() +
                 ggplot2::geom_point(data = ci$ci,
                                     ggplot2::aes(x = s, y = `l(s,m)`,
                                                  color = method)) +
                 ggplot2::geom_line(data = ci$ci,
                                    ggplot2::aes(x = s, y = `l(s,m)`,
                                                 color = method)) +
                 ggplot2::xlab(expression(italic(s))) +
                 ggplot2::ylab(expression(italic(l), "(", italic(s), ",",
                                          italic(n), ")")) +
                 ggthemes::scale_color_ptol("method") +
    theme_singlearm()

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
    theme_singlearm()   +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = 45,
                                           hjust = 1))

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
    theme_singlearm()
  ggplot2::theme(legend.position = "none",
                 axis.text.x =
                   ggplot2::element_text(angle = 45,
                                         hjust = 1))

  pi      <- ci$pi
  des     <- ci$des$des
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

  plot_ci$`E(L|pi)`     <- ggplot2::ggplot()
  if (min(pi) < des$pi0) {
    plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "green4")
  }
  plot_ci$`E(L|pi)` <- plot_ci$`E(L|pi)` +
    ggplot2::geom_line(data = ci$perf,
                       ggplot2::aes(x = pi,
                                    y = `E(L|pi)`,
                                    color = method)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(E), "(", italic(L), "|",
                                   pi, ")", sep = ""))) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm()

  plot_ci$`Var(L|pi)`   <- ggplot2::ggplot()
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
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "orange")
  }
  if (max(ci$pi) > ci$des$des$pi1) {
    plot_ci$`Var(L|pi)` <- plot_ci$`Var(L|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "green4")
  }
  plot_ci$`Var(L|pi)`   <- plot_ci$`Var(L|pi)` +
    ggplot2::geom_line(data = ci$perf,
                       ggplot2::aes(x = pi,
                                    y = `Var(L|pi)`,
                                    color = method)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(Var), "(",
                                   italic(L), "|", pi,
                                   ")", sep = ""))) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm()

  plot_ci$`Cover(C|pi)`   <- ggplot2::ggplot()
  if (min(pi) < des$pi0) {
    plot_ci$`Cover(C|pi)` <- plot_ci$`Cover(C|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
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
  plot_ci$`Cover(C|pi)`   <- plot_ci$`Cover(C|pi)` +
    ggplot2::geom_hline(yintercept =
                          1 - ci$des$alpha,
                        linetype = 2) +
    ggplot2::geom_line(data = ci$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `Cover(C|pi)`,
                                    color =
                                      method)) +
    ggthemes::scale_color_ptol("method") +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(Cover), "(", italic(C), "|",
                                   pi, ")", sep = ""))) +

    theme_singlearm()
  print(plot_ci$`Cover(C|pi)`)

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_ci = plot_ci, ci = ci)
    return(output)
  }
}
