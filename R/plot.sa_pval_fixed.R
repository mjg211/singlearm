#' Plot the p-values, and the p-value calculation procedures performance, in a
#' single-stage single-arm trial design for a single binary endpoint
#'
#' Plots the p-values, and the performance of the p-value calculation
#' procedures, in a single-stage single-arm trial design for a single binary
#' endpoint determined using \code{pval_fixed()}. A range of plots are
#' available, of which the p-values and the expected p-value curve will be
#' printed by default.
#'
#' @param pval An object of class \code{"sa_pval_fixed"}, as returned by
#' \code{pval_fixed()}.
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
#' des  <- des_fixed()
#' # Determine the performance of the p-value calculation procedures for a range
#' # of possible response probabilities
#' pval <- pval_fixed(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the p-values and the calculation procedures performance
#' plot(pval)
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, \code{\link{ci_fixed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_pval_fixed <- function(pval, output = F) {

  ##### Input Checking #########################################################

  check_sa_pval_fixed(pval)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_pval  <- list()
  new_levels <- levels(pval$pval$method)
  for (i in 1:length(new_levels)) {
    split         <- strsplit(new_levels[i], split = "")[[1]]
    split[1]      <- toupper(split[1])
    new_levels[i] <- paste(split, collapse = "")
  }
  pval$pval$method <- plyr::mapvalues(pval$pval$method,
                                      from = levels(pval$pval$method),
                                      to = new_levels)

  plot_pval$pval <- ggplot2::ggplot() +
                      ggplot2::xlab(expression(italic(s))) +
                      ggplot2::ylab(expression(paste(italic(p), "(", italic(s),
                                                     ",", italic(n), "|", pi[0],
                                                     ")", sep = ""))) +
                      ggplot2::geom_hline(yintercept = pval$des$alpha,
                                          linetype = 2) +
                      ggplot2::geom_point(data = pval$pval,
                                          ggplot2::aes(x = s, y = `p(s,m)`,
                                                       color = method)) +
                      ggplot2::geom_line(data = pval$pval,
                                         ggplot2::aes(x = s, y = `p(s,m)`,
                                                      color = method)) +
                      ggthemes::scale_color_ptol("method") +
    theme_singlearm()
  print(plot_pval$pval)

  pval$perf$method <- plyr::mapvalues(pval$perf$method,
                                      from = levels(pval$perf$method),
                                      to = new_levels)
  pi      <- pval$pi
  des     <- pval$des$des
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

  plot_pval$`E(p|pi)`   <- ggplot2::ggplot()
  if (min(pi) < des$pi0) {
    plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_pval$`E(p|pi)` <- plot_pval$`E(p|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1, fill = "green4")
  }
  plot_pval$`E(p|pi)`   <- plot_pval$`E(p|pi)` +
    ggplot2::geom_hline(yintercept = des$alpha,
                        linetype = 2) +
    ggplot2::geom_line(data = pval$perf,
                       ggplot2::aes(x = pi,
                                    y = `E(p|pi)`,
                                    color =
                                      method)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(E), "(", italic(p),
                                   "|", pi, ")",
                                   sep = ""))) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm()
  print(plot_pval$`E(p|pi)`)

  plot_pval$`Var(p|pi)`   <- ggplot2::ggplot()
  if (min(pi) < des$pi0) {
    plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_pval$`Var(p|pi)` <- plot_pval$`Var(p|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin = start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "green4")
  }
  plot_pval$`Var(p|pi)`   <- plot_pval$`Var(p|pi)` +
    ggplot2::geom_line(data = pval$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `Var(p|pi)`,
                                    color =
                                      method)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(Var),
                                   "(", italic(p),
                                   "|", pi, ")",
                                   sep = ""))) +
    ggthemes::scale_color_ptol("method") +
    theme_singlearm()

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_pval = plot_pval, pval = pval)
    return(output)
  }
}
