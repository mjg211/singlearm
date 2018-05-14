#' Plot the point estimates, and the point estimation procedures performance, in
#' a single-stage single-arm trial design for a single binary endpoint
#'
#' Plots the point estimates, and the performance of the point estimation
#' procedure, in a single-stage single-arm trial design for a single binary
#' endpoint determined using \code{est_fixed()}. A range of plots are available,
#' of which the point estimates and the root mean squared error curve will be
#' printed by default.
#'
#' @param est An object of class \code{"sa_est_fixed"}, as returned by
#' \code{est_fixed()}.
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
#' # Determine the performance of the point estimation procedure for a range of
#' # possible response probabilities
#' est <- est_fixed(des, pi = seq(from = 0, to = 1, by = 0.01))
#' # Plot the point estimates and the estimation procedure performance
#' plot(est)
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, \code{\link{ci_fixed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_est_fixed <- function(est, output = F) {

  ##### Input Checking #########################################################

  check_sa_est_fixed(est)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_est     <- list()
  plot_est$est <- ggplot2::ggplot() +
                    ggplot2::geom_hline(yintercept = est$des$pi0,
                                        linetype = 2) +
                    ggplot2::geom_hline(yintercept = est$des$pi1,
                                        linetype = 2) +
                    ggplot2::geom_point(data = est$est,
                                        ggplot2::aes(x = s, y = `hat(pi)(s,m)`)) +
                    ggplot2::geom_line(data = est$est,
                                       ggplot2::aes(x = s, y = `hat(pi)(s,m)`)) +
                    ggplot2::xlab(expression(italic(s))) +
                    ggplot2::ylab(expression(paste(hat(pi), "(", italic(s), ",",
                                                   italic(n), ")", sep = ""))) +
                    theme_singlearm()
  print(plot_est$est)
  pi      <- est$pi
  des     <- est$des$des
  if (min(pi) < des$pi0) {
    red   <- tibble::tibble(start = min(pi),
                            end = min(des$pi0, max(pi)))
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    amber <- tibble::tibble(start = max(des$pi0, min(pi)),
                            end = min(des$pi1, max(pi)))
  }
  if (max(pi) > des$pi1) {
    green <- tibble::tibble(start = max(des$pi1,
                                        min(pi)), end = max(pi))
  }
  plot_est$`E(hat(pi)|pi)`   <- ggplot2::ggplot()
  if (min(est$pi) < est$des$des$pi0) {
    plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "orange")
  }
  if (max(est$pi) > est$des$des$pi1) {
    plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax = end,
                                      ymin = -Inf,
                                      ymax = Inf),
                         alpha = 0.1,
                         fill = "green4")
  }
  plot_est$`E(hat(pi)|pi)`   <- plot_est$`E(hat(pi)|pi)` +
    ggplot2::geom_abline(slope = 1,
                         intercept = 0,
                         linetype = 2) +
    ggplot2::geom_hline(yintercept = des$pi0,
                        linetype = 2) +
    ggplot2::geom_hline(yintercept = des$pi1,
                        linetype = 2) +
    ggplot2::geom_line(data = est$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `E(hat(pi)|pi)`)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(E), "(", hat(pi),
                                   "|", pi, ")",
                                   sep = ""))) +
    theme_singlearm()
  plot_est$`Var(hat(pi)|pi)`     <- ggplot2::ggplot()
  if (min(est$pi) < est$des$des$pi0) {
    plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "green4")
  }
  plot_est$`Var(hat(pi)|pi)`   <- plot_est$`Var(hat(pi)|pi)` +
    ggplot2::geom_line(data = est$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `Var(hat(pi)|pi)`)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(Var),
                                   "(",
                                   hat(pi),
                                   "|", pi,
                                   ")",
                                   sep =
                                     ""))) +
    theme_singlearm()

  plot_est$`Bias(hat(pi)|pi)`   <- ggplot2::ggplot()
  if (min(est$pi) < est$des$des$pi0) {
    plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "green4")
  }
  plot_est$`Bias(hat(pi)|pi)`   <- plot_est$`Bias(hat(pi)|pi)` +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = 2) +
    ggplot2::geom_line(data = est$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `Bias(hat(pi)|pi)`)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(Bias),
                                   "(",
                                   hat(pi),
                                   "|", pi,
                                   ")",
                                   sep =
                                     ""))) +
    theme_singlearm()

  plot_est$`RMSE(hat(pi)|pi)`     <- ggplot2::ggplot()
  if (min(pi) < des$pi0) {
    plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
      ggplot2::geom_rect(data = red,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "firebrick2")
  }
  if (all(min(pi) <= des$pi1, max(pi) >= des$pi0)) {
    plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
      ggplot2::geom_rect(data = amber,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "orange")
  }
  if (max(pi) > des$pi1) {
    plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
      ggplot2::geom_rect(data = green,
                         ggplot2::aes(xmin =
                                        start,
                                      xmax =
                                        end,
                                      ymin =
                                        -Inf,
                                      ymax =
                                        Inf),
                         alpha = 0.1,
                         fill = "green4")
  }
  plot_est$`RMSE(hat(pi)|pi)`   <- plot_est$`RMSE(hat(pi)|pi)` +
    ggplot2::geom_line(data = est$perf,
                       ggplot2::aes(x = pi,
                                    y =
                                      `RMSE(hat(pi)|pi)`)) +
    ggplot2::xlab(expression(pi)) +
    ggplot2::ylab(expression(paste(italic(RMSE),
                                   "(",
                                   hat(pi),
                                   "|", pi,
                                   ")",
                                   sep =
                                     ""))) +
    theme_singlearm()
  print(plot_est$`RMSE(hat(pi)|pi)`)

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_est = plot_est, est = est)
    return(output)
  }
}
