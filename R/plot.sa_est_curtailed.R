#' Plot the point estimates, and the point estimation procedures performance, in
#' a curtailed group sequential single-arm trial design for a single binary
#' endpoint
#'
#' Plots the point estimates, and the performance of the point estimation
#' procedure, in a curtailed group sequential single-arm trial design for a
#' single binary endpoint determined using \code{est_curtailed()}. A range of
#' plots are available, of which the bias curve will be printed by default.
#'
#' @param x An object of class \code{"sa_est_curtailed"}, as returned by
#' \code{est_curtailed()}.
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
#' # Find the optimal curtailed group sequential design for the default
#' parameters
#' des <- des_curtailed()
#' # Determine the performance of the point estimation procedure for a range of
#' # possible response probabilities
#' est <- est_curtailed(des, pi = seq(0, 1, 0.01))
#' # Plot the point estimates and the estimation procedure performance
#' plot(est)
#' @seealso \code{\link{des_curtailed}}, \code{\link{opchar_curtailed}},
#' \code{\link{est_curtailed}}, \code{\link{pval_curtailed}},
#' \code{\link{ci_curtailed}}, and their associated \code{plot} family of
#' functions.
#' @export
plot.sa_est_curtailed <- function(x, ..., output = F) {

  est <- x

  ##### Input Checking #########################################################

  check_sa_est_curtailed(est)
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_est     <- list()
  if (!is.null(est$perf)) {
    new_levels <- levels(est$perf$method)
    for (i in 1:length(new_levels)) {
      if (new_levels[i] == "bias_adj") {
        new_levels[i] <- "Bias adjusted"
      } else if (new_levels[i] == "bias_sub") {
        new_levels[i] <- "Bias subtracted"
      } else if (new_levels[i] == "conditional") {
        new_levels[i] <- "Conditional"
      } else if (new_levels[i] == "naive") {
        new_levels[i] <- "Naive"
      } else if (new_levels[i] == "mue") {
        new_levels[i] <- "MUE"
      } else {
        new_levels[i] <- "UMVUE"
      }
    }
    est$perf$method <- plyr::mapvalues(est$perf$method,
                                       from = levels(est$perf$method),
                                       to = new_levels)

    if (min(est$pi) < est$des$des$pi0) {
      red   <- tibble::tibble(start = min(est$pi),
                              end = min(est$des$des$pi0, max(est$pi)))
    }
    if (all(min(est$pi) <= est$des$des$pi1,
            max(est$pi) >= est$des$des$pi0)) {
      amber <- tibble::tibble(start = max(est$des$des$pi0, min(est$pi)),
                              end = min(est$des$des$pi1, max(est$pi)))
    }
    if (max(est$pi) > est$des$des$pi1) {
      green <- tibble::tibble(start = max(est$des$des$pi1,
                                          min(est$pi)), end = max(est$pi))
    }

    plot_est$`E(hat(pi)|pi)` <- ggplot2::ggplot()
    if (min(est$pi) < est$des$des$pi0) {
      plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(est$pi) <= est$des$des$pi1,
            max(est$pi) >= est$des$des$pi0)) {
      plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(est$pi) > est$des$des$pi1) {
      plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_est$`E(hat(pi)|pi)` <- plot_est$`E(hat(pi)|pi)` +
      ggplot2::geom_abline(slope = 1,
                           intercept = 0,
                           linetype = 2) +
      ggplot2::geom_hline(yintercept = est$des$pi0,
                          linetype = 2) +
      ggplot2::geom_hline(yintercept = est$des$pi1,
                          linetype = 2) +
      ggplot2::geom_line(data = est$perf,
                         ggplot2::aes(x = pi,
                                      y = `E(hat(pi)|pi)`, colour = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste("E(", hat(pi),
                                     "|", pi, ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(est$perf$pi),
                                             max(est$perf$pi)))

    plot_est$`Var(hat(pi)|pi)` <- ggplot2::ggplot()
    if (min(est$pi) < est$des$des$pi0) {
      plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(est$pi) <= est$des$des$pi1,
            max(est$pi) >= est$des$des$pi0)) {
      plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(est$pi) > est$des$des$pi1) {
      plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_est$`Var(hat(pi)|pi)` <- plot_est$`Var(hat(pi)|pi)` +
      ggplot2::geom_line(data = est$perf,
                         ggplot2::aes(x = pi,
                                      y = `Var(hat(pi)|pi)`, colour = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Var), "(",
                                     hat(pi), "|",
                                     pi, ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(est$perf$pi),
                                             max(est$perf$pi)))

    plot_est$`Bias(hat(pi)|pi)` <- ggplot2::ggplot()
    if (min(est$pi) < est$des$des$pi0) {
      plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(est$pi) <= est$des$des$pi1,
            max(est$pi) >= est$des$des$pi0)) {
      plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(est$pi) > est$des$des$pi1) {
      plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_est$`Bias(hat(pi)|pi)` <- plot_est$`Bias(hat(pi)|pi)` +
      ggplot2::geom_hline(yintercept = 0,
                          linetype = 2) +
      ggplot2::geom_line(data = est$perf,
                         ggplot2::aes(x = pi,
                                      y = `Bias(hat(pi)|pi)`, colour = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(Bias), "(",
                                     hat(pi),
                                     "|", pi,
                                     ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(est$perf$pi),
                                             max(est$perf$pi)))
    print(plot_est$`Bias(hat(pi)|pi)`)

    plot_est$`RMSE(hat(pi)|pi)` <- ggplot2::ggplot()
    if (min(est$pi) < est$des$des$pi0) {
      plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
        ggplot2::geom_rect(data = red,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "firebrick2")
    }
    if (all(min(est$pi) <= est$des$des$pi1,
            max(est$pi) >= est$des$des$pi0)) {
      plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
        ggplot2::geom_rect(data = amber,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "orange")
    }
    if (max(est$pi) > est$des$des$pi1) {
      plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
        ggplot2::geom_rect(data = green,
                           ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = -Inf,
                                        ymax = Inf),
                           alpha = 0.1,
                           fill = "green4")
    }
    plot_est$`RMSE(hat(pi)|pi)` <- plot_est$`RMSE(hat(pi)|pi)` +
      ggplot2::geom_line(data = est$perf,
                         ggplot2::aes(x = pi,
                                      y = `RMSE(hat(pi)|pi)`, colour = method)) +
      ggplot2::xlab(expression(pi)) +
      ggplot2::ylab(expression(paste(italic(RMSE), "(",
                                     hat(pi),
                                     "|", pi,
                                     ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(min(est$perf$pi),
                                             max(est$perf$pi)))

  }

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_est = plot_est, est = est)
    return(output)
  }
}
