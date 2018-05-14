#' Plot the stopping boundaries of single-stage single-arm trial designs for a
#' single binary endpoint
#'
#' Plots the stopping boundaries of single-stage single-arm trial designs
#' determined using \code{des_fixed()}. The possible
#' \ifelse{html}{\out{(<i>s</i>,<i>m</i>)}}{\eqn{(s,m)}} states during the trial
#' are plotted in a colour coded manner, to indicate their associated decision
#' rules.
#'
#' Support is available to simultaneously plot the stopping boundaries of
#' multiple applicable designs, using faceting.
#'
#' @param des An object of class \code{"sa_des_fixed"}, as returned by
#' \code{des_fixed()}.
#' @param ... Additional objects of class \code{"sa_des_fixed"}. These will be
#' grouped in to a list named \code{"add_des"}.
#' @param output A logical variable indicating whether the outputs described
#' below should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is
#' returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item A tibble in the slot \code{$states} containing details of the possible
#' states for each of the designs, and their associated decision rules.
#' \item Each of the input variables as specified, subject to internal
#' modification.
#' }
#' @examples
#' # Find the optimal single-stage design for the default parameters
#' des    <- des_fixed()
#' # Plot the stopping boundaries
#' plot(des)
#' # Find the optimal single-stage design for a 10% type-I error rate
#' des_10 <- des_fixed(alpha = 0.1)
#' # Plot the stopping boundaries for both designs
#' plot(des, des_10)
#' @seealso \code{\link{des_fixed}}, \code{\link{opchar_fixed}},
#' \code{\link{est_fixed}}, \code{\link{pval_fixed}}, \code{\link{ci_fixed}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_des_fixed <- function(des, ..., output = F) {

  ##### Input Checking #########################################################

  check_sa_des_fixed(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_fixed(eval(add_des[[i]]), paste("add_des[[", i, "]]",
                                                   sep = ""))
    }
    for (i in 1:num_add_des) {
      if (eval(add_des[[i]])$des$pi0 != des$des$pi0) {
        stop("Each supplied design must have been designed for the same value of pi0")
      }
    }
  }
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_des <- list()
  if (num_add_des == 0) {
    states <- iterpc::getall(iterpc::iterpc(des$des$n + 1, 2, 0:des$des$n, F,
                                            T))[-1, ]
    states <- tibble::tibble(s = states[, 1], m = states[, 2])
    states <- dplyr::mutate(states,
                            status = factor(ifelse(s <= des$des$a &
                                                     m == des$des$n,
                                            "Do not reject",
                                            ifelse(s >= des$des$r &
                                                    m == des$des$n, "Reject",
                                                   "Continue")),
                                            c("Continue", "Do not reject",
                                              "Reject")))
    plot_des$states <- ggplot2::ggplot(states, ggplot2::aes(x = m, y = s,
                                                            colour = status,
                                                            shape = status)) +
                         theme_singlearm() +
                         ggplot2::geom_point() +
                         ggplot2::xlab(expression(italic(m))) +
                         ggplot2::ylab(expression(italic(s))) +
                         ggplot2::scale_x_continuous(minor_breaks = NULL) +
                         ggplot2::scale_y_continuous(minor_breaks = NULL) +
                         ggplot2::scale_color_manual(values = c("gray50",
                                                                "firebrick2",
                                                                "green4")) +
                         ggplot2::scale_shape_manual(values = c(1, 4, 3))
    print(plot_des$states)
    add_des <- NULL
  } else {
    all_des            <- list()
    all_des[[1]]       <- des
    for (i in 1:num_add_des) {
      all_des[[i + 1]] <- eval(add_des[[i]])
    }
    states             <- list()
    for (i in 1:(num_add_des + 1)) {
      states[[i]] <- iterpc::getall(iterpc::iterpc(all_des[[i]]$des$n + 1, 2,
                                                0:all_des[[i]]$des$n, F,
                                                T))[-1, ]
      states[[i]] <- tibble::tibble(Design = paste("Design ", i, ": (",
                                                   all_des[[i]]$des$a, ", ",
                                                   all_des[[i]]$des$r, ")/",
                                                   all_des[[i]]$des$n,
                                                   sep = ""),
                                    s = states[[i]][, 1], m = states[[i]][, 2])
      states[[i]] <-
        dplyr::mutate(states[[i]],
                      status = factor(ifelse(s <= all_des[[i]]$des$a &
                                               m == all_des[[i]]$des$n,
                                             "Do not reject",
                                             ifelse(s >= all_des[[i]]$des$r &
                                                      m == all_des[[i]]$des$n,
                                                    "Reject", "Continue")),
                                      c("Continue", "Do not reject", "Reject")))
    }
    states <- tibble::as_tibble(plyr::rbind.fill(states))
    plot_des$states <- ggplot2::ggplot(states,
                                       ggplot2::aes(x = m, y = s,
                                                    colour = status,
                                                    shape = status)) +
                         theme_singlearm() +
                         ggplot2::geom_point(size = 2/(num_add_des + 1)) +
                         ggplot2::xlab(expression(italic(m))) +
                         ggplot2::ylab(expression(italic(s))) +
                         ggplot2::scale_color_manual(values = c("gray50",
                                                                "firebrick2",
                                                                "green4")) +
                         ggplot2::scale_x_continuous(minor_breaks = NULL) +
                         ggplot2::scale_y_continuous(minor_breaks = NULL) +
                         ggplot2::facet_wrap(~Design) +
                         ggplot2::scale_shape_manual(values = c(1, 4, 3))
    print(plot_des$states)
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_des = plot_des, states = states, des = des,
                add_des = add_des))
  }
}
