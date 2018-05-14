#' Plot the second stage sample sizes of two-stage adaptive single-arm trial
#' designs for a single binary endpoint
#'
#' Plots the second stage sample sizes of two-stage adaptive single-arm trial
#' designs determined using \code{des_adaptive()}. The values of
#' \ifelse{html}{\out{<i>n</i><sub>2</sub>(<i>s</i><sub>n<sub>1</sub>
#' </sub>))}}{\eqn{n_2(\tilde{s}_1)}} are plotted.
#'
#' Support is available to simultaneously plot the second stage sample sizes of multiple
#' two-stage adaptive single-arm clinical trial designs for a single binary primary
#' endpoint.
#'
#' @param des An object of class \code{"sa_des_adaptive"}, as returned by \code{des_adaptive()}.
#' @param ... Additional objects of class \code{"sa_des_adaptive"}. These will be grouped
#' in to a list named \code{"add_des"}.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the optimal adaptive two-stage design for the default parameters
#' des    <- des_adaptive()
#' # Plot the stopping boundaries
#' plot(des)
#' # Find the optimal adaptive two-stage design for a 10% type-I error rate
#' des_10 <- des_adaptive(alpha = 0.1)
#' # Plot the second stage sample sizes for both designs
#' plot(des, des_10)
#' @seealso \code{\link{des_adaptive}}, \code{\link{opchar_adaptive}}, and their
#' associated \code{plot} family of functions.
#' @export
plot.sa_des_adaptive <- function(des, ..., output = F) {

  ##### Input Checking #########################################################

  check_sa_des_adaptive(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_adaptive(eval(add_des[[i]]), paste("add_des", i, sep = ""))
    }
    for (i in 1:num_add_des) {
      if (eval(add_des[[i]])$des$pi0 != des$des$pi0) {
        stop("Each supplied design must have been designed for the same value of pi0")
      }
    }
  }
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_des     <- list()
  if (num_add_des == 0) {
    add_des    <- NULL
    int_tibble <- tibble::tibble(s1 = 0:des$des$n1,
                                 n2 = des$des$n2)
    plot_des$`n2(tilde(s)1)` <- ggplot2::ggplot(int_tibble,
                                         ggplot2::aes(x = s1, y = n2)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab(expression(italic(tilde(s)[1]))) +
      ggplot2::ylab(expression(paste(italic(n[2]), "(", italic(tilde(s)[1]), ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(minor_breaks = NULL) +
      ggplot2::scale_y_continuous(minor_breaks = NULL)
    print(plot_des$`n2(tilde(s)1)`)
  } else {
    all_des            <- list()
    all_des[[1]]       <- des
    for (i in 1:num_add_des) {
      all_des[[i + 1]] <- eval(add_des[[i]])
    }
    num_des            <- 1 + num_add_des
    int_tibble         <- NULL
    for (i in 1:num_des) {
      int_tibble <- rbind(int_tibble, cbind(paste("Design", i),
                                            0:all_des[[i]]$des$n1,
                                            all_des[[i]]$des$n2))
    }
    int_tibble <- tibble::tibble(Design = factor(int_tibble[, 1],
                                                 levels = unique(int_tibble[, 1])),
                                 s1 = as.integer(int_tibble[, 2]),
                                 n2 = as.integer(int_tibble[, 3]))

    plot_des$`n2(tilde(s)1)` <- ggplot2::ggplot(int_tibble,
                                         ggplot2::aes(x = s1, y = n2,
                                                      colour = Design)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab(expression(italic(tilde(s)[1]))) +
      ggplot2::ylab(expression(paste(italic(n[2]), "(", italic(tilde(s)[1]), ")",
                                     sep = ""))) +
      theme_singlearm() +
      ggplot2::scale_x_continuous(minor_breaks = NULL) +
      ggplot2::scale_y_continuous(minor_breaks = NULL)
    print(plot_des$`n2(tilde(s)1)`)
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(plot_des = plot_des, des = des, add_des = add_des))
  }
}
