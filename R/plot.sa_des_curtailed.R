#' Plot the stopping boundaries of curtailed single-arm trial designs for a
#' single binary endpoint
#'
#' Plots the stopping boundaries of group sequential single-arm trial designs
#' determined using \code{des_gs()}. The possible
#' \ifelse{html}{\out{(<i>s</i>,<i>m</i>)}}{\eqn{(s,m)}} states during the trial
#' are plotted in a colour coded manner to indicate their associated decision
#' rules.
#'
#' Support is available to simultaneously plot the stopping boundaries of multiple
#' group sequential single-arm clinical trial designs for a single binary primary
#' endpoint, using faceting.
#'
#' @param x An object of class \code{"sa_des_gs"}, as returned by \code{des_gs()}.
#' @param ... Additional objects of class \code{"sa_des_gs"}. These will be grouped
#' in to a list named \code{"add_des"}.
#' @param output A logical variable indicating whether the outputs described below
#' should be returned.
#' @return If \code{output = TRUE}, a list containing the following elements is returned
#' \itemize{
#' \item A list in the slot \code{$plot_des} containing the available plots.
#' \item Each of the input variables as specified, subject to internal modification.
#' }
#' @examples
#' # Find the optimal two-stage design for the default parameters with NSC
#' des         <- des_curtailed()
#' # Plot the stopping boundaries
#' plot(des)
#' # Find the optimal two-stage design for a 10% type-I error rate with NSC
#' des_10      <- des_curtailed(alpha = 0.1)
#' # Plot the stopping boundaries for both designs
#' plot(des, des_10)
#' @seealso \code{\link{des_curtailed}}, \code{\link{opchar_curtailed}}, and
#' their associated \code{plot} family of functions.
#' @export
plot.sa_des_curtailed <- function(x, ..., output = F) {

  des <- x

  ##### Input Checking #########################################################

  check_sa_des_curtailed(des, "des")
  add_des     <- pryr::named_dots(...)
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    for (i in 1:num_add_des) {
      check_sa_des_curtailed(eval(add_des[[i]]), paste("add_des", i, sep = ""))
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
    J          <- des$des$J_curt
    a          <- des$des$a_curt
    r          <- des$des$r_curt
    n          <- des$des$n_curt
    states     <- tibble::as.tibble(expand.grid(s = 0:n[1],
                                                m = 1:n[1]))
    states     <- dplyr::filter(states, s <= m)
    states     <- dplyr::mutate(states,
                                status = ifelse(s <= a[1] & m == n[1],
                                                "Do not reject",
                                                ifelse(s >= r[1] &
                                                         m == n[1], "Reject",
                                                       "Continue")))
    cont       <- c(max(0, a[1] + 1), min(r[1] - 1, n[1]))
    for (j in 2:J) {
      vals_j     <- tibble::as.tibble(expand.grid(s = 0:n[j], m = 1:n[j]))
      vals_j     <- dplyr::filter(vals_j, s <= m)
      states_j   <- NULL
      for (sj in seq(from = cont[1], to = cont[2], by = 1)) {
        states_j <- rbind(states_j, dplyr::mutate(vals_j,  s = s + sj,
                                                  m = m +
                                                    cumsum(n[1:(j - 1)])[j - 1]))
      }
      states_j   <- dplyr::mutate(states_j,
                                  status = ifelse(s <= a[j] &
                                                    m == cumsum(n[1:j])[j],
                                                  "Do not reject",
                                                  ifelse(s >= r[j] &
                                                           m == cumsum(n[1:j])[j],
                                                         "Reject", "Continue")))
      cont       <- c(min(states_j$s[states_j$s > a[j]]),
                      max(states_j$s[states_j$s < r[j]]))
      states     <- rbind(states, states_j)
    }
    states$status <- factor(states$status, levels = c("Continue", "Do not reject",
                                                      "Reject"))
    plot_des$states <- ggplot2::ggplot(states, ggplot2::aes(x = m, y = s,
                                                            colour = status,
                                                            shape = status)) +
      ggplot2::geom_point() +
      ggplot2::xlab(expression(italic(m))) +
      ggplot2::ylab(expression(italic(s[m]))) +
      ggplot2::scale_color_manual(values = c("gray50",
                                             "firebrick2",
                                             "green4")) +
      ggplot2::scale_shape_manual(values = c(1, 4, 3)) +
      theme_singlearm()
    print(plot_des$states)
  } else {
    all_des            <- list()
    all_des[[1]]       <- des
    for (i in 1:num_add_des) {
      all_des[[i + 1]] <- eval(add_des[[i]])
    }
    num_des            <- 1 + num_add_des
    Js   <- NULL
    for (i in 1:num_des) {
      Js <- c(Js, all_des[[i]]$des$J_curt)
    }
    all_states         <- NULL
    for (i in 1:num_des) {
      states <- tibble::as.tibble(expand.grid(s = 0:all_des[[i]]$des$n_curt[1],
                                              m = 1:all_des[[i]]$des$n_curt[1]))
      states <- dplyr::filter(states, s <= m)
      states <- dplyr::mutate(states,
                              status = ifelse(s <= all_des[[i]]$des$a_curt[1] &
                                                m == all_des[[i]]$des$n_curt[1],
                                              "Do not reject",
                                              ifelse(s >= all_des[[i]]$des$r_curt[1] &
                                                       m == all_des[[i]]$des$n_curt[1],
                                                     "Reject", "Continue")))
      cont   <- c(max(0, all_des[[i]]$des$a_curt[1] + 1),
                  min(all_des[[i]]$des$r_curt[1] - 1, all_des[[i]]$des$n_curt[1]))
      if (Js[i] > 1) {
        for (j in 2:Js[i]) {
          vals_j     <- tibble::as.tibble(expand.grid(s = 0:all_des[[i]]$des$n_curt[j],
                                                      m = 1:all_des[[i]]$des$n_curt[j]))
          vals_j     <- dplyr::filter(vals_j, s <= m)
          states_j   <- NULL
          for (sj in seq(from = cont[1], to = cont[2], by = 1)) {
            states_j <- rbind(states_j, dplyr::mutate(vals_j,  s = s + sj,
                                                      m = m +
                                                        cumsum(all_des[[i]]$des$n_curt[1:(j - 1)])[j - 1]))
          }
          states_j   <- dplyr::mutate(states_j,
                                      status = ifelse(s <= all_des[[i]]$des$a_curt[j] &
                                                        m == cumsum(all_des[[i]]$des$n_curt[1:j])[j],
                                                      "Do not reject",
                                                      ifelse(s >= all_des[[i]]$des$r_curt[j] &
                                                               m == cumsum(all_des[[i]]$des$n_curt[1:j])[j],
                                                             "Reject", "Continue")))
          cont       <- c(min(states_j$s[states_j$s > all_des[[i]]$des$a_curt[j]]),
                          max(states_j$s[states_j$s < all_des[[i]]$des$r_curt[j]]))
          states     <- rbind(states, states_j)
        }
      }
      states$status <- factor(states$status,
                              levels = c("Continue", "Do not reject", "Reject"))
      all_states    <- rbind(all_states, cbind(paste("Design", i), states))
    }
    colnames(all_states) <- c("Design", "s", "m", "status")
    all_states           <- tibble::as_tibble(all_states)
    Js   <- NULL
    for (i in 1:num_des) {
      Js <- c(Js, all_des[[i]]$des$J)
    }
    new_levels      <- levels(all_states$Design)
    for (i in 1:num_des) {
      if (Js[i] > 1) {
        new_levels[i] <- paste(new_levels[i], ": ",
                               paste("(", all_des[[i]]$des$a[1:(Js[i] - 1)], ",",
                                     all_des[[i]]$des$r[1:(Js[i] - 1)], ")/",
                                     cumsum(all_des[[i]]$des$n[1:(Js[i] - 1)]), sep = "",
                                     collapse = ", "), ", ",
                               all_des[[i]]$des$a[Js[i]], "/", sum(all_des[[i]]$des$n), ", ",
                               paste("thetaF = (", all_des[[i]]$des$thetaF[1:(Js[i] - 1)], ",",
                                     sep = "", collapse = ", "),
                               all_des[[i]]$des$thetaF[Js[i]], "), thetaE = (",
                               paste(all_des[[i]]$des$thetaE[1:(Js[i] - 1)], ",",
                                     sep = "", collapse = ", "),
                               all_des[[i]]$des$thetaE[Js[i]], ")",
                               sep = "", collapse = "")
      } else {
        new_levels[i] <- paste(new_levels[i], ": ", all_des[[i]]$des$a, "/",
                               all_des[[i]]$des$n, ", thetaF = ", all_des[[i]]$des$thetaF,
                               ", thetaE = ", all_des[[i]]$des$thetaE)
      }
    }
    all_states$Design   <- plyr::mapvalues(all_states$Design,
                                           from = levels(all_states$Design),
                                           to = new_levels)
    plot_des$states <- ggplot2::ggplot(all_states, ggplot2::aes(x = m, y = s,
                                                                colour = status,
                                                                shape = status)) +
      ggplot2::geom_point(size = 1) +
      ggplot2::xlab(expression(italic(m))) +
      ggplot2::ylab(expression(italic(s[m]))) +
      ggplot2::scale_color_manual(values = c("gray50",
                                             "firebrick2",
                                             "green4")) +
      ggplot2::facet_grid(.~Design) +
      ggplot2::scale_x_continuous(minor_breaks = NULL) +
      ggplot2::scale_y_continuous(minor_breaks = NULL) +
      ggplot2::scale_shape_manual(values = c(1, 4, 3)) +
      theme_singlearm() +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 6))
    print(plot_des$states)

  }


  ##### Outputting #############################################################

  if (output) {
    return(list(plot_des = plot_des, des = des, add_des = add_des))
  }
}
