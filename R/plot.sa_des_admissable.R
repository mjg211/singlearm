#' Plot the stopping boundaries of admissable group sequential single-arm trial
#' designs for a single binary endpoint
#'
#' Plots the stopping boundaries of admissable group sequential single-arm trial
#' designs determined using \code{des_admissable()}. The possible
#' \ifelse{html}{\out{(<i>s</i>,<i>m</i>)}}{\eqn{(s,m)}} states during the trial
#' are plotted in a colour coded manner to indicate their associated decision
#' rules.
#'
#' In addition, admissable design plots are produced to depict in which region
#' of the weight-space each of the designs is optimal, and to illuminate the
#' expected sample size characteristics of the admissable designs.
#'
#' @param x An object of class \code{"sa_des_admissable"}, as returned by
#' \code{des_admissable()}.
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
#' # Find the admissable two-stage designs for the default parameters
#' des <- des_admissable()
#' # Plot the stopping boundaries and admissable design plot
#' plot(des)
#' @seealso \code{\link{des_admissable}}, \code{\link{opchar_admissable}},
#' and their associated \code{plot} family of functions.
#' @export
plot.sa_des_admissable <- function(x, ..., output = F) {

  x <- des

  ##### Input Checking #########################################################

  check_sa_des_admissable(des, "des")
  check_logical(output, "output")

  ##### Main Computations ######################################################

  plot_des       <- list()
  num_des        <- length(des$des)
  all_des        <- list()
  for (i in 1:num_des) {
    all_des[[i]] <- des$des[[i]]
  }
  all_states     <- NULL
  J              <- all_des[[1]]$J
  for (i in 1:num_des) {
    states <- tibble::as.tibble(expand.grid(s = 0:all_des[[i]]$n[1],
                                            m = 1:all_des[[i]]$n[1]))
    states <- dplyr::filter(states, s <= m)
    states <- dplyr::mutate(states,
                            status = ifelse(s <= all_des[[i]]$a[1] &
                                              m == all_des[[i]]$n[1],
                                            "Do not reject",
                                            ifelse(s >= all_des[[i]]$r[1] &
                                                     m == all_des[[i]]$n[1],
                                                   "Reject", "Continue")))
    cont   <- c(max(0, all_des[[i]]$a[1] + 1),
                min(all_des[[i]]$r[1] - 1, all_des[[i]]$n[1]))
    for (j in 2:J) {
      vals_j     <- tibble::as.tibble(expand.grid(s = 0:all_des[[i]]$n[j],
                                                  m = 1:all_des[[i]]$n[j]))
      vals_j     <- dplyr::filter(vals_j, s <= m)
      states_j   <- NULL
      for (sj in seq(from = cont[1], to = cont[2], by = 1)) {
        states_j <- rbind(states_j, dplyr::mutate(vals_j,  s = s + sj,
                                                  m = m +
                                                    cumsum(all_des[[i]]$n[1:(j - 1)])[j - 1]))
      }
      states_j   <- dplyr::mutate(states_j,
                                  status = ifelse(s <= all_des[[i]]$a[j] &
                                                    m == cumsum(all_des[[i]]$n[1:j])[j],
                                                  "Do not reject",
                                                  ifelse(s >= all_des[[i]]$r[j] &
                                                           m == cumsum(all_des[[i]]$n[1:j])[j],
                                                         "Reject", "Continue")))
      cont       <- c(min(states_j$s[states_j$s > all_des[[i]]$a[j]]),
                      max(states_j$s[states_j$s < all_des[[i]]$r[j]]))
      states     <- rbind(states, states_j)
    }
    states$status <- factor(states$status,
                            levels = c("Continue", "Do not reject", "Reject"))
    all_states    <- rbind(all_states, cbind(des$admissable[i, 1], states))
  }
  colnames(all_states) <- c("Design", "s", "m", "status")
  all_states           <- tibble::as_tibble(all_states)
  plot_des$states <- ggplot2::ggplot(all_states, ggplot2::aes(x = m, y = s,
                                                              colour = status,
                                                              shape = status)) +
    ggplot2::geom_point(size = 2/num_des) +
    ggplot2::xlab(expression(italic(m))) +
    ggplot2::ylab(expression(italic(s[m]))) +
    ggplot2::scale_color_manual(values = c("gray50",
                                           "firebrick2",
                                           "green4")) +
    ggplot2::scale_shape_manual(values = c(1, 4, 3)) +
    theme_singlearm() +
    ggplot2::facet_wrap(~Design)
  print(plot_des$states)

  plot_des$admissable_tri <- ggplot2::ggplot(des$weights, ggplot2::aes(x = w0, y = w1,
                                                              colour = Design)) +
    ggplot2::geom_point(shape = 20) +
    ggplot2::xlab(expression(italic(w[0]))) +
    ggplot2::ylab(expression(italic(w[1]))) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = ceiling(num_des/3), byrow = T)) +
    ggplot2::scale_shape_manual(values = c(1, 4, 3)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme_singlearm()
  print(plot_des$admissable_tri)

  n_range         <- sort(unique(des$feasible$`max(N)`))
  admissable_null <- des$feasible[1:length(n_range), ]
  admissable_alt  <- des$feasible[1:length(n_range), ]
  check_null      <- rep("Not admissable", length(n_range))
  check_alt       <- rep("Not admissable", length(n_range))
  for (n in 1:length(n_range)) {
    feasible_n           <- dplyr::filter(des$feasible, `max(N)` == n_range[n])
    admissable_null[n, ] <- feasible_n[which.min(feasible_n$`ESS(pi0)`), ]
    admissable_alt[n, ]  <- feasible_n[which.min(feasible_n$`ESS(pi1)`), ]
    if (!is.na(row_match(admissable_null[n, ], des$admissable[, -1]))) {
      check_null[n] <- "Admissable"
    }
    if (!is.na(row_match(admissable_alt[n, ], des$admissable[, -1]))) {
      check_alt[n]  <- "Admissable"
    }
  }
  if (which.min(admissable_null$`ESS(pi0)`) != 1) {
    check_null[which.min(admissable_null$`ESS(pi0)`)] <- "Null-optimal"
    check_null[1] <- "Minimax"
  } else {
    check_null[1] <- "Minimax/Null-optimal"
  }
  if (which.min(admissable_alt$`ESS(pi1)`) != 1) {
    check_alt[which.min(admissable_alt$`ESS(pi1)`)] <- "Alt-optimal"
    check_alt[1] <- "Minimax"
  } else {
    check_alt[1] <- "Minimax/Alt-optimal"
  }
  designs_null <- numeric(length(n_range))
  designs_alt  <- numeric(length(n_range))
  if (des$J == 2) {
    if (des$efficacy) {
      for (i in 1:length(n_range)) {
        designs_null[i] <- paste(paste(paste("(", admissable_null[i, ]$a1, ",",
                                             admissable_null[i, ]$r1, ")/",
                                             admissable_null[i, ]$n1, sep = "",
                                             collapse = ", "), sep = "", collapse = ""), ", ",
                                 admissable_null[i, ]$a2, "/",
                                 admissable_null[i, ]$`max(N)`, sep = "", collapse = "")
        designs_alt[i]  <- paste(paste(paste("(", admissable_alt[i, ]$a1, ",",
                                             admissable_alt[i, ]$r1, ")/",
                                             admissable_alt[i, ]$n1, sep = "",
                                             collapse = ", "), sep = "", collapse = ""), ", ",
                                 admissable_alt[i, ]$a2, "/",
                                 admissable_alt[i, ]$`max(N)`, sep = "", collapse = "")
      }
    } else {
      for (i in 1:length(n_range)) {
        designs_null[i] <- paste(paste(paste(admissable_null[i, ]$a1, "/",
                                             admissable_null[i, ]$n1, sep = "",
                                             collapse = ", "), sep = "", collapse = ""), ", ",
                                 admissable_null[i, ]$a2, "/",
                                 admissable_null[i, ]$`max(N)`, sep = "", collapse = "")
        designs_alt[i]  <- paste(paste(paste(admissable_alt[i, ]$a1, "/",
                                             admissable_alt[i, ]$n1, sep = "",
                                             collapse = ", "), sep = "", collapse = ""), ", ",
                                 admissable_alt[i, ]$a2, "/",
                                 admissable_alt[i, ]$`max(N)`, sep = "", collapse = "")
      }
    }
  } else {
    for (i in 1:length(n_range)) {
      designs_null[i] <- paste(paste(paste("(", c(admissable_null[i, ]$a1,
                                                  admissable_null[i, ]$a2), ",",
                                           c(admissable_null[i, ]$r1,
                                             admissable_null[i, ]$r2), ")/",
                                           c(admissable_null[i, ]$n1,
                                             admissable_null[i, ]$n1 +
                                               admissable_null[i, ]$n2), sep = "",
                                           collapse = ", "), sep = "", collapse = ""), ", ",
                               admissable_null[i, ]$a3, "/",
                               admissable_null[i, ]$`max(N)`, sep = "", collapse = "")
      designs_alt[i] <- paste(paste(paste("(", c(admissable_alt[i, ]$a1,
                                                  admissable_alt[i, ]$a2), ",",
                                           c(admissable_alt[i, ]$r1,
                                             admissable_alt[i, ]$r2), ")/",
                                           c(admissable_alt[i, ]$n1,
                                             admissable_alt[i, ]$n1 +
                                               admissable_alt[i, ]$n2), sep = "",
                                           collapse = ", "), sep = "", collapse = ""), ", ",
                               admissable_alt[i, ]$a3, "/",
                               admissable_alt[i, ]$`max(N)`, sep = "", collapse = "")
    }
  }
  admissable_null <- tibble::as_tibble(cbind(designs_null, check_null, admissable_null))
  admissable_alt  <- tibble::as_tibble(cbind(designs_alt, check_alt, admissable_alt))


  plot_des$admissable_pi0 <- ggplot2::ggplot() +
    ggplot2::geom_line(data = admissable_null, ggplot2::aes(x = `max(N)`,
                                                            y = `ESS(pi0)`)) +
    ggplot2::geom_point(data = dplyr::filter(admissable_null,
                                     check_null != "Not admissable"),
                       ggplot2::aes(x = `max(N)`, y = `ESS(pi0)`,
                                    colour = check_null), size = 2) +
    ggplot2::geom_point(data = dplyr::filter(admissable_null,
                                             check_null == "Not admissable"),
                        ggplot2::aes(x = `max(N)`, y = `ESS(pi0)`),
                        colour = "black") +
    ggrepel::geom_label_repel(data = dplyr::filter(admissable_null,
                                                   check_null %in% c("Minimax",
                                                                     "Admissable",
                                                                     "Null-optimal",
                                                                     "Minimax/Null-optimal")),
                              ggplot2::aes(x = `max(N)`, y = `ESS(pi0)`,
                                           label = designs_null,
                                           fill = check_null),
                              color = "grey", alpha = 0.9,
                              box.padding = ggplot2::unit(0.25, "lines"),
                              show.legend = F) +
    ggplot2::xlab(expression(italic(N[J]))) +
    ggplot2::ylab(expression(paste(italic(ESS), "(",
                                   pi[0], ")",
                                   sep = ""))) +
    ggplot2::ylim(0.97*min(admissable_null$`ESS(pi0)`),
                  1.03*max(admissable_null$`ESS(pi0)`)) +
    ggplot2::scale_x_continuous(minor_breaks = NULL) +
    theme_singlearm()

  plot_des$admissable_pi1 <- ggplot2::ggplot() +
    ggplot2::geom_line(data = admissable_alt, ggplot2::aes(x = `max(N)`,
                                                            y = `ESS(pi1)`)) +
    ggplot2::geom_point(data = dplyr::filter(admissable_alt,
                                             check_alt != "Not admissable"),
                        ggplot2::aes(x = `max(N)`, y = `ESS(pi1)`,
                                     colour = check_alt), size = 2) +
    ggplot2::geom_point(data = dplyr::filter(admissable_alt,
                                             check_null == "Not admissable"),
                        ggplot2::aes(x = `max(N)`, y = `ESS(pi1)`),
                        colour = "black") +
    ggrepel::geom_label_repel(data = dplyr::filter(admissable_alt,
                                                   check_null %in% c("Minimax",
                                                                     "Admissable",
                                                                     "Alt-optimal",
                                                                     "Minimax/Alt-optimal")),
                              ggplot2::aes(x = `max(N)`, y = `ESS(pi1)`,
                                           label = designs_alt,
                                           fill = check_alt),
                              color = "grey", alpha = 0.9,
                              box.padding = ggplot2::unit(0.25, "lines"),
                              show.legend = F) +
    ggplot2::xlab(expression(italic(N[J]))) +
    ggplot2::ylab(expression(paste(italic(ESS), "(",
                                   pi[1], ")",
                                   sep = ""))) +
    ggplot2::ylim(0.97*min(admissable_null$`ESS(pi1)`),
                  1.03*max(admissable_null$`ESS(pi1)`)) +
    ggplot2::scale_x_continuous(minor_breaks = NULL) +
    theme_singlearm()

  ##### Outputting #############################################################

  if (output) {
    output <- list(plot_des = plot_des, des = des)
    return(output)
  }

}
