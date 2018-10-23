theme_singlearm <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_size = base_size,
                      base_family = base_family) + ggplot2::theme(
                      panel.background =
                                     ggplot2::element_rect(fill = "white",
                                                           colour = NA),
                      panel.border = ggplot2::element_rect(fill = NA,
                                                           colour = "grey70",
                                                           size = 0.5),
                      panel.grid.major =
                        ggplot2::element_line(colour = "grey87", size = 0.25),
                      panel.grid.minor =
                        ggplot2::element_line(colour = "grey87", size = 0.125),
                      axis.ticks = ggplot2::element_line(colour = "grey70",
                                                         size = 0.25),
                      legend.key = ggplot2::element_rect(fill = "white",
                                                         colour = NA),
                      strip.background = ggplot2::element_rect(fill = "grey70",
                                                               colour = NA),
                      strip.text =
                        ggplot2::element_text(colour = "white",
                                              size = ggplot2::rel(0.8)),
                      legend.title = ggplot2::element_blank(),
                      legend.position = "bottom",
                      plot.margin = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3), "cm"),
                      complete = TRUE)
}

sub_num <- function(n) {
  codes <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084", "\u2085",
             "\u2086", "\u2087", "\u2088", "\u2089")
  if (n < 10) {
    code <- codes[n + 1]
  } else {
    code <- paste(codes[n%/%10 + 1], codes[n%%10 + 1], sep = "")
  }
  return(code)
}

check_integer_range <- function(value, name, range) {
  if (any(length(value) != 1, value%%1 != 0, !is.finite(value),
          value <= range[1], value >= range[2])){
    stop(name, " must be a single integer in (", range[1], ",", range[2], ")")
  }
}

check_real_range <- function(value, name, range) {
  if (any(length(value) != 1, value < range[1], value > range[2])){
    stop(name, " must be a single number in [", range[1], ",", range[2], "]")
  }
}

check_real_range_strict <- function(value, name, range, len) {
  if (is.numeric(len)) {
    if (any(length(value) != len, value <= range[1], value >= range[2])){
      stop(name, " must be a single number in (", range[1], ",", range[2], ")")
    }
  } else {
    if (any(value <= range[1], value >= range[2])){
      stop(name, " must be a vector containing numbers in (", range[1], ",", range[2], ")")
    }
  }
}

check_pi <- function(pi, len) {
  if (len == "1") {
    if (any(!is.numeric(pi), length(pi) != 1, pi < 0, pi > 1)){
      stop("pi must be a single number in [0,1]")
    }
  } else {
    if (any(!is.numeric(pi), pi < 0, pi > 1)){
      stop("pi must be a numeric vector containing values in [0,1]")
    }
    if (length(pi) != length(unique(pi))) {
      warning("pi contains duplicate values")
    }
  }
}

check_real_pair_range_strict <- function(value1, value2, name1, name2, range){
  if (any(length(value1) != 1, value1 <= range[1], value1 >= range[2])){
    stop(name1, " must be a single number in (", range[1], ",", range[2], ")")
  }
  if (any(length(value2) != 1, value2 <= range[1], value2 >= range[2])){
    stop(name2, " must be a single number in (", range[1], ",", range[2], ")")
  }
  if (value1 >= value2){
    stop(name1, " must be strictly less than ", name2)
  }
}

check_integer_pair_range <- function(value1, value2, name1, name2, range){
  if (any(length(value1) != 1, value1 <= range[1], value1 >= range[2])){
    stop(name1, " must be a single integer in (", range[1], ",", range[2], ")")
  }
  if (any(length(value2) != 1, value2 <= range[1], value2 >= range[2])){
    stop(name2, " must be a single integer in (", range[1], ",", range[2], ")")
  }
  if (value1 > value2){
    stop(name1, " must be less than or equal to ", name2)
  }
}

check_stopping <- function(futility, efficacy){
  if (!is.logical(futility)){
    stop("futility must be a logical variable")
  }
  if (!is.logical(efficacy)){
    stop("efficacy must be a logical variable")
  }
  if (all(!futility, !efficacy)){
    stop("At least one of futility and efficacy must be set to TRUE")
  }
}

check_optimality <- function(optimality, point_prior, beta_prior, futility,
                             efficacy){
  if (!(optimality %in% c("minimax", "null_ess", "alt_ess", "null_med",
                          "alt_med", "prior", "admissable"))){
    stop("optimality must be set to one of \"minimax\", \"null_ess\", \"alt_ess\", \"null_med\", \"alt_med\", \"prior\" or \"admissable\"")
  }
  if (optimality == "prior"){
    if (all(!is.null(point_prior), !is.null(beta_prior))){
      stop("point_prior and beta_prior cannot both be specified when optimality is set to \"prior\"")
    }
    if (all(is.null(point_prior), is.null(beta_prior))){
      stop("One of point_prior and beta_prior must be specified when optimality is set to \"prior\"")
    }
    if (!is.null(point_prior)){
      if (any(length(point_prior) != 1, !is.finite(point_prior),
              point_prior < 0, point_prior > 1)){
        stop("point_prior must be a single number in [0,1]")
      }
    }
    if (!is.null(beta_prior)){
      if (any(length(beta_prior) != 2, any(!is.finite(beta_prior)),
              any(beta_prior <= 0))){
        stop("beta_prior must be a vector of length 2, containing values in (0,Inf)")
      }
    }
  } else {
    if (!is.null(point_prior)){
      warning("point_prior is not missing. This will have no effect as optimality is not set to \"prior\"")
    }
    if (!is.null(beta_prior)){
      warning("beta_prior is not missing. This will have no effect as optimality is not set to \"prior\"")
    }
  }
  if (all(any(!futility, !efficacy), optimality == "delta_minimax")){
    warning("If futility or efficacy are set to FALSE the delta-minimax design is equivalent to the minimax design")
  }
}

check_logical <- function(value, name) {
  if (!is.logical(value)) {
    stop(name, " must be a logical variable")
  }
}

check_ensign <- function(ensign, futility){
  if (!is.logical(ensign)){
    stop("ensign must be a logical variable")
  }
  if (all(!futility, ensign)){
    stop("ensign cannot be set to TRUE and futility to FALSE")
  }
}

check_k <- function(k, des, add_des) {
  num_add_des <- length(add_des)
  if (num_add_des > 0) {
    Js <- c(des$des$J, rep(0, num_add_des))
    for (i in 1:num_add_des) {
      Js[i + 1] <- eval(add_des)$des$J
    }
    if (!all(k %in% 1:max(Js))) {
      stop("k must contain values in [1,max(J)]")
    }
  } else {
    if (!all(k %in% 1:des$des$J)) {
      stop("k must contain values in [1,des$des$J].")
    }
  }
}

check_sa_opchar_gs <- function(opchar, name) {
  if (!("sa_opchar_gs" %in% class(opchar))) {
    stop(name, " must be of class \"sa_opchar_gs\"")
  }
}

check_sa_est_gs <- function(est, name) {
  if (!("sa_est_gs" %in% class(est))) {
    stop(name, " must be of class \"sa_est_gs\"")
  }
}

check_sa_pval_gs <- function(pval, name) {
  if (!("sa_pval_gs" %in% class(pval))) {
    stop(name, " must be of class \"sa_pval_gs\"")
  }
}

check_sa_ci_gs <- function(ci, name) {
  if (!("sa_ci_gs" %in% class(ci))) {
    stop(name, " must be of class \"sa_ci_gs\"")
  }
}

check_sa_opchar_fixed <- function(opchar, name) {
  if (!("sa_opchar_fixed" %in% class(opchar))) {
    stop(name, " must be of class \"sa_opchar_fixed\"")
  }
}

check_sa_est_fixed <- function(est, name) {
  if (!("sa_est_fixed" %in% class(est))) {
    stop(name, " must be of class \"sa_est_fixed\"")
  }
}

check_sa_pval_fixed <- function(pval, name) {
  if (!("sa_pval_fixed" %in% class(pval))) {
    stop(name, " must be of class \"sa_pval_fixed\"")
  }
}

check_sa_ci_fixed <- function(ci, name) {
  if (!("sa_ci_fixed" %in% class(ci))) {
    stop(name, " must be of class \"sa_ci_fixed\"")
  }
}

check_sa_des_fixed <- function(des, name) {
  if (!("sa_des_fixed" %in% class(des))) {
    stop(name, " must be of class \"sa_des_fixed\"")
  }
  if (any(length(des$des$n) != 1, des$des$n%%1 != 0, !is.finite(des$des$n),
          des$des$n < 1)) {
    stop(name, "$des$n must be a single integer in [1,\u221e)")
  }
  if (any(length(des$des$a) != 1, des$des$a%%1 != 0, !is.finite(des$des$a),
          des$des$a < 0, des$des$a >= des$des$n)) {
    stop(name, "$des$a must be a single integer in [0, ", name, "$des$n - 1]")
  }
  if (des$des$r != des$des$a + 1) {
    stop(name, "$des$r must be equal to ", name, "$des$a + 1")
  }
  if (any(length(des$des$pi0) != 1, des$des$pi0 <= 0, des$des$pi0 >= 1)){
    stop(name, "$des$pi0 must be a single number in (0,1)")
  }
  if (any(length(des$des$pi1) != 1, des$des$pi1 <= 0, des$des$pi1 >= 1)){
    stop(name, "$des$pi1 must be a single number in (0,1)")
  }
  if (des$des$pi1 <= des$des$pi0){
    stop(name, "$des$pi1 must be strictly greater than ", name, "$des$pi0")
  }
  if (any(length(des$des$alpha) != 1, des$des$alpha <= 0, des$des$alpha >= 1)){
    stop(name, "$des$alpha must be a single number in (0,1)")
  }
  if (any(length(des$des$beta) != 1, des$des$beta <= 0, des$des$beta >= 1)){
    stop(name, "$des$beta must be a single number in (0,1)")
  }
}

check_sa_des_gs <- function(des, name) {
  if (!("sa_des_gs" %in% class(des))) {
    stop(name, " must be of class \"sa_des_gs\"")
  }
  if (any(length(des$des$J) != 1, des$des$J%%1 != 0, !is.finite(des$des$J),
          des$des$J < 2)) {
    stop(name, "$des$J must belong to \u2115\u2081")
  }
  if (any(length(des$des$n) != des$des$J, des$des$n%%1 != 0,
          !is.finite(des$des$n), des$des$n < 1)) {
    stop(name, "$des$n must be a numeric vector of length ", name, "$des$J containing integers in [1,\u221e)")
  }
  if (any(!is.numeric(des$des$a), length(des$des$a) != des$des$J)) {
    stop(name, "$des$a must be a numeric vector of length ", name, "$des$J")
  }
  if (any(!is.numeric(des$des$r), length(des$des$r) != des$des$J)) {
    stop(name, "$des$r must be a numeric vector of length ", name, "$des$J")
  }
  for (j in 1:des$des$J) {
    if (des$des$a[j] >= des$des$r[j]) {
      stop("Elements of ", name, "$des$a must be strictly less than their corresponding element in ", name, "$des$r")
    }
    if (any(des$des$a[j] > -Inf & des$des$a[j] < 0, des$des$a[j] == Inf)) {
      stop("Elements in ", name, "$des$a must belong to \u2115\u222a{-\u221e}")
    }
    if (is.finite(des$des$a[j])) {
      if (des$des$a[j]%%1 != 0) {
        stop("Elements in ", name, "$des$a must belong to \u2115\u222a{-\u221e}")
      }
    }
    if (des$des$r[j] <= 0) {
      stop("Elements in ", name, "$des$r must belong to \u2115\u222a{\u221e}")
    }
    if (is.finite(des$des$r[j])) {
      if (des$des$r[j]%%1 != 0) {
        stop("Elements in ", name, "$des$r must belong to \u2115\u222a{\u221e}")
      }
    }
  }
  if (any(!is.finite(des$des$a[des$des$J]), !is.finite(des$des$r[des$des$J]))) {
    stop(name, "$des$a[", name, "$des$J]) and ", name, "$des$r[", name, "$des$J]) must be finite")
  }
  if (des$des$r[des$des$J] != des$des$a[des$des$J] + 1) {
    stop(name, "$des$r[", name, "$des$J] must be equal to ", name, "$des$a[", name, "$des$J] + 1")
  }
  if (any(length(des$des$pi0) != 1, des$des$pi0 <= 0, des$des$pi0 >= 1)){
    stop(name, "$des$pi0 must be a single number in (0,1)")
  }
  if (any(length(des$des$pi1) != 1, des$des$pi1 <= 0, des$des$pi1 >= 1)){
    stop(name, "$des$pi1 must be a single number in (0,1)")
  }
  if (any(length(des$des$alpha) != 1, des$des$alpha <= 0, des$des$alpha >= 1)){
    stop(name, "$des$alpha must be a single number in (0,1)")
  }
  if (any(length(des$des$beta) != 1, des$des$beta <= 0, des$des$beta >= 1)){
    stop(name, "$des$beta must be a single number in (0,1)")
  }
}

check_sa_des_adaptive <- function(des, name) {
  if (!("sa_des_adaptive" %in% class(des))) {
    stop(name, " must be of class \"sa_des_adaptive\"")
  }
}

check_sa_des_gehan <- function(des, name) {
  if (!("sa_des_gehan" %in% class(des))) {
    stop(name, " must be of class \"sa_des_gehan\"")
  }
}

check_sa_opchar_adaptive <- function(des, name) {
  if (!("sa_opchar_adaptive" %in% class(des))) {
    stop(name, " must be of class \"sa_opchar_adaptive\"")
  }
}

check_sa_des_admissable <- function(des, name) {
  if (!("sa_des_admissable" %in% class(des))) {
    stop(name, " must be of class \"sa_des_admissable\"")
  }
}

check_sa_opchar_admissable <- function(des, name) {
  if (!("sa_opchar_admissable" %in% class(des))) {
    stop(name, " must be of class \"sa_opchar_admissable\"")
  }
}

check_sa_des_curtailed <- function(des, name) {
  if (!("sa_des_curtailed" %in% class(des))) {
    stop(name, " must be of class \"sa_des_curtailed\"")
  }
}

check_sa_est_curtailed <- function(est, name) {
  if (!("sa_est_curtailed" %in% class(est))) {
    stop(name, " must be of class \"sa_est_curtailed\"")
  }
}

check_sa_opchar_curtailed <- function(des, name) {
  if (!("sa_opchar_curtailed" %in% class(des))) {
    stop(name, " must be of class \"sa_opchar_curtailed\"")
  }
}

check_sa_des_bayesfreq <- function(des, name) {
  if (!("sa_des_bayesfreq" %in% class(des))) {
    stop(name, " must be of class \"sa_des_bayesfreq\"")
  }
}

check_sa_opchar_bayesfreq <- function(des, name) {
  if (!("sa_opchar_bayesfreq" %in% class(des))) {
    stop(name, " must be of class \"sa_opchar_bayesfreq\"")
  }
}

check_belong <- function(value, name, allowed, len) {
  if (len == "any") {
    for (i in 1:length(value)) {
      if (!(value[i] %in% allowed)){
        stop(name, " must contain values only in ", paste(allowed, collapse = ", "))
      }
    }
  } else {
    if (any(length(value) != 1, !(value %in% allowed))) {
      stop(name, " must be set to one of ", paste(allowed, collapse = ", "))
    }
  }
}

row_match <- function(vec, mat) {
  cvec <- do.call("paste", c(vec[, , drop = FALSE], sep = "\r"))
  cmat <- do.call("paste", c(mat[, , drop = FALSE], sep = "\r"))
  return(match(cvec, cmat, nomatch = NA_integer_))
}

check_gs_boundaries <- function(J, n, a, r) {
  if (any(!is.numeric(n), length(n) != J, n%%1 != 0)) {
    stop("n must be a numeric vector of length ", J,
         ", containing integer elements")
  }
  if (any(!is.numeric(a), length(a) != J)) {
    stop("a must be a numeric vector of length ", J)
  }
  if (any(!is.numeric(r), length(r) != J)) {
    stop("r must be a numeric vector of length ", J)
  }
  if (J > 1) {
    for (j in 1:(J - 1)) {
      if (all(is.finite(a[j]), a[j]%%1 != 0)) {
        stop("Finite elements of a must be integers")
      }
      if (all(is.finite(r[j]), r[j]%%1 != 0)) {
        stop("Finite elements of r must be integers")
      }
    }
  }
  if (!is.finite(a[J])) {
    stop("a[J] must be finite")
  }
  if (!is.finite(r[J])) {
    stop("r[J] must be finite")
  }
  if (a[J] != r[J] - 1) {
    stop("a[J] must be one less than r[J]")
  }
}
