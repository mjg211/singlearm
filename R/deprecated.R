check_sm_gs <- function(s, m, des) {
  if (any(!is.numeric(s), !is.numeric(m))) {
    stop("s and m must be numeric vectors")
  }
  if (all(length(s) > 1, length(m) > 1, length(s) != length(m))) {
    stop("If s and m both have length greater than one, then they must have equal length")
  }
  if (all(length(s) > 1, length(m) == 1)) {
    m <- rep(m, length(s))
  } else if (all(length(m) > 1, length(s) == 1)) {
    s <- rep(s, length(m))
  }
  terminal <- terminal_states_gs(des$J, des$a, des$r, des$n)
  for (i in 1:length(s)) {
    if (sum(apply(terminal, 1, function(row, s, m){
      all(as.numeric(row[1:2]) ==
          as.numeric(c(s, m)))},
      s = s[i], m = m[i])) == 0) {
      stop("States (sᵢ,mᵢ) defined by s = (s₁,s₂,…) and m = (m₁,m₂,…) must belong to the set of possible terminal states for design")
    }
  }
  return(list(s = s, m = m))
}
