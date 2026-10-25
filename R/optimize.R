# two-dose comparison, exact probabilities

## P(pick L): xL/n >= xH/n - delta
p_select_low <- function(n, pH, pL, delta = 0.10) {
  x <- 0:n
  keep <- outer(x, x, function(xl, xh) xl >= xh - delta * n - 1e-9)
  sum(outer(dbinom(x, n, pL), dbinom(x, n, pH)) * keep)
}

oc_grid <- function(ns, pH, diffs, delta) {
  g <- expand.grid(n = ns, diff = diffs, pH = pH, delta = delta)
  g$p_low <- mapply(function(n, d, h, dl) p_select_low(n, h, h - d, dl), g$n, g$diff, g$pH, g$delta)
  g
}

## smallest n, both 80% targets
smallest_n <- function(g, pH, delta, target = 0.80) {
  g <- as.data.frame(g)  # pH: also a column name
  a <- g[g$pH == pH & g$delta == delta & abs(g$diff) < 1e-9, c("n", "p_low")]
  b <- g[g$pH == pH & g$delta == delta & abs(g$diff - 0.20) < 1e-9, c("n", "p_low")]
  m <- merge(a, b, by = "n", suffixes = c("_flat", "_gap"))
  ok <- m$p_low_flat >= target & (1 - m$p_low_gap) >= target
  stays <- rev(cumprod(rev(ok))) == 1  # holds for every larger n
  c(first = if (any(ok)) min(m$n[ok]) else NA, from = if (any(stays)) min(m$n[stays]) else NA)
}
