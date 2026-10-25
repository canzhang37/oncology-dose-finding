# fixed + random scenarios, target 0.25

FIXED <- list(
  S1 = c(0.25, 0.40, 0.50, 0.60, 0.70, 0.80),
  S2 = c(0.10, 0.25, 0.40, 0.50, 0.60, 0.70),
  S3 = c(0.05, 0.12, 0.25, 0.40, 0.50, 0.60),
  S4 = c(0.03, 0.08, 0.15, 0.25, 0.40, 0.50),
  S5 = c(0.02, 0.05, 0.10, 0.15, 0.25, 0.40),
  S6 = c(0.01, 0.03, 0.06, 0.10, 0.15, 0.25),
  S7 = c(0.05, 0.08, 0.10, 0.45, 0.55, 0.65),
  S8 = c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
  S9 = c(0.45, 0.55, 0.65, 0.75, 0.85, 0.90)
)
CORRECT <- c(S1 = 1, S2 = 2, S3 = 3, S4 = 4, S5 = 5, S6 = 6, S7 = 3, S8 = 3, S9 = 0)
OVER <- 0.35

## random curve, correct dose M
rand_scen <- function(K = 6) {
  M <- sample.int(K, 1)
  p <- numeric(K)
  p[M] <- runif(1, 0.20, 0.30)
  if (M > 1) for (d in (M - 1):1) p[d] <- p[d + 1] - if (d == M - 1) runif(1, 0.12, 0.25) else runif(1, 0.03, 0.12)
  if (M < K) for (d in (M + 1):K) p[d] <- p[d - 1] + if (d == M + 1) runif(1, 0.12, 0.25) else runif(1, 0.05, 0.20)
  list(p = pmin(pmax(p, 0.01), 0.95), correct = M)
}
