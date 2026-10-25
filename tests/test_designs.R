library(testthat)
library(data.table)
for (f in list.files(file.path("..", "R"), full.names = TRUE)) source(f)

test_that("3+3 on no toxicity and on certain toxicity", {
  a <- tpt(rep(0, 6), rep(0.5, 36))
  expect_equal(a$mtd, 6L)
  expect_equal(a$n, c(3L, 3L, 3L, 3L, 3L, 6L))
  b <- tpt(rep(1, 6), rep(0.5, 36))
  expect_equal(b$mtd, 0L)
  expect_equal(b$n, c(3L, 0L, 0L, 0L, 0L, 0L))
})

test_that("exact 3+3 sums to one and matches a one-dose hand count", {
  e <- tpt_exact(FIXED$S3)
  expect_equal(sum(e$sel), 1)
  p <- 0.3
  q <- dbinom(0:3, 3, p)
  expect_equal(tpt_exact(p)$sel[2], q[1] * (q[1] + q[2]) + q[2] * q[1])
})

test_that("BOIN boundaries match the closed form", {
  ph <- 0.25; p1 <- 0.6 * ph; p2 <- 1.4 * ph
  le <- log((1 - p1) / (1 - ph)) / log(ph * (1 - p1) / (p1 * (1 - ph)))
  ld <- log((1 - ph) / (1 - p2)) / log(p2 * (1 - ph) / (ph * (1 - p2)))
  expect_equal(unname(boin_bounds()), c(le, ld), tolerance = 1e-6)
})

test_that("CRM posterior mean matches dfcrm", {
  sk <- dfcrm::getprior(0.05, 0.25, 3, 6)
  st <- crm_setup(sk)
  lev <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  tox <- c(0, 0, 0, 0, 1, 0, 1, 1, 0)
  n <- tabulate(lev, 6); y <- tabulate(lev[tox == 1], 6)
  f <- dfcrm::crm(prior = sk, target = 0.25, tox = tox, level = lev, model = "empiric", method = "bayes")
  expect_equal(crm_post(st, n, y)$bhat, f$estimate, tolerance = 1e-4)
})

test_that("two-dose rule is symmetric at delta 0 and favours L with a margin", {
  expect_gt(p_select_low(30, 0.4, 0.4, 0), 0.5)
  expect_gt(p_select_low(30, 0.4, 0.4, 0.1), p_select_low(30, 0.4, 0.4, 0))
})
