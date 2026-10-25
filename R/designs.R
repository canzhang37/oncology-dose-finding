# dose-escalation designs: 3+3, BOIN, CRM
# u: latent tolerances, DLT if u < p[d]

TARGET <- 0.25
CSIZE <- 3
NCOHORT <- 12

## 3+3; draw(d): DLTs in next cohort
tpt_core <- function(K, draw) {
  n <- y <- integer(K)
  treat <- function(d) {
    y[d] <<- y[d] + draw(d)
    n[d] <<- n[d] + 3L
  }
  # MTD search downward from j
  confirm <- function(j) {
    while (j >= 1) {
      if (n[j] >= 6) return(j)
      treat(j)
      if (y[j] <= 1) return(j)
      j <- j - 1
    }
    0L
  }
  d <- 1L
  repeat {
    treat(d)
    if (y[d] == 1) treat(d)
    if (y[d] >= 2) {
      mtd <- confirm(d - 1L)
      break
    }
    if (d == K) {
      if (n[d] < 6) treat(d)
      mtd <- if (y[d] <= 1) K else confirm(K - 1L)
      break
    }
    d <- d + 1L
  }
  list(mtd = mtd, n = n, y = y)
}

tpt <- function(p, u) {
  i <- 0L
  tpt_core(length(p), function(d) {
    k <- sum(u[i + 1:3] < p[d])
    i <<- i + 3L
    k
  })
}

## exact 3+3: every path, cohort outcomes 0-3
tpt_exact <- function(p) {
  K <- length(p)
  sel <- numeric(K + 1)
  en <- numeric(K)
  walk <- function(ks, pr) {
    j <- 0L
    res <- tryCatch(tpt_core(K, function(d) {
      j <<- j + 1L
      if (j > length(ks)) stop(structure(class = c("more", "condition"), list(message = "", call = NULL, d = d)))
      ks[j]
    }), more = function(e) e)
    if (inherits(res, "more")) {
      q <- dbinom(0:3, 3, p[res$d])
      for (k in 0:3) if (q[k + 1] > 0) walk(c(ks, k), pr * q[k + 1])
    } else {
      sel[res$mtd + 1] <<- sel[res$mtd + 1] + pr
      en <<- en + pr * res$n
    }
  }
  walk(integer(0), 1)
  list(sel = sel, en = en)
}

## BOIN
boin_bounds <- function(target = TARGET) {
  b <- BOIN::get.boundary(target = target, ncohort = NCOHORT, cohortsize = CSIZE)
  c(e = b$lambda_e, d = b$lambda_d)
}

boin <- function(p, u, lam, target = TARGET, cutoff = 0.95) {
  K <- length(p)
  n <- y <- integer(K)
  elim <- logical(K)
  d <- 1L
  i <- 0L
  stopped <- FALSE
  for (k in seq_len(NCOHORT)) {
    y[d] <- y[d] + sum(u[i + 1:CSIZE] < p[d])
    n[d] <- n[d] + CSIZE
    i <- i + CSIZE
    if (n[d] >= 3 && pbeta(target, y[d] + 1, n[d] - y[d] + 1, lower.tail = FALSE) > cutoff) {
      elim[d:K] <- TRUE
      if (d == 1) {
        stopped <- TRUE
        break
      }
    }
    ph <- y[d] / n[d]
    if (elim[d]) {
      d <- d - 1L
    } else if (ph <= lam[["e"]] && d < K && !elim[d + 1]) {
      d <- d + 1L
    } else if (ph >= lam[["d"]] && d > 1) {
      d <- d - 1L
    }
  }
  mtd <- if (stopped) 0L else BOIN::select.mtd(target = target, npts = n, ntox = y)$MTD
  if (mtd == 99) mtd <- 0L
  list(mtd = as.integer(mtd), n = n, y = y, stopped = stopped)
}

## CRM, empiric model on a grid of beta
crm_setup <- function(skel, scale2 = 1.34, lim = 8, m = 2401) {
  b <- seq(-lim, lim, length.out = m)
  lP <- outer(exp(b), log(skel))  # log p, stays finite at the grid ends
  list(skel = skel, b = b, lp = dnorm(b, 0, sqrt(scale2), log = TRUE),
       lP = lP, l1P = log(-expm1(lP)), P1 = exp(lP[, 1]))
}

crm_post <- function(st, n, y, target = TARGET) {
  ll <- st$lp + drop(st$lP %*% y + st$l1P %*% (n - y))
  w <- exp(ll - max(ll))
  w <- w / sum(w)
  bhat <- sum(w * st$b)
  list(bhat = bhat, ptox = st$skel^exp(bhat), p1_over = sum(w[st$P1 > target]))
}

crm <- function(p, u, st, target = TARGET, safety = 0.90) {
  K <- length(p)
  n <- y <- integer(K)
  d <- 1L
  i <- 0L
  stopped <- FALSE
  for (k in seq_len(NCOHORT)) {
    tox <- sum(u[i + 1:CSIZE] < p[d])
    y[d] <- y[d] + tox
    n[d] <- n[d] + CSIZE
    i <- i + CSIZE
    po <- crm_post(st, n, y, target)
    if (n[1] >= 3 && po$p1_over > safety) {
      stopped <- TRUE
      break
    }
    nxt <- which.min(abs(po$ptox - target))
    cap <- if (tox > 0) d else d + 1L  # no skipping, no escalation after a DLT
    d <- min(nxt, cap)
  }
  if (stopped) {
    mtd <- 0L
  } else {
    tried <- which(n > 0)
    mtd <- tried[which.min(abs(po$ptox[tried] - target))]
  }
  list(mtd = as.integer(mtd), n = n, y = y, stopped = stopped)
}

## slow check: dfcrm::crm every cohort
crm_dfcrm <- function(p, u, skel, target = TARGET, safety = 0.90, st = crm_setup(skel)) {
  K <- length(p)
  n <- y <- integer(K)
  lev <- tox_seq <- integer(0)
  d <- 1L
  i <- 0L
  for (k in seq_len(NCOHORT)) {
    t3 <- as.integer(u[i + 1:CSIZE] < p[d])
    lev <- c(lev, rep(d, CSIZE))
    tox_seq <- c(tox_seq, t3)
    y[d] <- y[d] + sum(t3)
    n[d] <- n[d] + CSIZE
    i <- i + CSIZE
    # safety rule: posterior tail
    if (n[1] >= 3 && crm_post(st, n, y, target)$p1_over > safety) return(0L)
    f <- dfcrm::crm(prior = skel, target = target, tox = tox_seq, level = lev,
                    model = "empiric", method = "bayes")
    cap <- if (sum(t3) > 0) d else d + 1L
    d <- min(f$mtd, cap)
  }
  tried <- which(n > 0)
  tried[which.min(abs(f$ptox[tried] - target))]
}
