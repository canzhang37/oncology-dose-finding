# paired simulation + summaries

DESIGNS <- c("3+3", "BOIN", "CRM")

## nsim trials, same tolerances for all designs
sim_scenario <- function(p, correct, nsim, seed, lam, st) {
  set.seed(seed)
  U <- matrix(runif(nsim * 36), nsim)
  over <- p > OVER
  rows <- vector("list", nsim * 3)
  r <- 0L
  for (t in seq_len(nsim)) {
    u <- U[t, ]
    fits <- list(tpt(p, u), boin(p, u, lam), crm(p, u, st))
    for (k in 1:3) {
      f <- fits[[k]]
      r <- r + 1L
      rows[[r]] <- list(trial = t, design = DESIGNS[k], mtd = f$mtd, N = sum(f$n),
                        n_over = sum(f$n[over]), n_correct = if (correct > 0) f$n[correct] else NA_integer_,
                        dlt = sum(f$y))
    }
  }
  data.table::rbindlist(rows)
}

summarise_trials <- function(tr, p, correct) {
  K <- length(p)
  over <- which(p > OVER)
  tr[, .(correct_sel = 100 * mean(mtd == correct),
         over_sel = 100 * mean(mtd %in% over),
         under_sel = if (correct > 1) 100 * mean(mtd >= 1 & mtd < correct) else NA_real_,
         stop = 100 * mean(mtd == 0),
         pts_over = 100 * mean(n_over / N),
         pts_correct = if (correct > 0) 100 * mean(n_correct / N) else NA_real_,
         N = mean(N), dlt = mean(dlt)),
     by = design]
}

sel_table <- function(tr, K = 6) {
  s <- tr[, as.list(100 * tabulate(mtd + 1, K + 1) / .N), by = design]
  data.table::setnames(s, c("design", "none", paste0("d", 1:K)))
  s
}
