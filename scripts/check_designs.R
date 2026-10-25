# code checks: exact 3+3, get.oc, dfcrm
library(data.table)
for (f in list.files("R", full.names = TRUE)) source(f)
trials <- as.data.table(readRDS("results/trials_fixed.rds"))
nsim <- trials[, uniqueN(trial)]

## 3+3 vs exact enumeration
tp <- rbindlist(lapply(names(FIXED), function(nm) {
  e <- tpt_exact(FIXED[[nm]])
  tr <- trials[scenario == nm & design == "3+3"]
  obs <- tabulate(tr$mtd + 1, 7) / nrow(tr)
  se <- sqrt(e$sel * (1 - e$sel) / nrow(tr))
  data.table(scenario = nm, outcome = c("none", paste0("d", 1:6)), exact = 100 * e$sel, sim = 100 * obs,
             z = ifelse(se > 0, (obs - e$sel) / se, 0),
             exact_N = sum(e$en), sim_N = mean(tr$N), z_N = (mean(tr$N) - sum(e$en)) / (sd(tr$N) / sqrt(nrow(tr))))
}))
fwrite(tp, "results/check_tpt.csv")
cat("3+3 vs exact: max |z| selection", round(max(abs(tp$z)), 2), " max |z| N", round(max(abs(tp$z_N)), 2), "\n")

## BOIN vs BOIN::get.oc
bo <- rbindlist(lapply(seq_along(FIXED), function(s) {
  nm <- names(FIXED)[s]
  o <- BOIN::get.oc(target = TARGET, p.true = FIXED[[nm]], ncohort = NCOHORT, cohortsize = CSIZE,
                    ntrial = nsim, seed = 6 + s)
  tr <- trials[scenario == nm & design == "BOIN"]
  mine <- tabulate(tr$mtd + 1, 7) / nrow(tr)
  ref <- c(o$percentstop, o$selpercent) / 100
  se <- sqrt(mine * (1 - mine) / nrow(tr) + ref * (1 - ref) / nsim)
  data.table(scenario = nm, outcome = c("none", paste0("d", 1:6)), get_oc = 100 * ref, mine = 100 * mine,
             z = ifelse(se > 0, (mine - ref) / se, 0))
}))
fwrite(bo, "results/check_boin.csv")
cat("BOIN vs get.oc: max |z|", round(max(abs(bo$z)), 2), "\n")

## CRM posterior vs dfcrm::crm on random data
sk <- dfcrm::getprior(halfwidth = 0.05, target = TARGET, nu = 3, nlevel = 6)
st <- crm_setup(sk)
set.seed(31)
cr <- rbindlist(lapply(1:200, function(r) {
  lev <- sort(sample(1:6, 3 * sample(1:12, 1), replace = TRUE, prob = c(6, 5, 4, 3, 2, 1)))
  tox <- rbinom(length(lev), 1, FIXED$S4[lev])
  n <- tabulate(lev, 6); y <- tabulate(lev[tox == 1], 6)
  f <- dfcrm::crm(prior = sk, target = TARGET, tox = tox, level = lev, model = "empiric", method = "bayes")
  po <- crm_post(st, n, y)
  data.table(set = r, b_dfcrm = f$estimate, b_mine = po$bhat, mtd_dfcrm = f$mtd,
             mtd_mine = which.min(abs(po$ptox - TARGET)))
}))
fwrite(cr, "results/check_crm_posterior.csv")
cat("CRM posterior: max |diff| in beta", signif(max(abs(cr$b_dfcrm - cr$b_mine)), 3),
    " same recommended dose", sum(cr$mtd_dfcrm == cr$mtd_mine), "of", nrow(cr), "\n")

## whole trials: fast vs dfcrm CRM
set.seed(77)
sc <- c("S3", "S5", "S7", "S8", "S9")
ct <- rbindlist(lapply(rep(sc, each = 100), function(nm) {
  u <- runif(36)
  data.table(scenario = nm, fast = crm(FIXED[[nm]], u, st)$mtd, slow = crm_dfcrm(FIXED[[nm]], u, sk, st = st))
}))
fwrite(ct, "results/check_crm_trials.csv")
cat("CRM trials: same selection in", sum(ct$fast == ct$slow), "of", nrow(ct), "\n")
