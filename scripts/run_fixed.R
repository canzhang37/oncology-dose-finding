# Part A, fixed scenarios S1-S9
library(data.table)
library(parallel)
for (f in list.files("R", full.names = TRUE)) source(f)

NSIM <- 4000
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) NSIM <- as.integer(args[1])

one <- function(s) {
  for (f in list.files("R", full.names = TRUE)) source(f)
  library(data.table)
  lam <- boin_bounds()
  st <- crm_setup(dfcrm::getprior(halfwidth = 0.05, target = TARGET, nu = 3, nlevel = 6))
  nm <- names(FIXED)[s]
  tr <- sim_scenario(FIXED[[s]], CORRECT[[nm]], NSIM, seed = 1000 + s, lam = lam, st = st)
  tr[, scenario := nm]
}

cl <- makeCluster(min(9, detectCores()))
clusterExport(cl, "NSIM")
t0 <- Sys.time()
trials <- rbindlist(parLapply(cl, seq_along(FIXED), one))
stopCluster(cl)
cat("elapsed", format(Sys.time() - t0), "\n")

dir.create("results", showWarnings = FALSE)
saveRDS(trials, "results/trials_fixed.rds", compress = "xz")
summ <- rbindlist(lapply(names(FIXED), function(nm) {
  s <- summarise_trials(trials[scenario == nm], FIXED[[nm]], CORRECT[[nm]])
  cbind(scenario = nm, s)
}))
sel <- rbindlist(lapply(names(FIXED), function(nm) cbind(scenario = nm, sel_table(trials[scenario == nm]))))
fwrite(summ, "results/fixed_summary.csv")
fwrite(sel, "results/fixed_selection.csv")
print(summ[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)], nrows = 100)
