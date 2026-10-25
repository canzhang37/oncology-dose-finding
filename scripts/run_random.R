# Part A, 1,000 random scenarios
library(data.table)
library(parallel)
for (f in list.files("R", full.names = TRUE)) source(f)

NSCEN <- 1000
NSIM <- 200
set.seed(2027)
scen <- lapply(seq_len(NSCEN), function(r) rand_scen())

one <- function(r) {
  library(data.table)
  for (f in list.files("R", full.names = TRUE)) source(f)
  lam <- boin_bounds()
  st <- crm_setup(dfcrm::getprior(halfwidth = 0.05, target = TARGET, nu = 3, nlevel = 6))
  s <- scen[[r]]
  tr <- sim_scenario(s$p, s$correct, NSIM, seed = 5000 + r, lam = lam, st = st)
  cbind(scen = r, correct = s$correct, summarise_trials(tr, s$p, s$correct))
}

cl <- makeCluster(min(16, detectCores()))
clusterExport(cl, c("scen", "NSIM"))
t0 <- Sys.time()
res <- rbindlist(parLapply(cl, seq_len(NSCEN), one))
stopCluster(cl)
cat("elapsed", format(Sys.time() - t0), "\n")

curves <- rbindlist(lapply(seq_len(NSCEN), function(r) as.list(c(scen = r, correct = scen[[r]]$correct, setNames(scen[[r]]$p, paste0("p", 1:6))))))
fwrite(curves, "results/random_scenarios.csv")
fwrite(res, "results/random_summary.csv")
print(res[, .(correct_sel = mean(correct_sel), over_sel = mean(over_sel), pts_over = mean(pts_over),
              under_sel = mean(under_sel, na.rm = TRUE), N = mean(N)), by = design])
