# Part B, exact + simulation check
library(data.table)
for (f in list.files("R", full.names = TRUE)) source(f)

g <- as.data.table(oc_grid(ns = 10:100, pH = c(0.25, 0.40, 0.60), diffs = c(0, 0.05, 0.10, 0.15, 0.20, 0.25),
                           delta = c(0.05, 0.10, 0.15)))
fwrite(g, "results/partb_oc.csv")

sn <- rbindlist(lapply(c(0.25, 0.40, 0.60), function(h) rbindlist(lapply(c(0.05, 0.10, 0.15), function(dl) {
  s <- smallest_n(g, h, dl)
  data.table(pH = h, delta = dl, first_n = s[["first"]], n_from = s[["from"]])
}))))
fwrite(sn, "results/partb_smallest_n.csv")
print(sn)

## 100,000 simulated comparisons at n = 20, 40, 60
set.seed(99)
ck <- rbindlist(lapply(c(20, 40, 60), function(n) rbindlist(lapply(c(0, 0.20), function(d) {
  xh <- rbinom(1e5, n, 0.40); xl <- rbinom(1e5, n, 0.40 - d)
  sim <- mean(xl / n >= xh / n - 0.10 - 1e-9)
  ex <- p_select_low(n, 0.40, 0.40 - d, 0.10)
  data.table(n = n, diff = d, exact = ex, sim = sim, z = (sim - ex) / sqrt(ex * (1 - ex) / 1e5))
}))))
fwrite(ck, "results/partb_check.csv")
print(ck)
print(g[pH == 0.40 & delta == 0.10 & n %in% c(20, 30, 40, 50, 60) & diff %in% c(0, 0.10, 0.20),
        .(n, diff, p_low = round(p_low, 3))])
