# P1-P5 verdicts + readout summaries
library(data.table)
for (f in list.files("R", full.names = TRUE)) source(f)

fx <- fread("results/fixed_summary.csv")
rs <- fread("results/random_summary.csv")
sn <- fread("results/partb_smallest_n.csv")
trials <- as.data.table(readRDS("results/trials_fixed.rds"))
s18 <- paste0("S", 1:8)

avg <- fx[scenario %in% s18, .(correct = mean(correct_sel), pts_over = mean(pts_over)), by = design]
under <- fx[scenario %in% paste0("S", 2:6) & design == "3+3", mean(under_sel)]
ravg <- rs[, .(correct = mean(correct_sel), se = sd(correct_sel) / sqrt(.N)), by = design]
g <- function(d, x, tab) tab[design == d][[x]]

pred <- rbind(
  data.table(id = "P1", check = "BOIN - 3+3, S1-S8 correct selection (>= 10)",
             value = g("BOIN", "correct", avg) - g("3+3", "correct", avg)),
  data.table(id = "P1", check = "CRM - 3+3, S1-S8 correct selection (>= 10)",
             value = g("CRM", "correct", avg) - g("3+3", "correct", avg)),
  data.table(id = "P2", check = "|BOIN - CRM|, S1-S8 (< 5)", value = abs(g("BOIN", "correct", avg) - g("CRM", "correct", avg))),
  data.table(id = "P2", check = "|BOIN - CRM|, random (< 5)", value = abs(g("BOIN", "correct", ravg) - g("CRM", "correct", ravg))),
  data.table(id = "P3", check = "3+3 - BOIN, % patients at overdoses S1-S8 (< 0)",
             value = g("3+3", "pts_over", avg) - g("BOIN", "pts_over", avg)),
  data.table(id = "P3", check = "3+3 - CRM, % patients at overdoses S1-S8 (< 0)",
             value = g("3+3", "pts_over", avg) - g("CRM", "pts_over", avg)),
  data.table(id = "P3", check = "3+3 under-dosing, S2-S6 (>= 30)", value = under),
  data.table(id = "P4", check = "3+3 minus best other, random correct selection (< 0)",
             value = g("3+3", "correct", ravg) - max(g("BOIN", "correct", ravg), g("CRM", "correct", ravg))),
  data.table(id = "P5", check = "smallest n per arm, pH 0.40, delta 0.10 (30-50)",
             value = sn[pH == 0.40 & delta == 0.10, first_n])
)
pred[, met := c(value[1] >= 10, value[2] >= 10, value[3] < 5, value[4] < 5,
                value[5] < 0, value[6] < 0, value[7] >= 30, value[8] < 0, value[9] >= 30 & value[9] <= 50)]
fwrite(pred, "results/predictions.csv")
print(pred[, .(id, check, value = round(value, 1), met)])

cat("\nS1-S8 averages\n"); print(avg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
cat("\nrandom averages\n"); print(ravg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

## random: by dose position, paired diffs
bypos <- rs[, .(correct_sel = mean(correct_sel), over_sel = mean(over_sel), pts_over = mean(pts_over), n = .N),
            by = .(correct, design)][order(correct, design)]
fwrite(bypos, "results/random_by_position.csv")
w <- dcast(rs, scen ~ design, value.var = "correct_sel")
pd <- data.table(pair = c("BOIN - 3+3", "CRM - 3+3", "CRM - BOIN"),
                 mean = c(mean(w$BOIN - w$`3+3`), mean(w$CRM - w$`3+3`), mean(w$CRM - w$BOIN)),
                 se = c(sd(w$BOIN - w$`3+3`), sd(w$CRM - w$`3+3`), sd(w$CRM - w$BOIN)) / sqrt(nrow(w)),
                 share_pos = c(mean(w$BOIN > w$`3+3`), mean(w$CRM > w$`3+3`), mean(w$CRM > w$BOIN)))
fwrite(pd, "results/random_paired.csv")
cat("\nrandom paired differences\n"); print(pd[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## post hoc: patients at overdoses as a count
cnt <- trials[scenario %in% s18, .(n_over = mean(n_over)), by = .(scenario, design)][, .(n_over = mean(n_over)), by = design]
fwrite(cnt, "results/posthoc_overdose_count.csv")
cat("\npost hoc, patients at overdoses per trial (count, S1-S8)\n"); print(cnt)
