# readout figures
library(data.table)
library(ggplot2)

COL <- c("3+3" = "#2a78d6", "BOIN" = "#eb6834", "CRM" = "#1baf7a")
RAMP <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        legend.position = "top", legend.title = element_blank(),
        axis.text = element_text(colour = "#5f5e5a"), strip.text = element_text(face = "bold", hjust = 0))

## fixed scenarios
fx <- fread("results/fixed_summary.csv")
fx[, design := factor(design, names(COL))]
long <- melt(fx, id.vars = c("scenario", "design"), measure.vars = c("correct_sel", "pts_over"))
long[, variable := factor(variable, c("correct_sel", "pts_over"),
                          c("Trials selecting the correct dose (%)", "Patients treated at overdoses, p > 0.35 (%)"))]
p1 <- ggplot(long, aes(scenario, value, fill = design)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.74) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = COL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = NULL, caption = "S9: every dose is too toxic, so the correct decision is to stop") + th
ggsave("figures/fixed_scenarios.png", p1, width = 8, height = 6, dpi = 150)

## random scenarios by dose position
bp <- fread("results/random_by_position.csv")
bp[, design := factor(design, names(COL))]
lab <- bp[correct == 6]
p2 <- ggplot(bp, aes(correct, correct_sel, colour = design)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  geom_text(data = lab, aes(label = design), hjust = -0.3, size = 3.6, show.legend = FALSE) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = 1:6, limits = c(1, 6.5)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Position of the correct dose (1,000 random curves, 200 trials each)",
       y = "Trials selecting the correct dose (%)") + th
ggsave("figures/random_by_position.png", p2, width = 7, height = 4.2, dpi = 150)

## Part B
g <- fread("results/partb_oc.csv")[pH == 0.40 & delta == 0.10 & diff <= 0.20 + 1e-9]
g[, diff := factor(sprintf("%.2f", diff))]
end <- g[n == 100]
p3 <- ggplot(g, aes(n, p_low, colour = diff)) +
  geom_hline(yintercept = c(0.2, 0.8), colour = "#b4b2a9", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 33, colour = "#b4b2a9", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_text(data = end, aes(label = paste("diff", diff)), hjust = -0.1, size = 3.3, show.legend = FALSE) +
  annotate("text", x = 34, y = 0.04, label = "n = 33 per arm", hjust = 0, size = 3.3, colour = "#5f5e5a") +
  scale_colour_manual(values = RAMP, name = "Response rate at H minus at L") +
  scale_x_continuous(limits = c(10, 115), breaks = seq(10, 100, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Patients per arm", y = "P(select the lower dose)") +
  th + theme(legend.title = element_text(size = 9))
ggsave("figures/dose_optimization.png", p3, width = 7, height = 4.4, dpi = 150)
cat("figures written\n")
