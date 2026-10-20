## 03 - Cox models
##  (a) adjusted model  (b) subgroup forest plot Lev+5FU vs Obs  (c) nodes as a spline
library(survival)
library(dplyr)
library(ggplot2)
library(broom)
source("R/00_theme.R")

d <- readRDS("output/colon_wide.rds")
set.seed(2024)
dcc <- d %>% filter(!is.na(nodes), !is.na(differ))   # complete case, n should be ~888

m0 <- coxph(Surv(t_death_yr, status_death) ~ rx, data = d)
m1 <- coxph(Surv(t_death_yr, status_death) ~ rx + age + sex + obstruct + perfor + adhere +
              node4 + differ + extent + surg, data = dcc)
summary(m1)

# PH check
zph <- cox.zph(m1)
print(zph)
# global test is significant (p<0.001) -- mostly differ, some obstruct.
# rx itself is fine (p~0.26) so the trt HR should be ok. sensitivity: stratify on differ
m1s <- update(m1, . ~ . - differ + strata(differ))
round(exp(cbind(HR = coef(m1s), confint(m1s)))[1:2, ], 3)   # basically same as m1
# png("figures/zph_node4.png"); plot(zph["node4"]); dev.off()

tab <- tidy(m1, exponentiate = TRUE, conf.int = TRUE)
write.csv(tab, "output/cox_adjusted.csv", row.names = FALSE)


### subgroup forest -------------------------------------------------
dd <- d %>% filter(rx != "Lev") %>% mutate(rx = droplevels(rx))

subgrp <- list(
  "Sex"              = "sex",
  "Age"              = "agegrp",
  "Positive nodes"   = "node4",
  "Differentiation"  = "differ",
  "Tumour extent"    = "extent",
  "Obstruction"      = "obstruct",
  "Adherence"        = "adhere"
)
dd$obstruct <- factor(dd$obstruct, 0:1, c("No", "Yes"))
dd$adhere   <- factor(dd$adhere, 0:1, c("No", "Yes"))
# extent: too few in submucosa/contiguous, collapse
dd$extent <- forcats::fct_collapse(dd$extent, "Submucosa/Muscle" = c("Submucosa", "Muscle"),
                                   "Serosa/Contiguous" = c("Serosa", "Contiguous"))

res <- list()
for (g in names(subgrp)) {
  v <- subgrp[[g]]
  # interaction p
  f_int <- as.formula(paste0("Surv(t_death_yr, status_death) ~ rx * ", v))
  f_main <- as.formula(paste0("Surv(t_death_yr, status_death) ~ rx + ", v))
  p_int <- anova(coxph(f_main, data = dd), coxph(f_int, data = dd))[2, "Pr(>|Chi|)"]
  for (lv in levels(dd[[v]])) {
    sub <- dd[!is.na(dd[[v]]) & dd[[v]] == lv, ]
    fit <- coxph(Surv(t_death_yr, status_death) ~ rx, data = sub)
    ci <- exp(confint(fit))
    res[[length(res) + 1]] <- data.frame(group = g, level = lv, n = nrow(sub),
                                         ev = sum(sub$status_death),
                                         hr = exp(coef(fit)), lo = ci[1], hi = ci[2], p_int = p_int)
  }
}
fr <- bind_rows(res)
overall <- coxph(Surv(t_death_yr, status_death) ~ rx, data = dd)
fr <- rbind(data.frame(group = "All patients", level = "All patients", n = nrow(dd), ev = sum(dd$status_death),
                       hr = exp(coef(overall)), lo = exp(confint(overall))[1], hi = exp(confint(overall))[2],
                       p_int = NA), fr)
write.csv(fr, "output/subgroup_hr.csv", row.names = FALSE)

# build y positions w/ group header rows
rows <- list(); y <- 0
for (g in unique(fr$group)) {
  if (g != "All patients") { y <- y + 1; rows[[length(rows) + 1]] <- data.frame(y = y, lab = g, hdr = TRUE) }
  sub <- fr[fr$group == g, ]
  for (i in seq_len(nrow(sub))) {
    y <- y + 1
    rows[[length(rows) + 1]] <- data.frame(y = y, lab = sub$level[i], hdr = FALSE)
  }
}
rows <- bind_rows(rows)
fr$y <- rows$y[!rows$hdr]
hdr <- rows[rows$hdr, ]
hdr$p_int <- fr$p_int[match(hdr$lab, fr$group)]

fr$txt <- sprintf("%.2f (%.2f-%.2f)", fr$hr, fr$lo, fr$hi)
fr$nev <- paste0(fr$ev, "/", fr$n)
fr$sig <- fr$hi < 1

fig2 <- ggplot(fr, aes(y = y)) +
  # zebra
  geom_rect(data = hdr, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#EFEBE4", inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = fr$hr[1], colour = pal_rx["Lev+5FU"], alpha = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.8, colour = "grey40") +
  geom_point(aes(x = hr, size = n, fill = sig), shape = 22, colour = "grey20") +
  geom_text(aes(x = 0.06, label = level), hjust = 0, size = 3.4, colour = "grey15") +
  geom_text(data = hdr, aes(x = 0.06, y = y, label = lab), hjust = 0, size = 3.6, fontface = "bold",
            inherit.aes = FALSE) +
  geom_text(data = hdr, aes(x = 5.5, y = y, label = paste0("p-int = ", sprintf("%.2f", p_int))),
            size = 3.1, colour = "grey35", inherit.aes = FALSE, hjust = 1) +
  geom_text(aes(x = 2.65, label = nev), size = 3.2, colour = "grey35") +
  geom_text(aes(x = 3.6, label = txt), size = 3.2, hjust = 0, colour = "grey25") +
  annotate("text", x = c(0.06, 2.65, 3.6), y = -0.4, hjust = c(0, 0.5, 0), fontface = "bold", size = 3.3,
           label = c("Subgroup", "Deaths/N", "HR (95% CI)")) +
  annotate("text", x = c(0.55, 1.8), y = max(fr$y) + 1.4, label = c("<- favours Lev+5FU", "favours observation ->"),
           size = 3.2, colour = c(pal_rx["Lev+5FU"], "grey40"), fontface = "italic") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2), limits = c(0.05, 6)) +
  scale_y_reverse(expand = expansion(add = c(1, 1.5))) +
  scale_fill_manual(values = c(`TRUE` = pal_rx[["Lev+5FU"]], `FALSE` = "white"), guide = "none") +
  scale_size_area(max_size = 5, guide = "none") +
  labs(x = "Hazard ratio, death (log scale)", y = NULL,
       title = "Lev+5FU benefit holds across most subgroups",
       subtitle = "Lev+5FU vs observation. Filled = upper CI below 1, size ~ n. Sex interaction (p = 0.04) is probably chance w/ 7 tests") +
  theme_sc() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
fig2
ggsave("figures/fig2_forest.png", fig2, width = 10, height = 9, dpi = 200)


### nodes spline -------------------------------------------------
m_sp <- coxph(Surv(t_death_yr, status_death) ~ rx + pspline(nodes, df = 3) + age + sex, data = dcc)
# termplot gives the partial effect, grab it instead of plotting
tp <- termplot(m_sp, term = 2, se = TRUE, plot = FALSE)$nodes
ref <- tp$y[tp$x == 1]   # reference = 1 node
tp$hr <- exp(tp$y - ref); tp$lo <- exp(tp$y - 1.96 * tp$se - ref); tp$hi <- exp(tp$y + 1.96 * tp$se - ref)
tp <- tp[tp$x <= 20, ]   # very sparse after ~20

fig3 <- ggplot(tp, aes(x, hr)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#BC4B51", alpha = 0.15) +
  geom_line(linewidth = 1.3, colour = "#BC4B51") +
  geom_hline(yintercept = 1, linetype = 3) +
  # jitter by hand, position_jitter wants a y for rugs
  geom_rug(data = transform(dcc[dcc$nodes <= 20, ], nj = nodes + runif(sum(dcc$nodes <= 20), -0.3, 0.3)),
           aes(x = nj), inherit.aes = FALSE, alpha = 0.08, length = unit(0.04, "npc")) +
  annotate("segment", x = 4.5, xend = 4.5, y = 0.8, yend = 4.6, colour = "grey60", linetype = 2) +
  annotate("text", x = 4.7, y = 4.4, hjust = 0, size = 3.3, colour = "grey35",
           label = "usual >4 nodes cut-off") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  labs(x = "Number of positive lymph nodes", y = "HR vs 1 node (log scale)",
       title = "Risk rises steeply up to ~10 nodes, then flattens",
       subtitle = "Penalized spline (df = 3) in a Cox model, adjusted for treatment, age, sex") +
  theme_sc()
ggsave("figures/fig3_nodes_spline.png", fig3, width = 9, height = 5.5, dpi = 200)
