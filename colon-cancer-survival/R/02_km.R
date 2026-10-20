## 02 - KM curves, overall survival + recurrence-free
library(survival)
library(dplyr)
library(ggplot2)
library(patchwork)
source("R/00_theme.R")

d <- readRDS("output/colon_wide.rds")

fit_os <- survfit(Surv(t_death_yr, status_death) ~ rx, data = d)
print(fit_os, rmean = "common")
lr <- survdiff(Surv(t_death_yr, status_death) ~ rx, data = d)
p_lr <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

# tidy the survfit obj by hand (didnt want survminer's default look)
km_df <- function(fit) {
  s <- summary(fit, censored = TRUE)
  strata <- gsub("rx=", "", rep(names(fit$strata), fit$strata))
  out <- data.frame(time = fit$time, surv = fit$surv, lo = fit$lower, hi = fit$upper,
                    n.censor = fit$n.censor, rx = strata)
  # add t=0 row for each arm so the steps start at 1
  start <- data.frame(time = 0, surv = 1, lo = 1, hi = 1, n.censor = 0, rx = unique(strata))
  out <- rbind(start, out)
  out$rx <- factor(out$rx, levels = names(pal_rx))
  out
}

kd <- km_df(fit_os)

# 5-yr survival labels for the right side
s5 <- summary(fit_os, times = 5)
lab5 <- data.frame(rx = gsub("rx=", "", s5$strata), surv = s5$surv)
lab5$rx <- factor(lab5$rx, levels = names(pal_rx))
lab5$txt <- sprintf("%.0f%%", 100 * lab5$surv)

# arm name at the end of each curve
lab_end <- kd %>% group_by(rx) %>% slice_max(time, n = 1, with_ties = FALSE) %>% ungroup()

## stair-step CI. geom_ribbon draws diagonals between event times, so duplicate the points
# first try did index gymnastics, this is simpler: copy each row w/ the next time
stairify <- function(df) {
  df <- df[order(df$time), ]
  n <- nrow(df)
  a <- df; b <- df
  b$time <- c(df$time[-1], NA)
  out <- rbind(a, b)
  out <- out[order(out$time, -out$surv), ]
  out[!is.na(out$time), ]
}
stairs <- do.call(rbind, lapply(split(kd, kd$rx), stairify))

p_km <- ggplot(kd, aes(time, surv, colour = rx)) +
  annotate("rect", xmin = 4.8, xmax = 5.2, ymin = 0, ymax = 1, fill = "grey93") +
  geom_ribbon(data = stairs, aes(x = time, ymin = lo, ymax = hi, fill = rx),
              alpha = 0.13, inherit.aes = FALSE) +
  geom_step(linewidth = 1.1) +
  geom_point(data = subset(kd, n.censor > 0), shape = 124, size = 2.2, show.legend = FALSE) +
  geom_point(data = lab5, aes(x = 5, y = surv), size = 2.6, shape = 21, fill = "white", stroke = 1.2) +
  ggrepel::geom_text_repel(data = lab5, aes(x = 5, y = surv, label = txt), nudge_x = -0.45,
                           direction = "y", size = 3.5, fontface = "bold", segment.colour = NA,
                           show.legend = FALSE, seed = 1) +
  ggrepel::geom_text_repel(data = lab_end, aes(x = time + 0.15, y = surv, label = rx), hjust = 0,
                           direction = "y", nudge_x = 0.2, size = 4, fontface = "bold",
                           segment.colour = "grey70", xlim = c(9.3, 10.6), show.legend = FALSE, seed = 2) +
  annotate("text", x = 5, y = 0.04, label = "5 yr", colour = "grey45", size = 3.2) +
  annotate("label", x = 0.3, y = 0.12, hjust = 0, size = 3.4, label.size = 0, fill = "white",
           label = paste0("Log-rank p = ", format.pval(p_lr, digits = 2, eps = 0.001))) +
  scale_colour_manual(values = pal_rx) +
  scale_fill_manual(values = pal_rx) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = 0:9, limits = c(-0.3, 10.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Overall survival", colour = NULL, fill = NULL,
       title = "Levamisole + 5-FU improves survival after resection of stage C colon cancer",
       subtitle = "Kaplan-Meier estimates with 95% CI, NCCTG/Intergroup trial (n = 929)") +
  theme_sc() + theme(legend.position = "none")

# number at risk table
tt <- 0:9
nr <- summary(fit_os, times = tt, extend = TRUE)
risk <- data.frame(time = nr$time, n = nr$n.risk, rx = factor(gsub("rx=", "", nr$strata), levels = names(pal_rx)))

p_risk <- ggplot(risk, aes(time, rx, label = n, colour = rx)) +
  geom_text(size = 3.4) +
  scale_colour_manual(values = pal_rx, guide = "none") +
  scale_x_continuous(breaks = 0:9, limits = c(-0.3, 10.6), expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(risk$rx))) +
  labs(x = "Years since randomization", y = NULL, subtitle = "Number at risk") +
  theme_sc(11) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(colour = rev(pal_rx), face = "bold"),
        plot.subtitle = element_text(size = 10, face = "bold"))

fig1 <- p_km / p_risk + plot_layout(heights = c(4, 1)) &
  theme(plot.background = element_rect(fill = "#FCFBF8", colour = NA))
ggsave("figures/fig1_km_os.png", fig1, width = 10, height = 7.2, dpi = 200)


## RFS, faceted vs OS -- quick side by side
fit_rfs <- survfit(Surv(t_rec_yr, status_rec) ~ rx, data = d)
both <- rbind(cbind(km_df(fit_os), endpoint = "Overall survival"),
              cbind(km_df(fit_rfs), endpoint = "Recurrence-free"))
fig1b <- ggplot(both, aes(time, surv, colour = rx)) +
  geom_step(linewidth = 1) +
  facet_wrap(~endpoint) +
  scale_colour_manual(values = pal_rx) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Years", y = NULL, colour = NULL,
       title = "Most of the benefit shows up early, as fewer recurrences") +
  theme_sc() + theme(strip.text = element_text(face = "bold", size = 12, hjust = 0))
ggsave("figures/fig1b_os_vs_rfs.png", fig1b, width = 10, height = 5, dpi = 200)
