## 04 - competing risks + multi-state
## recurrence and death aren't independent, KM for "recurrence" w/ death censored overestimates
## -> Aalen-Johansen via survfit() w/ a factor status (no extra pkg needed)
library(survival)
library(dplyr)
library(ggplot2)
library(patchwork)
source("R/00_theme.R")

d <- readRDS("output/colon_wide.rds")

# check: ppl who died w/o recurrence have time_rec == time_death?
with(d, table(same = time_rec == time_death, d$cr_event))
# yes for all 38 "death" rows. good

## 1) CIF for first event ----------------------------------------
aj <- survfit(Surv(t_rec_yr, cr_event) ~ rx, data = d)
# aj$pstate cols: (s0) event-free, recurrence, death   -- check aj$states
aj$states

cif <- data.frame(time = aj$time, aj$pstate, rx = rep(names(aj$strata), aj$strata))
names(cif)[2:4] <- c("free", "recurrence", "death")
cif$rx <- factor(sub("rx=", "", cif$rx), levels = names(pal_rx))

# fine-gray just to have the subdistribution HRs
library(cmprsk)
cr <- with(d, crr(t_rec_yr, as.integer(cr_event) - 1, cov1 = model.matrix(~rx, d)[, -1], failcode = 1))
summary(cr)


## 2) multi-state: event-free -> recurrence -> death, plus event-free -> death ----------
# build counting process data. recurrence and death on the same day happens (n=7), nudge recurrence 0.5 day
ms <- d %>%
  mutate(t_rec2 = ifelse(status_rec == 1 & time_rec >= time_death, time_death - 0.5, time_rec)) %>%
  select(id, rx, status_rec, status_death, t_rec2, time_death)
sum(d$status_rec == 1 & d$time_rec >= d$time_death)   # how many same-day

r1 <- ms %>% transmute(id, rx, tstart = 0,
                       tstop = ifelse(status_rec == 1, t_rec2, time_death),
                       state = case_when(status_rec == 1 ~ "recurrence",
                                         status_death == 1 ~ "death",
                                         TRUE ~ "censor"))
r2 <- ms %>% filter(status_rec == 1) %>%
  transmute(id, rx, tstart = t_rec2, tstop = time_death,
            state = ifelse(status_death == 1, "death", "censor"))
msd <- bind_rows(r1, r2) %>% arrange(id, tstart)
msd$state <- factor(msd$state, levels = c("censor", "recurrence", "death"))
msd$tstart <- msd$tstart / 365.25; msd$tstop <- msd$tstop / 365.25

mfit <- survfit(Surv(tstart, tstop, state) ~ rx, data = msd, id = id)
mfit$transitions   # sanity
# NB: here "death" is absorbing from both states, and pstate gives P(in state) at t

occ <- data.frame(time = mfit$time, mfit$pstate, rx = rep(names(mfit$strata), mfit$strata))
names(occ)[2:4] <- c("Alive, no recurrence", "Alive after recurrence", "Dead")
occ$rx <- factor(sub("rx=", "", occ$rx), levels = names(pal_rx))
occ_l <- tidyr::pivot_longer(occ, 2:4, names_to = "state", values_to = "p")
occ_l$state <- factor(occ_l$state, levels = c("Dead", "Alive after recurrence", "Alive, no recurrence"))
occ_l <- occ_l %>% filter(time <= 8)

st_col <- c("Alive, no recurrence" = "#7FB7BE", "Alive after recurrence" = "#F4A259", "Dead" = "#3D405B")

# labels inside the areas at t=2.5 (orange band is widest there) (stack goes free on bottom, then recurrence, dead on top)
o6 <- occ %>% filter(rx == "Obs") %>% slice_min(abs(time - 2.5), n = 1, with_ties = FALSE)
lab <- data.frame(state = names(st_col),
                  y = c(o6[[2]] / 2, o6[[2]] + o6[[3]] / 2, o6[[2]] + o6[[3]] + o6[[4]] / 2),
                  rx = factor("Obs", levels = names(pal_rx)))

fig4 <- ggplot(occ_l, aes(time, p, fill = state)) +
  geom_area(position = "stack", alpha = 0.95, colour = "white", linewidth = 0.15) +
  geom_text(data = lab, aes(x = 2.5, y = y, label = state),
            colour = "white", size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~rx, nrow = 1) +
  scale_fill_manual(values = st_col, breaks = names(st_col)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 8, 2), expand = c(0, 0)) +
  labs(x = "Years", y = "Probability of being in state", fill = NULL,
       title = "Where patients are over time",
       subtitle = "Multi-state Aalen-Johansen estimate. Lev+5FU mainly shrinks the orange + navy bands (fewer recurrences)") +
  theme_sc() +
  theme(strip.text = element_text(face = "bold", size = 13, hjust = 0), panel.spacing = unit(1.4, "lines"),
        panel.grid = element_blank())

## CIF panel underneath
cif_l <- tidyr::pivot_longer(cif, c(recurrence, death), names_to = "event", values_to = "p")
fig4b <- ggplot(cif_l, aes(time, p, colour = rx, linetype = event)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(values = pal_rx) +
  scale_linetype_manual(values = c(recurrence = 1, death = 3),
                        labels = c(recurrence = "Recurrence", death = "Death w/o recurrence")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Years", y = "Cumulative incidence", colour = NULL, linetype = NULL,
       subtitle = sprintf("Cumulative incidence of first event. Fine-Gray sHR recurrence, Lev+5FU vs Obs: %.2f (%.2f-%.2f)",
                          exp(cr$coef[2]), exp(cr$coef[2] - 1.96 * sqrt(cr$var[2, 2])),
                          exp(cr$coef[2] + 1.96 * sqrt(cr$var[2, 2])))) +
  theme_sc()

ggsave("figures/fig4_state_occupancy.png", fig4 / fig4b + plot_layout(heights = c(1.3, 1)),
       width = 12, height = 9.5, dpi = 200)
