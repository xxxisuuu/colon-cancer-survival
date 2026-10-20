## 01 - data prep
## colon data from the survival pkg (Moertel et al. 1990 adjuvant chemo trial, stage B/C)
## each patient has 2 rows: etype==1 recurrence, etype==2 death
## -> reshape to one row per pt w/ both endpoints

library(survival)
library(dplyr)
library(tidyr)

raw <- survival::colon
cat("rows:", nrow(raw), " pts:", length(unique(raw$id)), "\n")

# wide format
colon_w <- raw %>%
  mutate(ev = ifelse(etype == 1, "rec", "death")) %>%
  select(id, rx, sex, age, obstruct, perfor, adhere, nodes, differ, extent, surg, node4,
         ev, time, status) %>%
  pivot_wider(names_from = ev, values_from = c(time, status))

# sanity check -- recurrence time should never be after death time
stopifnot(all(colon_w$time_rec <= colon_w$time_death))

colon_w <- colon_w %>%
  mutate(
    rx = factor(rx, levels = c("Obs", "Lev", "Lev+5FU")),
    sex = factor(sex, 0:1, c("Female", "Male")),
    differ = factor(differ, 1:3, c("Well", "Moderate", "Poor")),
    extent = factor(extent, 1:4, c("Submucosa", "Muscle", "Serosa", "Contiguous")),
    surg = factor(surg, 0:1, c("Short", "Long")),   # time from surgery to registration
    node4 = factor(node4, 0:1, c("1-4 nodes", ">4 nodes")),
    agegrp = cut(age, c(0, 50, 65, 100), labels = c("<50", "50-64", "65+"), right = FALSE),
    t_death_yr = time_death / 365.25,
    t_rec_yr   = time_rec / 365.25
  )

# competing risks outcome: first event = recurrence vs death w/o recurrence
# status_rec==1 -> recurrence (time_rec)
# status_rec==0 & status_death==1 -> died first. time_rec==time_death in that case (checked below)
colon_w$cr_event <- with(colon_w, ifelse(status_rec == 1, 1, ifelse(status_death == 1, 2, 0)))
colon_w$cr_event <- factor(colon_w$cr_event, 0:2, c("censor", "recurrence", "death"))
# table(colon_w$cr_event)

# missing: nodes (18), differ (23)  -- small, just complete-case for the models
colSums(is.na(colon_w))

saveRDS(colon_w, "output/colon_wide.rds")
