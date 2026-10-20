# run everything in order. takes ~20s
dir.create("figures", showWarnings = FALSE); dir.create("output", showWarnings = FALSE)
for (f in c("R/01_data_prep.R", "R/02_km.R", "R/03_cox.R", "R/04_multistate.R")) {
  message(">> ", f)
  source(f, echo = FALSE)
}
