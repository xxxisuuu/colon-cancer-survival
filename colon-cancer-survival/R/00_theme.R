# shared plotting stuff. palette picked by hand
library(ggplot2)

pal_rx <- c("Obs" = "#8E9AAF", "Lev" = "#F4A259", "Lev+5FU" = "#BC4B51")

theme_sc <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.title = element_text(face = "bold", size = base + 4, colour = "#1d1d1f"),
      plot.subtitle = element_text(colour = "grey35", margin = margin(b = 10)),
      plot.caption = element_text(colour = "grey50", size = base - 3, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      axis.title = element_text(colour = "grey30"),
      legend.position = "top",
      legend.justification = "left",
      plot.background = element_rect(fill = "#FCFBF8", colour = NA),
      plot.margin = margin(15, 20, 10, 15)
    )
}
