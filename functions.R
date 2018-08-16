# Functions

## Save plots
save_plots <- function(filename, plot, dpi) {
  ggsave(filename,
         plot,
         device = "tiff",
         units = "cm",
         height = 20,
         width = 25,
         dpi = dpi)
}