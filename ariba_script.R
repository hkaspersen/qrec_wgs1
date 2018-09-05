# Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,tidyr,gridExtra,grid,cowplot,forcats,R.utils,purrr,stringr,kableExtra,knitr,IRdisplay,reprex)

# Import functions
source("src/ariba_functions.R")

# Set input data path
megares_results <- "data/ariba/megares_results"
resfinder_results <- "data/ariba/resfinder_results"

# Specify genes of interest
genes <- c("gyr", "par", "sox", "mar", "acr", "tol", "omp")
acquired_genes <- c("qnr", "bla", "tet", "sul", "oqx")

# Import data
mut_data <- get_ariba_data(megares_results)
acquired_data <- get_ariba_data(resfinder_results)

# Clean data
clean_mut_data <- fix_gene_names(mut_data)
clean_acquired_data <- fix_gene_names(acquired_data)

# Wrangle
mut_table <- create_mut_table(clean_mut_data)
mut_filtered <- filter_mut_table(mut_table)

acquired_table <- create_acquired_table(clean_acquired_data)
acquired_filtered <- filter_acquired_table(acquired_table)

mut_complete <- mut_filtered %>%
  spread(gene, result) %>%
  mutate(ref = gsub("(.*?)\\_.*", "\\1", ref)) %>%
  rename("saksnr" = ref)

acquired_complete <- acquired_filtered %>%
  spread(gene, result) %>%
  mutate(ref = gsub("(.*?)\\_.*", "\\1", ref)) %>%
  rename("saksnr" = ref)

write.table(mut_complete, "data/megares_results.txt", sep = "\t", row.names = F)
write.table(acquired_complete, "data/resfinder_results.txt", sep = "\t", row.names = F)
