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

filtered_mut_names <- unique(mut_filtered$gene)
filtered_acquired_names <- unique(acquired_filtered$gene)

mut_complete <- mut_filtered %>%
  spread(gene, result) %>%
  rename("plate_id" = ref) %>%
  mutate(plate_id = gsub("(.*?)\\_.+", "\\1", plate_id),
         plate_id = sub("^\\d*-?(\\d{4}-.*)", "\\1", plate_id),
         plate_id = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", plate_id)) %>%
  left_join(isolate_data2) %>%
  mutate(dupl = duplicated(id)) %>%
  filter(dupl == FALSE) %>%
  select_at(.vars = vars(id, type, filtered_mut_names))

acquired_complete <- acquired_filtered %>%
  spread(gene, result) %>%
  rename("plate_id" = ref) %>%
  mutate(plate_id = gsub("(.*?)\\_.+", "\\1", plate_id),
         plate_id = sub("^\\d*-?(\\d{4}-.*)", "\\1", plate_id),
         plate_id = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", plate_id)) %>%
  left_join(isolate_data2) %>%
  mutate(dupl = duplicated(id)) %>%
  filter(dupl == FALSE) %>%
  select_at(.vars = vars(id, type, filtered_acquired_names))

write.table(mut_complete, "data/megares_results.txt", sep = "\t", row.names = F)
write.table(acquired_complete, "data/resfinder_results.txt", sep = "\t", row.names = F)
