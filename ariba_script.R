#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

megares_report_loc <- "data/ariba/megares_results"
resfinder_report_loc <- "data/ariba/resfinder_results"
ac_genes <- c("qnr","aac","oqx")
mut_genes <- c("gyr","par","marR","soxS","tolC","acrR")

# adjust parameters for filtering
if (length(ac_genes) == 1) {
  if (grepl("all", ac_genes, ignore.case = TRUE) == TRUE) {
    ac_genes <- "ALL"
    } else {
      ac_genes <- ac_genes
    }
}

if (length(mut_genes) == 1) {  
  if (grepl("all", mut_genes, ignore.case = TRUE) == TRUE) {
    mut_genes <- "ALL"
    } else {
      mut_genes <- mut_genes
    }
}
# -------------------------------------- Libraries

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,tidyr,gridExtra,grid,
               forcats,purrr,stringr,kableExtra,
               knitr,IRdisplay,reprex,svglite)

# -------------------------------------- Functions
# Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

# Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

# Identifies filenames in input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath, pattern = "amr_report.tsv")
  return(files)
}

# Import ariba data from report.tsv from chosen database used in ariba
get_ariba_data <- function(filepath) {
  files <- file_names(filepath)
  
  data_list <- lapply(files,
                      FUN = function(file) {
                        read.delim(
                          paste0(filepath, "/", file),
                          stringsAsFactors = F,
                          header = TRUE,
                          sep = "\t"
                        )
                      })
  
  names(data_list) <- files
  data <- bind_rows(lapply(data_list, function(x) map(x, as.character)), .id = "ref")
  return(data)
}

# Corrects the gene names found in the "cluster" column
fix_gene_names <- function(df) {
  genes <- unique(df$ref_name)
  new_names <- gsub("^(.*?)\\..*", "\\1", genes)
  new_names <- gsub("_", "", new_names, fixed = T)
  new_names <- gsub("-", "", new_names, fixed = T)
  
  gene_names <- c()
  for (i in new_names) {
    p <- paste(tolower(substring(i, 1,3)), substring(i, 4), sep = "", collapse = " ")
    gene_names <- c(gene_names,p)
  }
  df2 <- data.frame(genes,gene_names) %>%
    mutate(genes = as.character(genes)) %>%
    rename(ref_name = genes)
  
  df <- df %>%
    left_join(df2, by = "ref_name") %>%
    mutate(gene_names = as.character(gene_names),
           ref = gsub("(.*?)_amr_report.tsv", "\\1", ref))
  
  return(df)
}

# Function for selecting genes of interest and filtering the 
# columns in the data frame on the resulting vector
select_genes <- function(df, gene_string) {
  names_df <- names(df)
  x <- "ref"
  for (i in names_df) {
    for (j in gene_string) {
      if (str_detect(i, regex(j, ignore_case = T)) == TRUE) {
        x <- c(x, i)
      }
    }
  }
  df <- df %>%
    select(one_of(x))
  return(df)
}

# Allowed flags from the ARIBA report
flag_selection <- c("19","27","147","155","403",
                    "411","915","923","787","795",
                    "531","539","659","667","787","795")  

# Function that returns info on flag selection
check_flags <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(flag_result = flag %in% flag_selection,
           flag_result = as.integer(flag_result),
           ref = gsub("(.*?)\\_.+", "\\1", ref),
           ref = sub("^\\d*-?(\\d{4}-.*)", "\\1", ref),
           ref = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", ref)) %>%
    rename("gene" = gene_names)
  return(df)
}

# Function that handles megares data and returns a data frame with information on whether a gene
# is mutated or not. Includes control for QRDR in gyrA.
create_mut_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    filter(flag %in% flag_selection) %>%
    mutate(id = 1:n()) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, flag)) %>%
    gather(gene, mut, -ref) %>%
    mutate(mut = ifelse(mut == "" | mut == "." | is.na(mut) == TRUE, 0, mut),
           result = ifelse(mut != 0, 1, 0),
           result = as.integer(result),
           type = "mut")
  return(df)
}                    

# Corrects or mutations in QRDR for gyrA, gyrB, parC and parE genes.
# The function returns columns of 1/0 values for whether or not the 
# mutations reported by ARIBA is within the QRDR in the gene:
# gyrA: AA 67 - 106
# gyrB: AA 333 - 481
# ParC: AA 51 - 170
# parE: AA 366 - 523
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142146/pdf/cjvr68pg229.pdf
fix_gyr_par_results <- function(df) {
  df <- df %>%
    mutate(gyrA_result = mut %>% # control for mutation within QRDR for gyrA
             str_extract_all("\\d+") %>% # from reprex package
             map(as.integer) %>% # converts all to integer
             map_lgl(~ any(.x >= 67L & .x <= 106L)), # returns TRUE/FALSE whether value is within range or not
           gyrA_result = if_else(gene != "gyrA", NA, gyrA_result), # filters out results for all other genes
           gyrA_result = as.integer(gyrA_result),
           gyrB_result = mut %>% # control for mutation within QRDR for gyrB
             str_extract_all("\\d+") %>%
             map(as.integer) %>%
             map_lgl(~ any(.x >= 333L & .x <= 481L)),
           gyrB_result = if_else(gene != "gyrB", NA, gyrB_result),
           gyrB_result = as.integer(gyrB_result),
           parC_result = mut %>% # control for mutation within QRDR for parC
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 51L & .x <= 170L)),
           parC_result = if_else(gene != "parC", NA, parC_result),
           parC_result = as.integer(parC_result),
           parE_result = mut %>% # control for mutation within QRDR for parE
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 366L & .x <= 523L)),
           parE_result = if_else(gene != "parE", NA, parE_result),
           parE_result = as.integer(parE_result)) %>%
    mutate(result_gyr_par = case_when(gene == "gyrA" ~ gyrA_result,
                                      gene == "gyrB" ~ gyrB_result,
                                      gene == "parC" ~ parC_result,
                                      gene == "parE" ~ parE_result)) %>%
    mutate(result_total = ifelse(gene %in% c("gyrA","gyrB","parC","parE"), result_gyr_par, result)) %>%
    select(-c(gyrA_result,
              gyrB_result,
              parC_result,
              parE_result,
              result_gyr_par,
              result))
  return(df)
}

# Function that handles resfinder data and returns a data frame with 
# presence/absence for acquired genes                           
create_acquired_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_amr_report.tsv$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total),
           type = "gene") %>%
    select(-result)
  return(df)
}

# Function that returns a filtered dataframe based on the string "genes"
filter_mut_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in mut_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = as.character(gene))
  return(df)
}

# Function that returns a filtered dataframe based on the string "acquired_genes"
filter_acquired_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in ac_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = gsub("_", "", gene),
           gene = as.character(gene))
  return(df)
}

# calculates percentage of present/absent mutations and genes
calc_stats <- function(df) {
  df <- df %>%
    group_by(gene, result_total) %>%
    count() %>%
    ungroup() %>%
    mutate(result_total = if_else(result_total == 1, "Present", "Absent")) %>%
    spread(result_total, n, fill = 0) %>%
    rowwise() %>%
    mutate(Total = Present + Absent,
           Percent = round(Present/Total*100, 1),
           lwr = round(get_binCI(Present, Total)[1], 1),
           upr = round(get_binCI(Present, Total)[2], 1))
  return(df)
}

# calculates how many mutations are present in the genes
# gyrA, gyrB, parC and parE
calc_no_of_mut <- function(df) {
  df1 <- df %>%
    filter(gene %in% c("gyrA","parC","gyrB","parE")) %>%
    mutate(entries = sapply(strsplit(.$mut, ","), FUN = function(x) {length(x)})) %>%
    separate(mut, into = as.character(c(1:max(.$entries)))) %>%
    select(-entries) %>%
    gather(id, mut, -c(ref, gene, result, type))
  
  df2 <- fix_gyr_par_results(df1)
  
  df3 <- df2 %>%
    mutate(test = ifelse(mut != "0" & result_total == 0, 0, 1)) %>%
    filter(test == 1) %>%
    spread(gene, mut) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate_at(.vars = vars(c("gyrA","gyrB","parC","parE")),
              .funs = function(x) ifelse(x == "0", "", x)) %>%
    mutate(mut_gyrA = sapply(strsplit(.$gyrA, ","), FUN = function(x) {length(x)}),
           mut_gyrB = sapply(strsplit(.$gyrB, ","), FUN = function(x) {length(x)}),
           mut_parC = sapply(strsplit(.$parC, ","), FUN = function(x) {length(x)}),
           mut_parE = sapply(strsplit(.$parE, ","), FUN = function(x) {length(x)})) %>%
    select(-c(type, id, result_total, test)) %>%
    mutate_at(.vars = vars(c("gyrA","gyrB","parC","parE")),
              .funs = function(x) ifelse(x == "", "0", x))
  return(df3)
}

# creates a data frame with one row per sample, and 1/0 results 
# for mutations in respective genes in columns
create_mut_report <- function(df) {
  df <- df %>%
    select(-mut) %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type))
  return(df)
}

# creates a data frame with one row per sample, and 1/0 results
# for acquired genes in columns
create_acquired_report <- function(df) {
  df <- df %>%
    group_by(ref) %>%
    mutate(id = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, type))
  return(df)
}

# prepares data frame for number of mutations heatmap
create_heatmap_df <- function(df) {
  qnr <- names(df)
  
  qnr_cols <- c()
  
  for (i in qnr) {
    if (grepl("qnr", i, ignore.case = TRUE) == TRUE) {
      qnr_cols <- c(qnr_cols, i)
    }
  }
  
  qnr_df <- df %>%
    select_at(vars(ref, qnr_cols)) %>%
    gather(gene, result, -ref) %>%
    mutate(type = "Gene")
  
  qnr_mut_quant <- mut_quant %>%
    select(-c("gyrA","gyrB","parC","parE")) %>%
    mutate(type = "Mut") %>%
    gather(gene, result, -c(ref, type)) %>%
    mutate(gene = case_when(gene == "mut_gyrA" ~ "gyrA",
                            gene == "mut_gyrB" ~ "gyrB",
                            gene == "mut_parC" ~ "parC",
                            gene == "mut_parE" ~ "parE",
                            TRUE ~ gene)) %>%
    rbind(., qnr_df)
  return(qnr_mut_quant)
}

# -------------------------------------- Analysis
# Import data
mut_data <- get_ariba_data(megares_report_loc)
acquired_data <- get_ariba_data(resfinder_report_loc)

# Clean data
clean_mut_data <- fix_gene_names(mut_data)
clean_acquired_data <- fix_gene_names(acquired_data)

# Check flags
mut_flags <- check_flags(clean_mut_data)
acquired_flags <- check_flags(clean_acquired_data)

# Wrangle
mut_table <- create_mut_table(clean_mut_data)
mut_table_fixed <- fix_gyr_par_results(mut_table)
mut_filtered <- filter_mut_table(mut_table_fixed)
acquired_table <- create_acquired_table(clean_acquired_data)
acquired_filtered <- filter_acquired_table(acquired_table)

# Stats
if (length(mut_genes) == 1) {
  if (mut_genes == "ALL") {
    mut_stats <- calc_stats(mut_table_fixed)
    mut_report <- create_mut_report(mut_table_fixed)
    }
} else {
  mut_stats <- calc_stats(mut_filtered)
  mut_report <- create_mut_report(mut_filtered)
}

if (length(ac_genes) == 1) {
  if (ac_genes == "ALL") {
    acquired_stats <- calc_stats(acquired_table)
    acquired_report <- create_acquired_report(acquired_table)
  }
} else {
  acquired_stats <- calc_stats(acquired_filtered)
  acquired_report <- create_acquired_report(acquired_filtered)
}

mut_quant <- suppressWarnings(calc_no_of_mut(mut_table))
heatmap_df <- create_heatmap_df(acquired_report)

# ---------------------------------- Plotting

# Total acquired genes
p1 <- ggplot(acquired_stats, aes(gene, Percent))+
  geom_col(color = "black")+
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0.6)+
  labs(title = "Percent occurrence of acquired resistance genes",
       y = "Percent (%) of isolates")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size = 30))

# Mutations in selected intrinsic genes
p2 <- ggplot(mut_stats, aes(gene, Percent))+
  geom_col(color = "black")+
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0.6)+
  labs(title = "Percent occurrence of mutations in intrinsic genes",
       y = "Percent (%) of isolates")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.text = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        plot.title = element_text(size = 30))

# Heatmap for gyr/par genes and qnr genes
p3 <- ggplot(heatmap_df, aes(gene, reorder(ref,as.numeric(result)), fill = type, alpha = factor(result)))+
  geom_tile()+
  geom_vline(xintercept = 4.5,
             alpha = 0.3)+
  scale_fill_manual(values = c("Mut" = "#1f78b4",
                               "Gene" = "#33a02c"),
                    limits = c("Mut",
                               "Gene"),
                    labels = c(" Intrinsic Genes ",
                               " Acquired Genes"))+
  guides(alpha = FALSE)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 40),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 40),
        legend.spacing = unit(1, "cm"),
        legend.position = "top")+
  coord_fixed(0.1)

# ---------------------------------- Tables

# number of mutations table
mut_no <- mut_quant %>%
  group_by(mut_gyrA, mut_gyrB, mut_parC, mut_parE) %>%
  count()

# number of isolates per mutation combination
mut_comb <- mut_quant %>%
  group_by(gyrA, gyrB, parC, parE) %>%
  count()

# write tables to file
write.table(acquired_report, paste0(output_dir, "acquired_report.txt"), sep = "\t", row.names = FALSE)
write.table(mut_report, paste0(output_dir, "mut_report.txt"), sep = "\t", row.names = FALSE)
write.table(mut_flags, paste0(output_dir, "mut_flags.txt"), sep = "\t", row.names = FALSE)
write.table(acquired_flags, paste0(output_dir, "acquired_flags.txt"), sep = "\t", row.names = FALSE)
write.table(mut_quant, paste0(output_dir, "mut_quant.txt"), sep = "\t", row.names = FALSE)
write.table(mut_no, paste0(output_dir, "no_of_mut.txt"), sep = "\t", row.names = FALSE)
write.table(mut_comb, paste0(output_dir, "mut_combinations.txt"), sep = "\t", row.names = FALSE)