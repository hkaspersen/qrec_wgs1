### Collapses data frame to unique lines while ignoring NA's
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

# Confidence interval function
get_binCI <- function(x, n) as.numeric(setNames(binom.test(x,n)$conf.int*100,
                                                c("lwr", "upr")))

# Identifies folder names in input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath, pattern = "amr_report")
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
    left_join(df2, by = "ref_name")
  
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

# Function that handles megares data and returns a data frame with information on whether a gene
# is mutated or not. Includes control for QRDR in gyrA.
create_mut_table <- function(df) {
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
    mutate(gyrA_result = result %>% # control for mutation within QRDR for gyrA
             str_extract_all("\\d+") %>% # from reprex package
             map(as.integer) %>% # converts all to integer
             map_lgl(~ any(.x >= 67L & .x <= 106L)), # returns TRUE/FALSE whether value is within range or not
           gyrA_result = if_else(gene != "gyrA", NA, gyrA_result), # filters out results for all other genes
           gyrA_result = as.integer(gyrA_result)) %>% # converts TRUE/FALSE to 1/0
    mutate(result_total = ifelse(result == "", NA, ifelse(result == ".", 0, 1)),
           result_total = as.integer(result_total),
           result_total = if_else(gene == "gyrA", gyrA_result, result_total),
           type = "mut") %>%
    select(-c(gyrA_result, result)) %>%
    rename("result" = result_total)
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
    mutate(result = ifelse(result == "", 0, 1),
           result = as.character(result),
           type = "gene")
  return(df)
}

# Function that returns a filtered dataframe based on the string "genes"
filter_mut_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes)
  return(df)
}

# Function that returns a filtered dataframe based on the string "acquired_genes"
filter_acquired_table <- function(df) {
  report_genes <- unique(df$gene)
  grep_genes <- c()
  for (gene in acquired_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = gsub("_", "", gene))
  return(df)
}