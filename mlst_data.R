total_df_offline <- read_excel("data/komplett_oversikt_isolater_art2.xlsx") %>%
  rename(species = NORMart,
         CIP = T_CIP,
         NAL = U_NAL) %>%
  mutate(species = case_when(species == "avi" ~ "Wild bird",
                             species == "chi" ~ "Broiler",
                             species == "fox" ~ "Red fox",
                             species == "pig" ~ "Pig"),
         species = factor(species),
         method = factor(method),
         qubit = as.numeric(qubit),
         `260/280` = as.numeric(`260/280`),
         `260/230` = as.numeric(`260/230`),
         antres = as.numeric(antres),
         read_size = factor(as.character(read_size)),
         group = factor(group),
         CIP = factor(CIP,
                      levels = c("0.03",
                                 "0.06",
                                 "0.12",
                                 "0.25",
                                 "0.5",
                                 "1",
                                 "2",
                                 "4",
                                 "8",
                                 "16"),
                      ordered = TRUE),
         NAL = factor(NAL,
                      levels = c("4",
                                 "8",
                                 "16",
                                 "32",
                                 "64",
                                 "128",
                                 "256"),
                      ordered = TRUE))

mlst_data <- read.table("data/mlst/mlst_results.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(plate_id = saksnr)



isolate_data <- total_df_offline %>%
  select(saksnr, plate_id, species, method)


mlst_clean <- mlst_data %>%
  left_join(isolate_data, by = "plate_id") %>%
  left_join(isolate_data, by = c("saksnr.x" = "saksnr")) %>%
  mutate(saksnr.y = if_else(is.na(saksnr.y) == TRUE, saksnr.x, saksnr.y),
         species.x = if_else(is.na(species.x) == TRUE, species.y, species.x),
         method.x = if_else(is.na(method.x) == TRUE, method.y, method.x)) %>%
  select(saksnr.y, plate_id.x, species.x, method.x, ST, adk, fumC, gyrB, icd, mdh, purA, recA) %>%
  rename("saksnr" = saksnr.y,
         "plate_id" = plate_id.x,
         "species" = species.x,
         "method" = method.x)

write.table(test, "data/mlst/mlst_results.txt", sep = "\t", row.names = FALSE)

test <- mlst_results %>%
  mutate(id = paste(saksnr, method, sep = "-")) %>%
  mutate(dupl = duplicated(id)) %>%
  filter(dupl == FALSE)



test <- megares_results %>%
  gather(key, value, -c(saksnr, type)) %>%
  group_by(key, value) %>%
  count()

test$key


