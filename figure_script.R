## Script for generating figures for the article

gene_names <- names(mut_report)

gene_names <- gene_names[-1]

ggsave(
  "C:/Users/VI1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/percent_mut_per_species.tiff",
  mut_report %>%
    mutate(
      ref = gsub("(.*?)\\_.+", "\\1", ref),
      ref = sub("^\\d*-?(\\d{4}-.*)", "\\1", ref),
      ref = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", ref)
    ) %>%
    arrange(desc(ref)) %>%
    mutate(id = id) %>%
    left_join(isolate_data, by = "id") %>%
    gather(key, value, gene_names) %>%
    group_by(key, value, species) %>%
    count() %>%
    ungroup() %>%
    mutate(value = if_else(value == 1, "Present", "Absent")) %>%
    spread(value, n, fill = 0) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    ) %>%
    ggplot(aes(key, Percent, fill = species)) +
    geom_col(color = "black",
             position = position_dodge(1)) +
    geom_errorbar(
      aes(ymin = lwr, ymax = upr),
      width = 0.6,
      position = position_dodge(1)
    ) +
    scale_x_discrete(limits = c(
      "gyrA", "gyrB", "parC", "parE",
      "marR", "soxS", "tolC"
    )) +
    scale_fill_brewer(palette = "Paired") +
    labs(y = "Percent (%) of isolates") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 12,
        hjust = 1,
        vjust = 0.4
      ),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12)
    ) +
    facet_wrap(~ species),
  device = "tiff",
  units = "cm",
  dpi = 300,
  height = 25,
  width = 30
)



qnr_genes <- c("qnrA1","qnrB19","qnrB6","qnrB60","qnrS1","qnrS2","qnrS4")


ggsave(
  "C:/Users/VI1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/percent_qnr_per_species.tiff",
  acquired_report %>%
    mutate(
      ref = gsub("(.*?)\\_.+", "\\1", ref),
      ref = sub("^\\d*-?(\\d{4}-.*)", "\\1", ref),
      ref = sub("^(\\d{4}-\\d{2}-\\d*)-1", "\\1", ref)
    ) %>%
    arrange(desc(ref)) %>%
    mutate(id = id) %>%
    left_join(isolate_data, by = "id") %>%
    gather(key, value, qnr_genes) %>%
    group_by(key, value, species) %>%
    count() %>%
    ungroup() %>%
    mutate(value = if_else(value == 1, "Present", "Absent")) %>%
    spread(value, n, fill = 0) %>%
    rowwise() %>%
    mutate(
      Total = Present + Absent,
      Percent = round(Present / Total * 100, 1),
      lwr = round(get_binCI(Present, Total)[1], 1),
      upr = round(get_binCI(Present, Total)[2], 1)
    ) %>%
    ggplot(aes(key, Percent, fill = species)) +
    geom_col(color = "black",
             position = position_dodge(0.9)) +
    geom_errorbar(
      aes(ymin = lwr, ymax = upr),
      width = 0.6,
      position = position_dodge(0.9)
    ) +
    scale_fill_brewer(palette = "Paired") +
    labs(y = "Percent (%) of isolates",
         fill = NULL) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.title.x = element_blank(),
      legend.justification = c(0, 1),
      legend.position = c(0.82, 0.97)
    ),
  device = "tiff",
  units = "cm",
  dpi = 300,
  height = 25,
  width = 30
)







