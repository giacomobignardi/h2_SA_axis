#Author: Giacomo Bignardi
#Description: compare parcel-wise results across G1 versions (Procrustes-aligned vs unrotated order-/sign-matched; Response 2.2)
#Program: Alignment sensitivity------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

# clear workspace
rm(list = ls())

# set not Open Access working directories
wd_noa_output = file.path(dirname(getwd()), "11_noa_output")

# G1 versions (names are used as run labels)
tags = c(main = "", ordermatched = "_ordermatched", signmatched = "_signmatched")

# FUNCTIONS####
# compare each sensitivity run with the main run, per label
compare = function(d) {
  d %>%
    pivot_wider(id_cols = c(Parcel, label), names_from = run, values_from = c(est, sig)) %>%
    group_by(label) %>%
    summarise(
      n_main     = sum(sig_main, na.rm = TRUE),
      n_order    = sum(sig_ordermatched, na.rm = TRUE),
      n_sign     = sum(sig_signmatched, na.rm = TRUE),
      r_order    = cor(est_main, est_ordermatched, use = "complete.obs"),
      r_sign     = cor(est_main, est_signmatched,  use = "complete.obs"),
      both_order = sum(sig_main & sig_ordermatched, na.rm = TRUE),
      same_order = sum(sig_main & sig_ordermatched & sign(est_main) == sign(est_ordermatched), na.rm = TRUE),
      both_sign  = sum(sig_main & sig_signmatched, na.rm = TRUE),
      same_sign  = sum(sig_main & sig_signmatched & sign(est_main) == sign(est_signmatched), na.rm = TRUE),
      .groups = "drop") %>%
    mutate(kept_order   = both_order / n_main,   # share of main parcels still significant
           nested_order = both_order / n_order,  # share of sensitivity parcels also significant in main
           kept_sign    = both_sign / n_main,
           nested_sign  = both_sign / n_sign)
}

# 07: PHENOTYPIC####
# Bonferroni over 400 parcels per label (as in 07)
read_07 = function(tag) {
  read_csv(sprintf("%s/07_SMEM_output%s.csv", wd_noa_output, tag), show_col_types = FALSE) %>%
    filter(label %in% c("rP1P2", "rP1P3", "rP2P3")) %>%
    distinct(Parcel, label, .keep_all = TRUE) %>%
    group_by(label) %>% mutate(sig = p.adjust(pvalue, "bonferroni") < .05) %>% ungroup() %>%
    transmute(Parcel, label, est = est.std, sig)
}

d_07   = imap_dfr(tags, ~ read_07(.x) %>% mutate(run = .y))
cmp_07 = compare(d_07)
cmp_07

# 09: TWIN####
# stored Bonferroni (as in 09); out-of-bounds h2 parcels and inadmissible |r| > 1 dropped per run
files_09 = sprintf("%s/09_CPMEM_output%s.csv", wd_noa_output, tags)

if (all(file.exists(files_09))) {
  load_09 = function(tag) read_csv(sprintf("%s/09_CPMEM_output%s.csv", wd_noa_output, tag), show_col_types = FALSE)
  exc_09  = function(x) x %>% filter(grepl("^h2", label), est > 1) %>% pull(Parcel) %>% unique()
  
  runs_09 = map(tags, load_09)
  
  # excluded parcels (h2 > 1) per run
  print(map(runs_09, exc_09))
  
  # inadmissible correlations per run
  inadm_09 = imap_dfr(runs_09, ~ .x %>%
                        filter(grepl("^r[AE]_", label), !Parcel %in% exc_09(.x)) %>%
                        distinct(Parcel, label, .keep_all = TRUE) %>%
                        group_by(label) %>% summarise(n_inadmissible = sum(!admissible), .groups = "drop") %>%
                        mutate(run = .y)) %>%
    pivot_wider(names_from = run, values_from = n_inadmissible)
  print(inadm_09)
  
  # mean heritability per run
  h2_09 = imap_dfr(runs_09, ~ .x %>%
                     filter(grepl("^h2", label), !Parcel %in% exc_09(.x)) %>%
                     distinct(Parcel, label, .keep_all = TRUE) %>%
                     group_by(label) %>% summarise(mean = mean(est), sd = sd(est), n = n(), .groups = "drop") %>%
                     mutate(run = .y)) %>%
    pivot_wider(names_from = run, values_from = c(mean, sd, n))
  print(h2_09)
  
  # genetic and environmental correlations
  d_09 = imap_dfr(runs_09, ~ .x %>%
                    filter(grepl("^r[AE]_", label), !Parcel %in% exc_09(.x), admissible) %>%
                    distinct(Parcel, label, .keep_all = TRUE) %>%
                    transmute(Parcel, label, est, sig = p_adjust < .05, run = .y))
  cmp_09 = compare(d_09)
  print(cmp_09)
  
  write_csv(h2_09, sprintf("%s/13_alignment_sensitivity_h2.csv", wd_noa_output))
} else {
  message("09 skipped, missing: ", paste(basename(files_09[!file.exists(files_09)]), collapse = ", "))
  cmp_09 = NULL
}

## save####
write_csv(bind_rows(cmp_07 %>% mutate(script = "07"),
                    cmp_09 %>% mutate(script = "09")),
          sprintf("%s/13_alignment_sensitivity.csv", wd_noa_output))