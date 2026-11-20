#Author: Giacomo Bignardi
#Description: collect parcel-wise sensitivity analyses across G1 versions (Procrustes-aligned vs unrotated order-/sign-matched) and write SI tables
#Program: Sensitivity analyses (Responses R1.6, R2.2, R3.3)------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

# clear workspace
rm(list = ls())

# set Open Access working directories
wd_oa = getwd()

# set not Open Access working directories
wd_noa_output = file.path(dirname(getwd()), "11_noa_output")
wd_noa_data   = file.path(dirname(getwd()), "03_rawData/noa")

# G1 versions (names = run ids, values = file suffixes)
runs = c(main = "", ordermatched = "_ordermatched", signmatched = "_signmatched")
run_labels = c(main         = "Procrustes-aligned (main)",
               ordermatched = "Unrotated, best-matching component",
               signmatched  = "Unrotated, first component (sign-matched)")

# readable labels for SI
measure_labels = c(rP1P2  = "r microstructure - G1",  rP1P3  = "r GD - G1",  rP2P3  = "r microstructure - GD",
                   h2_P1  = "h2 G1",                  h2_P2  = "h2 microstructure", h2_P3 = "h2 GD",
                   rA_P12 = "rA microstructure - G1", rA_P13 = "rA GD - G1", rA_P23 = "rA microstructure - GD",
                   rE_P12 = "rE microstructure - G1", rE_P13 = "rE GD - G1", rE_P23 = "rE microstructure - GD")

# Yeo-Krienen 7 network per parcel
networks = c("Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default")
annotations = read_csv(sprintf("%s/cortical_types.csv", wd_noa_data), show_col_types = FALSE, name_repair = "unique_quiet") %>%
  transmute(Parcel = parcel_Yeo, network = str_split_i(label_Yeo_7, "_", 3)) %>%
  mutate(network = recode(network, Vis = "Visual", SomMot = "Somatomotor", DorsAttn = "Dorsal attention",
                          SalVentAttn = "Ventral attention", Limbic = "Limbic", Cont = "Frontoparietal", Default = "Default"),
         network = factor(network, levels = networks))

# READ####
# 07: phenotypic (Bonferroni and BH-FDR over 400 parcels per label, as in 07)
read_07 = function(tag) {
  read_csv(sprintf("%s/07_SMEM_output%s.csv", wd_noa_output, tag), show_col_types = FALSE) %>%
    filter(label %in% c("rP1P2", "rP1P3", "rP2P3")) %>%
    distinct(Parcel, label, .keep_all = TRUE) %>%
    group_by(label) %>%
    mutate(p_bonf = p.adjust(pvalue, "bonferroni"),
           p_fdr  = p.adjust(pvalue, "BH")) %>%
    ungroup() %>%
    transmute(Parcel, label, est = est.std, pvalue, p_bonf, p_fdr)
}

# 09: twin (Bonferroni and FDR as stored in 09; out-of-bounds h2 parcels and inadmissible |r| > 1 dropped per run)
read_09 = function(tag) {
  x   = read_csv(sprintf("%s/09_CPMEM_output%s.csv", wd_noa_output, tag), show_col_types = FALSE, name_repair = "unique_quiet")
  exc = x %>% filter(grepl("^h2", label), est > 1) %>% pull(Parcel) %>% unique()
  x %>%
    filter(!Parcel %in% exc, admissible) %>%
    distinct(Parcel, label, .keep_all = TRUE) %>%
    transmute(Parcel, label, est, pvalue, p_bonf = p_adjust, p_fdr)
}

d_07 = imap_dfr(runs, ~ read_07(.x) %>% mutate(run = .y))
d_09 = imap_dfr(runs, ~ read_09(.x) %>% mutate(run = .y))

# SFILE6: counts per analysis, G1 version and correction####
count_sig = function(d) {
  d %>%
    group_by(run, label) %>%
    summarise(n_tested   = n(),
              n_bonf     = sum(p_bonf < .05),
              n_bonf_pos = sum(p_bonf < .05 & est > 0),
              n_bonf_neg = sum(p_bonf < .05 & est < 0),
              n_fdr      = sum(p_fdr < .05),
              n_uncorr   = sum(pvalue < .05),
              mean_est   = mean(est),
              sd_est     = sd(est),
              median_est = median(est),
              .groups = "drop")
}

sens_counts = bind_rows(count_sig(d_07) %>% mutate(analysis = "Phenotypic (measurement-error model)"),
                        count_sig(d_09) %>% mutate(analysis = "Twin (AE measurement-error model)")) %>%
  mutate(label = factor(label, levels = names(measure_labels)),
         run   = factor(run,   levels = names(runs))) %>%
  arrange(analysis, label, run) %>%
  mutate(measure = measure_labels[as.character(label)],
         G1      = run_labels[as.character(run)]) %>%
  select(analysis, measure, G1, n_tested:median_est) %>%
  mutate(across(where(is.double), \(x) round(x, 3)))

# SFILE7: overlap with the main analysis####
compare = function(d) {
  d %>%
    pivot_wider(id_cols = c(Parcel, label), names_from = run, values_from = c(est, p_bonf)) %>%
    mutate(across(starts_with("p_bonf"), \(p) p < .05, .names = "sig_{sub('p_bonf_', '', .col)}")) %>%
    group_by(label) %>%
    summarise(n_main       = sum(sig_main, na.rm = TRUE),
              n_order      = sum(sig_ordermatched, na.rm = TRUE),
              n_sign       = sum(sig_signmatched, na.rm = TRUE),
              r_order      = cor(est_main, est_ordermatched, use = "complete.obs"),
              r_sign       = cor(est_main, est_signmatched,  use = "complete.obs"),
              shared_order = sum(sig_main & sig_ordermatched, na.rm = TRUE),
              same_order   = sum(sig_main & sig_ordermatched & sign(est_main) == sign(est_ordermatched), na.rm = TRUE),
              shared_sign  = sum(sig_main & sig_signmatched, na.rm = TRUE),
              same_sign    = sum(sig_main & sig_signmatched & sign(est_main) == sign(est_signmatched), na.rm = TRUE),
              .groups = "drop")
}

sens_overlap = bind_rows(compare(d_07) %>% mutate(analysis = "Phenotypic (measurement-error model)"),
                         compare(d_09 %>% filter(!grepl("^h2", label))) %>% mutate(analysis = "Twin (AE measurement-error model)")) %>%
  mutate(label = factor(label, levels = names(measure_labels))) %>%
  arrange(analysis, label) %>%
  mutate(measure = measure_labels[as.character(label)]) %>%
  select(analysis, measure, n_main:same_sign) %>%
  mutate(across(where(is.double), \(x) round(x, 3)))

# SFILE8: network tests for GD - G1 (phenotypic), all parcels per network####
tt = function(z) {
  res = t.test(z)
  tibble(t = unname(res$statistic), df = unname(res$parameter), p = res$p.value)
}

sens_networks = d_07 %>%
  filter(label == "rP1P3") %>%
  left_join(annotations, by = "Parcel") %>%
  group_by(run, network) %>%
  summarise(n = n(),
            mean_r = psych::fisherz2r(mean(psych::fisherz(est))),
            tt(psych::fisherz(est)),
            .groups = "drop") %>%
  group_by(run) %>%
  mutate(p_bonf = p.adjust(p, "bonferroni")) %>%
  ungroup() %>%
  mutate(run = factor(run, levels = names(runs))) %>%
  arrange(run, network) %>%
  mutate(G1 = run_labels[as.character(run)]) %>%
  select(G1, network, n, mean_r, t, df, p, p_bonf) %>%
  mutate(across(c(mean_r, t), \(x) round(x, 3)),
         across(c(p, p_bonf), \(x) signif(x, 3)))

# SFILE9: alignment diagnostics per session (from 03_FC_G.py)####
diag_files = c(`Session 1` = "03_GFC_i_d1_diag.csv",
               `Session 2` = "03_GFC_i_d2_diag.csv",
               `Session average` = "03_GFC_i_diag.csv")
diag_paths = set_names(file.path(wd_noa_output, diag_files), names(diag_files))

if (all(file.exists(diag_paths))) {
  mq = function(x) sprintf("%.3f [%.3f, %.3f]", median(x, na.rm = TRUE),
                           quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
  sens_alignment = imap_dfr(diag_paths, ~ read_csv(.x, show_col_types = FALSE, name_repair = "unique_quiet") %>%
                              summarise(n              = sum(!is.na(best_raw_comp)),
                                        pct_best_comp1 = round(100 * mean(best_raw_comp == 1, na.rm = TRUE), 1),
                                        pct_best_comp2 = round(100 * mean(best_raw_comp == 2, na.rm = TRUE), 1),
                                        pct_best_comp3plus = round(100 * mean(best_raw_comp >= 3, na.rm = TRUE), 1),
                                        r_raw_g1       = mq(r_raw_g1),
                                        r_raw_best     = mq(r_raw_best),
                                        r_aligned_g1   = mq(r_aligned_g1),
                                        expl_g1        = mq(expl_g1),
                                        residual_10c   = mq(resid)) %>%
                              mutate(session = .y, .before = 1))
  print(sens_alignment)
} else {
  message("SFILE9 skipped, missing: ", paste(basename(diag_paths[!file.exists(diag_paths)]), collapse = ", "))
  sens_alignment = NULL
}

# print####
print(sens_counts, n = Inf)
print(sens_overlap, n = Inf)
print(sens_networks, n = Inf)

## save####
write_csv(sens_counts,   sprintf("%s/SI/SFILE6.csv", wd_oa))
write_csv(sens_overlap,  sprintf("%s/SI/SFILE7.csv", wd_oa))
write_csv(sens_networks, sprintf("%s/SI/SFILE8.csv", wd_oa))
if (!is.null(sens_alignment)) write_csv(sens_alignment, sprintf("%s/SI/SFILE9.csv", wd_oa))