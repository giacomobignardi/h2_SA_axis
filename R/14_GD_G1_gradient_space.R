#Author: Giacomo Bignardi
#Description: parcel-wise GD - G1 associations (phenotypic, measurement-error model) in the space of the first two template gradients (Response R1.5)
#Program: GD - G1 associations in gradient space------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(patchwork)

# clear workspace
rm(list = ls())

# set Open Access working directories
wd_oa = getwd()

# set not Open Access working directories
wd_noa_output = file.path(dirname(getwd()), "11_noa_output")
wd_noa_data   = file.path(dirname(getwd()), "03_rawData/noa")

# network colours (Yeo-Krienen 7, as in other figures)
yeo7_cols = c("Visual"            = "#9F53AA",
              "Somatomotor"       = "#7A9ABD",
              "Dorsal attention"  = "#3d8043",
              "Ventral attention" = "#b584cf",
              "Limbic"            = "#F4FEC8",
              "Frontoparietal"    = "#e8a633",
              "Default"           = "#D8707A")

# LOAD####
# parcel-wise GD - G1 associations (main run of 07), Bonferroni over 400 parcels as in 07
smem = read_csv(sprintf("%s/07_SMEM_output.csv", wd_noa_output), show_col_types = FALSE) %>%
  filter(label == "rP1P3") %>%
  distinct(Parcel, .keep_all = TRUE) %>%
  mutate(sig = p.adjust(pvalue, "bonferroni") < .05) %>%
  select(Parcel, est.std, pvalue, sig)

# template gradients (rows = parcels 1-400; first two gradient columns = G1, G2)
grad = read_csv(sprintf("%s/03_FC_m_G.csv", wd_noa_output), show_col_types = FALSE, name_repair = "unique_quiet") %>%
  select(-any_of("...1")) %>%
  select(1:2) %>%
  set_names(c("G1", "G2")) %>%
  mutate(Parcel = row_number())
stopifnot(nrow(grad) == 400)

# Yeo-Krienen 7 network per parcel
net = read_csv(sprintf("%s/cortical_types.csv", wd_noa_data), show_col_types = FALSE, name_repair = "unique_quiet") %>%
  transmute(Parcel = parcel_Yeo, network = str_split_i(label_Yeo_7, "_", 3)) %>%
  mutate(network = recode(network, Vis = "Visual", SomMot = "Somatomotor", DorsAttn = "Dorsal attention",
                          SalVentAttn = "Ventral attention", Limbic = "Limbic", Cont = "Frontoparietal", Default = "Default"))

d = smem %>%
  left_join(grad, by = "Parcel") %>%
  left_join(net,  by = "Parcel")

# CHECKS####
# n significant parcels (main 07: 230) and gradient sign convention (G1: Default > Somatomotor; G2: Visual > Somatomotor)
print(sum(d$sig))
net_means = d %>% group_by(network) %>% summarise(G1 = mean(G1), G2 = mean(G2), n = n(), .groups = "drop")
print(net_means)
stopifnot(with(net_means, G1[network == "Default"] > G1[network == "Somatomotor"]),
          with(net_means, G2[network == "Visual"]  > G2[network == "Somatomotor"]))

# ASSOCIATIONS WITH GRADIENT POSITION####
# descriptive: p-values do not account for spatial autocorrelation
rho_G1 = cor.test(d$est.std, d$G1, method = "spearman", exact = FALSE)
rho_G2 = cor.test(d$est.std, d$G2, method = "spearman", exact = FALSE)

# robustness: excluding somatomotor parcels
d_nosm = d %>% filter(network != "Somatomotor")
rho_G1_nosm = cor.test(d_nosm$est.std, d_nosm$G1, method = "spearman", exact = FALSE)
rho_G2_nosm = cor.test(d_nosm$est.std, d_nosm$G2, method = "spearman", exact = FALSE)

# joint contribution of G1 and G2
fit = lm(est.std ~ G1 + G2, data = d)

grad_assoc = tibble(
  test = c("rho G1", "rho G2", "rho G1 (no somatomotor)", "rho G2 (no somatomotor)", "R2 G1 + G2"),
  estimate = c(rho_G1$estimate, rho_G2$estimate, rho_G1_nosm$estimate, rho_G2_nosm$estimate, summary(fit)$r.squared),
  n = c(nrow(d), nrow(d), nrow(d_nosm), nrow(d_nosm), nrow(d))) %>%
  mutate(estimate = round(estimate, 2))
print(grad_assoc)

# FIGURE####
# A: networks in gradient space; B: GD - G1 associations, non-significant parcels in grey
p_net = ggplot(d, aes(G1, G2, fill = network)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.2, colour = "grey30") +
  scale_fill_manual(values = yeo7_cols, breaks = names(yeo7_cols)) +
  labs(x = "G1 (template)", y = "G2 (template)", fill = "Network") +
  theme_classic()

p_r = ggplot(d, aes(G1, G2)) +
  geom_point(data = filter(d, !sig), shape = 21, size = 2.5, stroke = 0.2,
             fill = "grey85", colour = "grey60") +
  geom_point(data = filter(d, sig), aes(fill = est.std), shape = 21, size = 2.5, stroke = 0.2, colour = "grey30") +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "r (GD, G1)") +
  labs(x = "G1 (template)", y = "G2 (template)") +
  theme_classic()

p_r = ggplot(d, aes(G1, G2, fill = est.std)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.2, colour = "grey30") +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "r (GD, G1)") +
  labs(x = "G1 (template)", y = "G2 (template)") +
  theme_classic()

SFig_GD_G1 = p_net + p_r + plot_annotation(tag_levels = "A")
SFig_GD_G1

## save####
ggsave(file.path(wd_oa, "Figures", "14_SFig_GD_G1_gradient_space.pdf"), SFig_GD_G1, width = 11, height = 4.5)