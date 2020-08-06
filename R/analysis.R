source("R/packages.R")
source("R/functions.R")
source("R/generic_functions.R")


# Load data ---------------------------------------------------------------
NMR_compounds <- c("alkylC", "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC")
# load("Data/Pyrolysis_Rosinedal_Oct2018_humus_spectrum_prop.RData")
# load("Data/Pyrolysis_Rosinedal_Oct2018_litter_spectrum_prop.RData")
load("Data/Rosinedal_Oct2018_NMR_integrated_standardised.RData")
load("Data/Pyrolysis_Rosinedal_Oct2018_spectrum_area.RData")
load("Data/Rosinedal_2018_IRMS_transect.RData")

# plots to be used (only a fraction of plots from the control was used for NMR)
use_id <- as.character(unique(nmr_intg_raw$id))
transect_dd <- filter(transect_cor, id %in% use_id)
irms_lingon_dd <- irms_lingon_raw %>% 
  filter(id %in% use_id) %>% 
  left_join(transect_dd)
irms_lingon_ful <- irms_lingon_raw %>% 
  left_join(transect_cor)
irms_dd <- irms_res %>% 
  filter(id %in% use_id & layer %in% c("Litter", "Humus")) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")))
pyr_litter_raw <- pyr_litter_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         category = ifelse(category %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                      ifelse(category %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                             as.character(category)))) %>% 
  left_join(irms_dd)
pyr_humus_raw <- pyr_humus_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         category = ifelse(category %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                      ifelse(category %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                             as.character(category)))) %>% 
  left_join(irms_dd)
nmr_intg_raw <- nmr_intg_raw %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")))

# Pyrolysis ---------------------------------------------------------------

# unknown and biproducts
pyr_unk_bip_prop <- bind_rows(pyr_litter_raw, pyr_humus_raw) %>% 
  group_by(category, horizon) %>% 
  summarise(area = sum(area)) %>% 
  group_by(horizon) %>% 
  mutate(prop = (area / sum(area)) * 100,
         prop = round(prop, 2)) %>% 
  filter(category %in% c("CO2", "Unknown")) %>% 
  select(-area) %>% 
  spread(horizon, prop)
knitr::kable(pyr_unk_bip_prop, caption = "Proportion of all product (%)")


# Proportion without biproducts or unknown
pyr_litter_spec <- pyr_litter_raw %>% 
  filter(!(category %in% c("CO2", "Unknown"))) %>% 
  group_by(id) %>% 
  mutate(value = area / sum(area)) %>% 
  ungroup() %>% 
  select(-area)

pyr_humus_spec <- pyr_humus_raw %>% 
  filter(!(category %in% c("CO2", "Unknown"))) %>% 
  group_by(id) %>% 
  mutate(value = area / sum(area)) %>% 
  ungroup() %>% 
  select(-area)


# . Proportion of each compound group at each location --------------------

pyr_comp_ed <- c("Aliphatic deriv.", "Aromatic deriv.", "Carbohydrate", "G-lignin", "N-compound", "Phenol", "S-lignin", "Others")

# All
py_prop <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(category, compound, id, value, layer, location2, horizon) %>% 
  group_by(category, layer, horizon, location2, id) %>% 
  summarise(value = sum(value) * 100) %>% 
  group_by(category, layer, location2, horizon) %>% 
  summarise_at(.vars = vars(value), .funs = funs(M = mean, CI = get_ci)) %>% 
  ungroup() %>% 
  mutate(value = paste0(round(M, 2), "[", round(CI, 2), ']'),
         category = factor(category, levels = pyr_comp_ed)) %>% 
  select(horizon, location2, category, value) %>% 
  spread(category, value)
knitr::kable(filter(py_prop, horizon  == "F/H horizon")[, -1], caption = "Propotion of pyrolysis product for humus")
knitr::kable(filter(py_prop, horizon  == "L horizon")[, -1], caption = "Propotion of pyrolysis product for litter")
write.csv(py_prop, "Output/Tables/Pyrolysis_product.csv", row.names = FALSE)

# Carbohydrate
py_prop_CH <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  filter(category == "Carbohydrate") %>% 
  mutate(compound = ifelse(compound %in% c("Levosugars1", "Levosugars2", "Levosugars3"), "levosugar", compound),
         compound = mapvalues(compound, "Levoglucosenone and Maltol", "Levoglucosenone")) %>% 
  select(compound, value, layer, id) %>% 
  group_by(layer, id, compound) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  group_by(layer, id) %>% 
  mutate(value = value/sum(value)) %>% 
  ungroup() %>% 
  group_by(compound, layer) %>% 
  summarise(value = mean(value)) %>% 
  spread(layer, value) %>% 
  arrange(compound)
knitr::kable(py_prop_CH, caption = "Propotion in carbohydrate")


# N compound
py_Ncomp <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(category, compound, value, layer, location2, id, leaf_d15N, CNratio) %>% 
  group_by(category, layer, location2, id, leaf_d15N, CNratio) %>% 
  summarise(prop = sum(value) * 100) %>% 
  ungroup() %>% 
  filter(category == "N-compound")


d_ply(py_Ncomp, .(layer), function(x){
  print(unique(x$layer))
  print(summary(lm(prop ~ leaf_d15N, x)))
})
ggplot(py_Ncomp, aes(x = leaf_d15N, y = prop))+
  geom_point(aes(col = location2))+
  geom_smooth(method = lm)+
  facet_grid(layer ~ ., scale = "free_y")+
  labs(y = "N compounds (%)")

d_ply(py_Ncomp, .(layer), function(x){
  print(unique(x$layer))
  print(summary(lm(prop ~ CNratio, x)))
})
ggplot(py_Ncomp, aes(x = CNratio, y = prop))+
  geom_point(aes(col = location2))+
  geom_smooth(method = lm)+
  facet_wrap(layer ~ ., nrow = 2, scale = "free")+
  labs(y = "N compounds (%)")

py_Ncomp_ind <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  filter(category == "N-compound") %>% 
  select(category, compound, id, value, layer) %>% 
  group_by(id, layer) %>% 
  mutate(Nprop = value * 100 / sum(value)) %>% 
  group_by(compound, layer) %>% 
  summarise(value = mean(Nprop)) %>% 
  ungroup() %>% 
  select(compound, layer, value) %>% 
  spread(layer, value) %>% 
  arrange(-Humus, -Litter)
knitr::kable(py_Ncomp_ind, caption = "Propotion in N compouds")




# . RDA -------------------------------------------------------------------


# .. Litter ----------------------------------------------------------------

# All
pyr_litter_grp <- unique(pyr_litter_spec$category)
pyr_litter <- pyr_litter_spec %>%
  group_by(category, id, leaf_d15N, location2, CNratio, layer) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(category, value)
pyr_litter_sp   <- decostand(select(pyr_litter, one_of(pyr_litter_grp)), method = "hellinger")
pyr_litter_rda  <- rda(pyr_litter_sp ~ leaf_d15N, pyr_litter)
anova(pyr_litter_rda, nperm = 4999)

# only fertilised plots
pyr_litter_f <- filter(pyr_litter, location2 != "control")
pyr_litter_sp_f   <- decostand(select(pyr_litter_f, one_of(pyr_litter_grp)), method = "hellinger")
pyr_litter_rda_f  <- rda(pyr_litter_sp_f ~ leaf_d15N, pyr_litter_f)
anova(pyr_litter_rda_f, nperm = 4999)


# Species loading
pyr_litter_sp_score <- ldply(list(full = pyr_litter_rda, fertilised = pyr_litter_rda_f),
                            function(x){
                              d <- data.frame(scores(x, display = "species", choices = 1, scaling = 3)) %>% 
                                mutate(compound = row.names(.)) %>% 
                                select(compound, RDA1)
                            }, .id = "dataset") %>% 
  spread(dataset, RDA1) %>% 
  mutate(compound = fct_relevel(compound, "Others", after = 7)) %>% 
  arrange(compound)
write.csv(pyr_litter_sp_score, 
          file = "Output/Tables/RDA_score/RDA_score_pyrolysys_litter.csv",
          row.names = FALSE)



# .. Humus ----------------------------------------------------------------
pyr_humus_grp <- unique(pyr_humus_spec$category)
pyr_humus <- pyr_humus_spec %>%
  group_by(category, id, leaf_d15N, location2, CNratio, layer) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(category, value)
pyr_humus_sp   <- decostand(select(pyr_humus, one_of(pyr_humus_grp)), method = "hellinger")
pyr_humus_rda  <- rda(pyr_humus_sp ~ leaf_d15N, pyr_humus)
anova(pyr_humus_rda, nperm = 4999)

# only fertilised
pyr_humus_f <- filter(pyr_humus, location2 != "control")
pyr_humus_sp_f   <- decostand(select(pyr_humus_f, one_of(pyr_humus_grp)), method = "hellinger")
pyr_humus_rda_f  <- rda(pyr_humus_sp_f ~ leaf_d15N, pyr_humus_f)
anova(pyr_humus_rda_f, nperm = 4999)

# Species loading
pyr_humus_sp_score <- ldply(list(full = pyr_humus_rda, fertilised = pyr_humus_rda_f),
                            function(x){
                              d <- data.frame(scores(x, display = "species", choices = 1, scaling = 3)) %>% 
                                mutate(compound = row.names(.)) %>% 
                                select(compound, RDA1)
                            }, .id = "dataset") %>% 
  spread(dataset, RDA1) %>% 
  mutate(compound = fct_relevel(compound, "Others", after = 7)) %>% 
  arrange(compound)
write.csv(pyr_humus_sp_score, 
          file = "Output/Tables/RDA_score/RDA_score_pyrolysys_humus.csv",
          row.names = FALSE)




# NMR ---------------------------------------------------------------------
nmr_comp_ed <- c("Alkyl C", "Methoxy/N-alkyl C", "O-alkyl C", "Di-O-alkylC", "Aromatic C", "O-aromatic C", "Carbonyl C")

# . Proportion of each compound group -------------------------------------
nmr_intg_smmry <- nmr_intg_raw %>% 
  gather(key = compound, prop, one_of(NMR_compounds)) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         compound = factor(compound, levels = NMR_compounds)) %>% 
  group_by(layer, horizon, location2, treatment, compound) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup()

nmr_prop <- nmr_intg_raw %>% 
  gather(key = compound, prop, one_of(NMR_compounds)) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         compound = factor(compound, levels = NMR_compounds)) %>% 
  group_by(layer, horizon, location2, treatment, compound) %>% 
  summarise_at(.vars = vars(prop), .funs = funs(M = mean, CI = get_ci)) %>% 
  ungroup() %>% 
  mutate(value = paste0(round(M, 2), "[", round(CI, 2), ']'),
         compound2 = mapvalues(compound, 
                              c("alkylC", "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                              nmr_comp_ed),
         compound2 = factor(compound2, levels = nmr_comp_ed)) %>% 
  select(horizon, location2, compound2, value) %>% 
  spread(compound2, value)
knitr::kable(filter(nmr_prop, horizon == "F/H horizon")[, -1], caption = "Propotion of NMR product for the F/H horizon")
knitr::kable(filter(nmr_prop, horizon == "L horizon")[, -1], caption = "Propotion of NMR product for the H horizon")
write.csv(nmr_prop, "Output/Tables/NMR_product.csv", row.names = FALSE)


# . RDA -------------------------------------------------------------------

# .. Litter ----------------------------------------------------------------

# all
nmr_litter <- filter(nmr_intg_raw, layer == "Litter")
nmr_litter_sp  <- decostand(select(nmr_litter, one_of(NMR_compounds)), method = "hellinger")
nmr_litter_rda  <- rda(nmr_litter_sp ~ leaf_d15N, nmr_litter)
anova(nmr_litter_rda, nperm = 4999)

# only fertilised plots
nmr_litter_f <- filter(nmr_intg_raw, layer == "Litter" & treatment == "fertilised")
nmr_litter_sp_f  <- decostand(select(nmr_litter_f, one_of(NMR_compounds)), method = "hellinger")
nmr_litter_rda_f  <- rda(nmr_litter_sp_f ~ leaf_d15N, nmr_litter_f)
anova(nmr_litter_rda_f, nperm = 4999)

# Species loading
nmr_litter_sp_score <- ldply(list(full = nmr_litter_rda, fertilised = nmr_litter_rda_f),
                            function(x){
                              d <- data.frame(scores(x, display = "species", choices = 1, scaling = 3)) %>% 
                                mutate(compound = row.names(.)) %>% 
                                select(compound, RDA1)
                            }, .id = "dataset") %>% 
  spread(dataset, RDA1) %>% 
  arrange(compound)
write.csv(nmr_litter_sp_score, 
          file = "Output/Tables/RDA_score/RDA_score_NMR_litter.csv",
          row.names = FALSE)



# .. Humus ----------------------------------------------------------------

# all
nmr_humus <- filter(nmr_intg_raw, layer == "Humus")
nmr_humus_sp  <- decostand(select(nmr_humus, one_of(NMR_compounds)), method = "hellinger")
nmr_humus_rda  <- rda(nmr_humus_sp ~ leaf_d15N, nmr_humus)
anova(nmr_humus_rda, nperm = 4999)


# only fertilised plots
nmr_humus_f <- filter(nmr_intg_raw, layer == "Humus" & treatment == "fertilised")
nmr_humus_sp_f  <- decostand(select(nmr_humus_f, one_of(NMR_compounds)), method = "hellinger")
nmr_humus_rda_f  <- rda(nmr_humus_sp_f ~ leaf_d15N, nmr_humus_f)
anova(nmr_humus_rda_f, nperm = 4999)


# Species loading
nmr_humus_sp_score <- ldply(list(full = nmr_humus_rda, fertilised = nmr_humus_rda_f),
                             function(x){
                               d <- data.frame(scores(x, display = "species", choices = 1, scaling = 3)) %>% 
                                 mutate(compound = row.names(.)) %>% 
                                 select(compound, RDA1)
                             }, .id = "dataset") %>% 
  spread(dataset, RDA1) %>% 
  arrange(compound)
write.csv(nmr_humus_sp_score, 
          file = "Output/Tables/RDA_score/RDA_score_NMR_humus.csv",
          row.names = FALSE)


# Others ------------------------------------------------------------------

pyr_all_raw <- bind_rows(pyr_litter, pyr_humus) %>% 
  mutate(LCratio = (`G-lignin` + `S-lignin` + Phenol)/Carbohydrate,
         layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon"))) %>% 
  left_join(irms_dd)

# . lignin:carbohydrate vs leaf d15N --------------------------------------
ggplot(pyr_all_raw, aes(x = leaf_d15N, y = log(LCratio), col = layer))+
  geom_point() +
  geom_smooth(method = "lm")

lc_leafd15N_m1 <- lmer(log(LCratio) ~ leaf_d15N * layer + (1|id), pyr_all_raw)
Anova(lc_leafd15N_m1, test.statistic = "F")
qqresidPlot(lc_leafd15N_m1)

# Outlisrs are suggested
lc_leafd15N_m2 <- update(lc_leafd15N_m1, subset=-c(11, 39))
Anova(lc_leafd15N_m2, test.statistic = "F")
qqresidPlot(lc_leafd15N_m2)
# interactions were derived by the outlisers

# try transformation
lc_leafd15N_m3 <- lmer(1/LCratio ~ leaf_d15N * layer + (1|id), pyr_all_raw)
Anova(lc_leafd15N_m3, test.statistic = "F")
qqresidPlot(lc_leafd15N_m3)

# only fertilised plot
lc_leafd15N_fm1 <- lmer(1/LCratio ~ leaf_d15N * layer + (1|id), 
                        filter(pyr_all_raw, treatment == "fertilised"))
Anova(lc_leafd15N_fm1, test.statistic = "F")
qqresidPlot(lc_leafd15N_fm1)

# Summary of lignin:carbohydrate ratios
LC_smmry <- pyr_all_raw %>% 
  group_by(horizon, location2) %>% 
  summarise_at(.vars = vars(LCratio), .funs = funs(M = mean, CI = get_ci)) %>% 
  ungroup() %>% 
  mutate(value = paste0(round(M, 3), "[", round(CI, 3), ']')) %>% 
  select(horizon, location2, value) %>% 
  spread(location2, value)


# . lignin:carbohydrte vs. C mass ----------------------------------------

ggplot(pyr_all_raw, aes(x = log(LCratio), y = Cmass, col = layer))+
  geom_point() +
  geom_smooth(method = "lm")
Cmass_lcr_m1 <- lmer(Cmass ~ log(LCratio) * layer + (1|id), pyr_all_raw)
Anova(Cmass_lcr_m1, test.statistic = "F")
qqresidPlot(Cmass_lcr_m1)
visreg(Cmass_lcr_m1, xvar = "LCratio", by = "layer", overlay = TRUE)

# only fertilised plot
ggplot(filter(pyr_all_raw, treatment == "fertilised"), 
       aes(x = log(LCratio), y = Cmass, color = layer))+
  geom_point()+
  geom_smooth(method = "lm")
Cmass_lcr_fm1 <- lmer(Cmass ~ log(LCratio) + layer + (1|id), 
                      filter(pyr_all_raw, treatment == "fertilised"))
Anova(Cmass_lcr_fm1, test.statistic = "F")
qqresidPlot(Cmass_lcr_fm1)



# Summary of each layer ---------------------------------------------------
irms_lingon_smmry <- irms_lingon_ful %>% 
  group_by(location2) %>% 
  summarise_at(.vars = vars(leaf_d15N), .funs = funs(d15N_M = mean, d15N_CI = get_ci, N = get_n)) %>% 
  mutate(horizon = "leaf")

irms_smmry <- irms_dd %>% 
  group_by(location2, horizon) %>% 
  summarise_at(.vars = vars(Cmass, Nmass, CNratio, d15N), .funs = funs(M = mean, CI = get_ci, N = get_n)) %>% 
  select(-Cmass_N, -Nmass_N, -CNratio_N) %>%
  rename(N = d15N_N) %>% 
  select(horizon, location2, starts_with("Cmass"), starts_with("Nmass"), starts_with("CN"), starts_with("d15N"), N) %>% 
  arrange(horizon) %>% 
  bind_rows(irms_lingon_smmry) 


# Figures ----------------------------------------------------------------
source("R/figs.R")

# Save --------------------------------------------------------------------
save.image("Output/Data/Analysis_results.RData")

