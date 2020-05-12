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
  mutate(layer = factor(layer, levels = c("Litter", "Humus")))
pyr_litter_raw <- pyr_litter_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         grp = ifelse(grp %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                      ifelse(grp %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                             as.character(grp)))) %>% 
  left_join(irms_dd)
pyr_humus_raw <- pyr_humus_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         grp = ifelse(grp %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                      ifelse(grp %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                             as.character(grp)))) %>% 
  left_join(irms_dd)
nmr_intg_raw <- nmr_intg_raw %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")))

# Pyrolysis ---------------------------------------------------------------

# unknown and biproducts
pyr_unk_bip_prop <- bind_rows(pyr_litter_raw, pyr_humus_raw) %>% 
  group_by(grp, layer) %>% 
  summarise(area = sum(area)) %>% 
  group_by(layer) %>% 
  mutate(prop = (area / sum(area)) * 100,
         prop = round(prop, 2)) %>% 
  filter(grp %in% c("co2", "Acetic acid", "unknown")) %>% 
  select(-area) %>% 
  spread(layer, prop)
knitr::kable(pyr_unk_bip_prop, caption = "Proportion of all product (%)")


# Proportion without biproducts or unknown
pyr_litter_spec <- pyr_litter_raw %>% 
  filter(!(grp %in% c("co2", "Acetic acid", "unknown"))) %>% 
  group_by(id) %>% 
  mutate(value = area / sum(area)) %>% 
  ungroup() %>% 
  select(-area)

pyr_humus_spec <- pyr_humus_raw %>% 
  filter(!(grp %in% c("co2", "Acetic acid", "unknown"))) %>% 
  group_by(id) %>% 
  mutate(value = area / sum(area)) %>% 
  ungroup() %>% 
  select(-area)


# . Proportion of each compound group at each location --------------------

pyr_comp_ed <- c("Aliphatic deriv.", "Aromatic", "Carbohydrate", "G-lignin", "N comp.", "Phenol", "S-lignin", "Others")

# All
py_prop <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(grp, comp, id, value, layer, location2) %>% 
  group_by(grp, layer, location2, id) %>% 
  summarise(value = sum(value) * 100) %>% 
  group_by(grp, layer, location2) %>% 
  summarise_at(.vars = vars(value), .funs = funs(M = mean, CI = get_ci)) %>% 
  ungroup() %>% 
  mutate(value = paste0(round(M, 2), "[", round(CI, 2), ']'),
         grp2 = mapvalues(grp, 
                          c("alphatic_derivative", "aromatic", "carbohydrate", 
                            "g_lignin", "N_comp", "Phenol", "s_lignin", "Others"),
                          pyr_comp_ed),
         grp2 = factor(grp2, levels = pyr_comp_ed)) %>% 
  select(layer, location2, grp2, value) %>% 
  spread(grp2, value)
knitr::kable(filter(py_prop, layer == "Humus")[, -1], caption = "Propotion of pyrolysis product for humus")
knitr::kable(filter(py_prop, layer == "Litter")[, -1], caption = "Propotion of pyrolysis product for litter")
write.csv(py_prop, "Output/Tables/Pyrolysis_product.csv", row.names = FALSE)

# Carbohydrate
py_prop_CH <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  filter(grp == "carbohydrate") %>% 
  mutate(comp = ifelse(comp %in% c("beta.-D-Glucopyranose, 1,6-anhydro-", "Levosugars1", "Levosugars2", "Levosugars3"), "levosugar", comp)) %>% 
  select(comp, value, layer, id) %>% 
  group_by(layer, id, comp) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  group_by(layer, id) %>% 
  mutate(value = value/sum(value)) %>% 
  ungroup() %>% 
  group_by(comp, layer) %>% 
  summarise(value = mean(value)) %>% 
  spread(layer, value) %>% 
  arrange(-Humus, -Litter)


knitr::kable(py_prop_CH, caption = "Propotion in carbohydrate")


# N compound
py_Ncomp <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(grp, comp, value, layer, location2, id, leaf_d15N, CNratio) %>% 
  group_by(grp, layer, location2, id, leaf_d15N, CNratio) %>% 
  summarise(prop = sum(value) * 100) %>% 
  ungroup() %>% 
  filter(grp == "N_comp")


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
  filter(grp == "N_comp") %>% 
  select(grp, comp, id, value, layer) %>% 
  group_by(id, layer) %>% 
  mutate(Nprop = value * 100 / sum(value)) %>% 
  group_by(comp, layer) %>% 
  summarise(value = mean(Nprop)) %>% 
  ungroup() %>% 
  select(comp, layer, value) %>% 
  spread(layer, value) %>% 
  arrange(-Humus, -Litter)
knitr::kable(py_Ncomp_ind, caption = "Propotion in N compouds")




# . RDA -------------------------------------------------------------------

# Litter
pyr_litter_grp <- unique(pyr_litter_spec$grp)
pyr_litter <- pyr_litter_spec %>%
  group_by(grp, id, leaf_d15N, location2, CNratio, layer) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(grp, value)
pyr_litter_sp   <- decostand(select(pyr_litter, one_of(pyr_litter_grp)), method = "hellinger")
pyr_litter_rda  <- rda(pyr_litter_sp ~ leaf_d15N, pyr_litter)
anova(pyr_litter_rda, nperm = 4999)


# Humus
pyr_humus_grp <- unique(pyr_humus_spec$grp)
pyr_humus <- pyr_humus_spec %>%
  group_by(grp, id, leaf_d15N, location2, CNratio, layer) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(grp, value)
pyr_humus_sp   <- decostand(select(pyr_humus, one_of(pyr_humus_grp)), method = "hellinger")
pyr_humus_rda  <- rda(pyr_humus_sp ~ leaf_d15N, pyr_humus)
anova(pyr_humus_rda, nperm = 4999)



# NMR ---------------------------------------------------------------------
nmr_comp_ed <- c("Alkyl C", "Methoxy/N-alkyl C", "O-alkyl C", "Di-O-alkylC", "Aromatic C", "O-aromatic C", "Carbonyl C")

# . Proportion of each compound group -------------------------------------
nmr_intg_smmry <- nmr_intg_raw %>% 
  gather(key = compound, prop, one_of(NMR_compounds)) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         compound = factor(compound, levels = NMR_compounds)) %>% 
  group_by(layer, location2, treatment, compound) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup()

nmr_prop <- nmr_intg_raw %>% 
  gather(key = compound, prop, one_of(NMR_compounds)) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         compound = factor(compound, levels = NMR_compounds)) %>% 
  group_by(layer, location2, treatment, compound) %>% 
  summarise_at(.vars = vars(prop), .funs = funs(M = mean, CI = get_ci)) %>% 
  ungroup() %>% 
  mutate(value = paste0(round(M, 2), "[", round(CI, 2), ']'),
         compound2 = mapvalues(compound, 
                              c("alkylC", "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                              nmr_comp_ed),
         compound2 = factor(compound2, levels = nmr_comp_ed)) %>% 
  select(layer, location2, compound2, value) %>% 
  spread(compound2, value)
knitr::kable(filter(nmr_prop, layer == "Humus")[, -1], caption = "Propotion of NMR product for humus")
knitr::kable(filter(nmr_prop, layer == "Litter")[, -1], caption = "Propotion of NMR product for litter")
write.csv(py_prop, "Output/Tables/NMR_product.csv", row.names = FALSE)


# . RDA -------------------------------------------------------------------

# Litter
nmr_litter <- filter(nmr_intg_raw, layer == "Litter")
nmr_litter_sp  <- decostand(select(nmr_litter, one_of(NMR_compounds)), method = "hellinger")
nmr_litter_rda  <- rda(nmr_litter_sp ~ leaf_d15N, nmr_litter)
anova(nmr_litter_rda, nperm = 4999)


# Humus
nmr_humus <- filter(nmr_intg_raw, layer == "Humus")
nmr_humus_sp  <- decostand(select(nmr_humus, one_of(NMR_compounds)), method = "hellinger")
nmr_humus_rda  <- rda(nmr_humus_sp ~ leaf_d15N, nmr_humus)
anova(nmr_humus_rda, nperm = 4999)


# Others ------------------------------------------------------------------


pyr_all_raw <- bind_rows(pyr_litter, pyr_humus) %>% 
  mutate(CLratio = carbohydrate/(g_lignin + s_lignin + Phenol),
         layer = factor(layer, levels = c("Litter", "Humus"))) %>% 
  left_join(irms_dd)

# . carbohydrate:lignin vs leaf d15N --------------------------------------
ggplot(pyr_all_raw, aes(x = leaf_d15N, y = CLratio, col = layer))+
  geom_point() +
  geom_smooth(method = "lm")


cl_leafd15N_m1 <- lmer(CLratio ~ leaf_d15N * layer + (1|id), pyr_all_raw)
Anova(cl_leafd15N_m1, test.statistic = "F")
qqresidPlot(cl_leafd15N_m1)

# only fertilised plot
cl_leafd15N_fm1 <- lmer(CLratio ~ leaf_d15N * layer + (1|id), 
                        filter(pyr_all_raw, treatment == "fertilised"))
Anova(cl_leafd15N_fm1, test.statistic = "F")
qqresidPlot(cl_leafd15N_fm1)



# . carbohydrate:lignin vs. C mass ----------------------------------------

ggplot(pyr_all_raw, aes(x = CLratio, y = Cmass, col = layer))+
  geom_point() +
  geom_smooth(method = "lm")
pyr_all_raw2 <- pyr_all_raw %>% 
  mutate(id = factor(id))
Cmass_lcr_m1 <- lmer(Cmass ~ CLratio * layer + (1|id), pyr_all_raw)
Anova(Cmass_lcr_m1, test.statistic = "F")
qqresidPlot(Cmass_lcr_m1)
visreg(Cmass_lcr_m1, xvar = "CLratio", by = "layer", overlay = TRUE)

# only fertilised plot
ggplot(filter(pyr_all_raw, treatment == "fertilised"), 
       aes(x = CLratio, y = Cmass, color = layer))+
  geom_point()+
  geom_smooth(method = "lm")
Cmass_lcr_fm1 <- lmer(Cmass ~ CLratio + layer + (1|id), 
                      filter(pyr_all_raw, treatment == "fertilised"))
Anova(Cmass_lcr_fm1, test.statistic = "F")
qqresidPlot(Cmass_lcr_fm1)


# Save --------------------------------------------------------------------
save.image("Output/Data/Analysis_results.RData")

