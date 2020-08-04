
# Pyrolysis raw data
pyr_litter_raw_data <- pyr_litter_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         category = ifelse(category %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                           ifelse(category %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                                  as.character(category)))) %>% 
  select(id, horizon, category, compound, area) %>% 
  arrange(horizon, as.numeric(as.character(id)))
write.csv(pyr_litter_raw_data, "Output/Data/Manuscript/Pyrolysis_chromatogram_Lhoriz.csv", row.names = FALSE)


pyr_humus_raw_data <- pyr_humus_raw_area %>% 
  filter(id %in% use_id) %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         category = ifelse(category %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "alphatic_derivative", 
                           ifelse(category %in% c("chlorophyll", "vitamin", "steroid", "hopanoid"), "Others",
                                  as.character(category)))) %>% 
  select(id, horizon, category, compound, area) %>%
  rename(peak_area = area) %>% 
  arrange(horizon, as.numeric(as.character(id)))
write.csv(pyr_humus_raw_data, "Output/Data/Manuscript/Pyrolysis_chromatogram_FHhoriz.csv", row.names = FALSE)

# NMR raw data
nmr_intg_raw_data <- nmr_intg_raw %>% 
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon"))) %>% 
  select(id, horizon, one_of(NMR_compounds)) %>% 
  arrange(horizon, as.numeric(as.character(id)))
write.csv(pyr_humus_raw_data, "Output/Data/Manuscript/NMR_spectra_relative_abund.csv", row.names = FALSE)

# IRMS data
irms_data <- irms_dd %>% 
  select(id, horizon, treatment, location2, Cmass, Nmass, CNratio) %>% 
  rename(location = location2) %>% 
  arrange(horizon, as.numeric(as.character(id)))
write.csv(irms_data, "Output/Data/Manuscript/Soil_CN.csv", row.names = FALSE)

# lingonberry leaf d15N
irms_lingon_ful_data <- irms_lingon_ful %>% 
  select(id, treatment, location2, distance, leaf_d15N) %>% 
  rename(location = location2)
write.csv(irms_lingon_ful_data, "Output/Data/Manuscript/Leaf_d15N.csv", row.names = FALSE)
