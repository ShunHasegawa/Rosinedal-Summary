
# Fig 1. Hiistory of N addition -------------------------------------------
Nadd_d <- data.frame(year = 2005:2019,
                     Nadd = c(0, rep(100, 6), rep(50, 7), 0)) %>% 
  mutate(total_addN = cumsum(Nadd))
Nadd_p <- ggplot(Nadd_d, aes(x = year, y = total_addN)) +
  geom_step(size = 1) +
  scale_x_continuous(breaks = 2006:2018, labels = 2006:2018) +
  labs(x = "", y = expression(Total~added~N~(kg~ha^'-1'))) +
  annotate("text", x = c(2006, 2012, 2017), y = c(100, 650, 950), label = c("+100", "+50", "Total: +950"), 
           col = "black", vjust = -.4, size = 5) +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13), 
        axis.title  = element_text(size = 15),
        panel.border = element_rect(size = 1.5)) +
  lims(y= c(0, 1100)) 
Nadd_p
ggsavePP("Output/Figs/Nadd_history", plot = Nadd_p, width = 5, height = 4)

# Fig 2. Site map and leaf d15N -------------------------------------------


# . Leaf d15N -------------------------------------------------------------
irms_lingon_trans_f <- filter(irms_lingon_ful, treatment == "fertilised") %>% 
  mutate(distance3 = distance + 100)
iv <- getInitial(leaf_d15N ~ SSweibull(distance3, Asym, Drop, lrc, pwr), data = irms_lingon_trans_f) # starting values
nm1 <- nls(leaf_d15N ~ Asym - Drop * exp(-exp(lrc) * distance3^pwr), start = iv, data = irms_lingon_trans_f)
summary(nm1)
# nf <- expression(y=-5.1+5.6*exp(-exp(-21.7)*(x-100)^4.9))
# plot(1, main = nf)
predd <- data.frame(distance3 = seq(0, 200, length.out = 100)) %>% 
  mutate(leaf_d15N = predict(nm1, list(distance3 = distance3)),
         distance  = distance3 - 100, 
         treatment = "Fertilised")
trt_lab_dc <- data.frame(distance = c(-65, 0, 65),
                         treatment = "Fertilised",
                         leaf_d15N = min(irms_lingon_ful$leaf_d15N),
                         location = c("Inside", "Edge", "Outside"))
irms_lingon_ful_ed <- irms_lingon_ful %>% 
  mutate(treatment = mapvalues(treatment, c("control", "fertilised"), c("Control", "Fertilised")))
leafd15N_p <- ggplot(irms_lingon_ful_ed, aes(x = distance, y = leaf_d15N))+
  geom_vline(data = filter(irms_lingon_ful_ed, treatment == "Fertilised"),
             aes(xintercept = -20), linetype = "dotted", col = "gray") +
  geom_vline(data = filter(irms_lingon_ful_ed, treatment == "Fertilised"),
             aes(xintercept = 20), linetype = "dotted", col = "gray") +
  geom_point(size = 2, alpha = .8, color = "gray30") +
  geom_line(data = predd, aes(x = distance, y = leaf_d15N), size = 1, col = "black")+
  labs(x = "Distance from the edge of plot (m)", y = expression(paste("Leaf ", delta^{15}, "N (\u2030)")))+
  theme(legend.position = c(.1, .15),
        legend.title    = element_blank())+
  facet_grid(. ~ treatment, scales = "free_x", space = "free_x") +
  geom_text(data = trt_lab_dc, aes(x = distance, y = leaf_d15N, label = location),
            fontface = "italic", nudge_y = -.3, size = 3)
leafd15N_p
cairo_pdf("Output/Figs/Leaf_d15N_distance_all.pdf", width = 6, height = 3.5)
leafd15N_p
dev.off()

# Fig 3. Bar grph of composition of NMR and Pyrolysis products ------------

# Pyrolysis
pyr_prop_dd <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(grp, comp, id, value, layer, location2) %>%
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         grp2 = mapvalues(grp, 
                          c("alphatic_derivative", "aromatic", "carbohydrate", 
                            "g_lignin", "N_comp", "Phenol", "s_lignin", "Others"),
                          pyr_comp_ed),
         grp2 = factor(grp2, levels = pyr_comp_ed)) %>% 
  group_by(grp2, layer, location2, id) %>% 
  summarise(value = sum(value) * 100) %>% 
  group_by(grp2, layer, location2) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

pyr_prop_p <- ggplot(pyr_prop_dd, aes(x = location2, y = value, fill = grp2))+
  geom_bar(stat = "identity", col = "black", size = .3) +
  facet_grid(. ~ layer)+
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666")) +
  labs(x = NULL, y = "Pyrolysis product\nproportion (%)") +
  science_theme+
  theme(legend.spacing.x  = unit(.1, "lines"),
        legend.key.width  = unit(1, "lines"),
        legend.position = "right",
        axis.text.x = element_blank(),
        legend.title = element_blank())

# NMR
nmr_intg_smmry$grp <- mapvalues(nmr_intg_smmry$compound, 
                                c("alkylC", "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                                nmr_comp_ed)
nmr_prop_p <- ggplot(nmr_intg_smmry, aes(x = location2, y = prop, fill = grp))+
  geom_bar(stat = "identity", col = "black", size = .3) +
  facet_grid(. ~ layer)+
  scale_fill_manual(values = c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"),
                    labels = c("Alkyl C", 
                               "Methoxy/N-alkyl C", 
                               expression(italic(O)*-alkyl~C),
                               expression(Di*-italic(O)*-alkyl~C), 
                               "Aromatic C", 
                               expression(italic(O)*-aromatic~C), "Carbonyl C")) +
  labs(x = NULL, y = "NMR product\nproortion (%)") +
  science_theme+
  theme(legend.spacing.x  = unit(.1, "lines"),
        legend.key.width  = unit(1, "lines"),
        legend.position   = "right",
        axis.text.x       = element_text(angle = 45, hjust = 1),
        strip.background  = element_blank(),
        strip.text.x      = element_blank(),
        legend.title      = element_blank(),
        legend.text.align = 0)
nmr_prop_p

# combine
comp_prop_p <- ggarrange(pyr_prop_p, nmr_prop_p, nrow = 2, labels = c("(a)", "(b)"),
                         label.args = list(gp = grid::gpar(cex = 1), hjust = -.3, vjust = 1.5))
comp_prop_p
ggsavePP("Output/Figs/Composition_Pyr_NMR", plot = comp_prop_p,
         width = 6.5, height = 5)


# Fig 4. RDA --------------------------------------------------------------


# . Litter ----------------------------------------------------------------

# Pyrolysis
pyr_litter_rda_site <- pyr_litter %>% 
  mutate(PYR_RDA1 = data.frame(scores(pyr_litter_rda)$sites)[, 1],
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, PYR_RDA1)
pyr_litter_rda_sp <- data.frame(scores(pyr_litter_rda)$species) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
         pyr_comp2 = mapvalues(pyr_comp, 
                   c("alphatic_derivative", "aromatic", "carbohydrate", 
                     "g_lignin", "N_comp", "Phenol", "s_lignin", "Others"),
                   pyr_comp_ed))

# NMR
nmr_litter_rda_site <- nmr_litter %>% 
  mutate(NMR_RDA1 = data.frame(scores(nmr_litter_rda)$sites)[, 1],
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, NMR_RDA1)
nmr_litter_rda_sp   <- data.frame(scores(nmr_litter_rda)$species) %>% 
  mutate(nmr_comp = row.names(.),
         nmr_comp = factor(nmr_comp, levels = nmr_comp[order(RDA1)]),
         nmr_comp2 = mapvalues(nmr_comp, 
                              c("alkylC" , "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                              c("Alkyl C", "Methoxy/N-\nalkyl C", "O-alkyl C", "Di-O-alkyl C", "Aromatic C", "O-aromatic C", "Carbonyl C")))

# Combine two
litter_rda_nmr_pyr <- left_join(nmr_litter_rda_site, pyr_litter_rda_site)

litter_pyrlim <- with(litter_rda_nmr_pyr, c(min(PYR_RDA1) * 1.01, max(PYR_RDA1) * 1.01))
litter_nmrlim <- with(litter_rda_nmr_pyr, c(min(NMR_RDA1) * 1.01, max(NMR_RDA1) * 1.01))


# NMR RDA1 vs Pyr RDA1
litter_rda_p <- ggplot(litter_rda_nmr_pyr, aes(x = PYR_RDA1, y = NMR_RDA1, fill = leaf_d15N))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(size = 3, alpha = .9, shape = 21)+
  science_theme +
  scale_fill_gradient2(name = expression(paste("Leaf ", delta^{15}, "N (\u2030)")),
                       low = "#F7FF00", high = "#8C0000", mid = "#DC4600", midpoint = -2.5)+
  lims(x = litter_pyrlim, y = litter_nmrlim) +
  labs(x = paste("Pyr", get_PCA_axislab(pyr_litter_rda)[1]),
       y = paste("NMR", get_PCA_axislab(nmr_litter_rda)[1]))+
  theme(legend.title = element_text(size = 7),
        legend.text =  element_text(size = 5),
        legend.key.width = unit(.6, units = "line"),
        legend.key.height = unit(.5, units = "line"),
        legend.position = c(.2, .8),
        legend.spacing.y = unit(0, "pt"),
        axis.title = element_text(size = 9),
        axis.text =  element_text(size = 7),
        plot.margin = margin(t = 0, r = 0, b = .5, l = .5, unit = "line")) 

# sp loading for NMR
litter_rda_nmrsp_p <- ggplot(nmr_litter_rda_sp, aes(x = 1, y = RDA1 * 3, label = nmr_comp2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = .92) +
  geom_text(aes(x = .92), size = 3, label = "-", hjust = 0) +
  geom_text(hjust = 0, size = 2.5, lineheight = .7, nudge_x = .1) +
  lims(x = c(.87, 2.77), y = litter_nmrlim) +
  labs(x = NULL, y = NULL) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = .5, l = 0, unit = "line"))

# sp loading for Pyrolysis
litter_rda_pyrsp_p <- ggplot(pyr_litter_rda_sp, aes(x = 1, y = RDA1 * 1.3, label = pyr_comp2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = .92) +
  geom_text(aes(x = .92), size = 3, label = "-", hjust = 0, angle = 90) +
  geom_text(angle = 45, hjust = 0, size = 2.5, nudge_x = .1) +
  labs(x = NULL, y = NULL) +
  lims(x = c(.87, 2.17), y = litter_pyrlim) +
  science_theme +
  coord_flip() +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = .5, unit = "line"))
# arrange(pyr_litter_rda_sp, RDA1)

# blank plot
blank_ggplot <- ggplot() + 
  geom_blank() +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "line"))

# merge the above plots
litter_rda_nmr_pyr_p <- ggarrange(litter_rda_pyrsp_p, blank_ggplot, litter_rda_p, litter_rda_nmrsp_p,  
                                  ncol = 2, widths = c(2, .9), heights = c(.63, 2))
cairo_pdf(filename = "Output/Figs/RDA_Litter.pdf", width = 3.5, height = 3)
litter_rda_nmr_pyr_p
dev.off()



# . Humus ----------------------------------------------------------------

# Pyrolysis
pyr_humus_rda_site <- pyr_humus %>% 
  mutate(PYR_RDA1 = data.frame(scores(pyr_humus_rda)$sites)[, 1],
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, PYR_RDA1)
pyr_humus_rda_sp <- data.frame(scores(pyr_humus_rda)$species) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
         pyr_comp2 = mapvalues(pyr_comp, 
                               c("alphatic_derivative", "aromatic", "carbohydrate", 
                                 "g_lignin", "N_comp", "Phenol", "s_lignin", "Others"),
                               pyr_comp_ed))

# NMR
nmr_humus_rda_site <- nmr_humus %>% 
  mutate(NMR_RDA1 = data.frame(scores(nmr_humus_rda)$sites)[, 1],
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, NMR_RDA1)
nmr_humus_rda_sp   <- data.frame(scores(nmr_humus_rda)$species) %>% 
  mutate(nmr_comp = row.names(.),
         nmr_comp = factor(nmr_comp, levels = nmr_comp[order(RDA1)]),
         nmr_comp2 = mapvalues(nmr_comp, 
                               c("alkylC" , "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                               c("Alkyl-C", "Methoxy/N-\nalkyl-C", "O-alkyl-C", "Di-O-alkyl-C", "Aromatic-C", "O-aromatic-C", "Carbonyl-C")))

# Combine two
humus_rda_nmr_pyr <- left_join(nmr_humus_rda_site, pyr_humus_rda_site)

# pyrlim <- with(humus_rda_nmr_pyr, c(min(PYR_RDA1) * 1.01, max(PYR_RDA1) * 1.03))
humus_pyrlim <- with(humus_rda_nmr_pyr, c(min(PYR_RDA1) * 1.01, max(PYR_RDA1) * 1.01))
humus_nmrlim <- with(humus_rda_nmr_pyr, c(min(NMR_RDA1) * 1.01, max(NMR_RDA1) * 1.01))


# NMR RDA1 vs Pyr RDA1
humus_rda_p <- ggplot(humus_rda_nmr_pyr, aes(x = PYR_RDA1, y = NMR_RDA1, fill = leaf_d15N))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(size = 3, alpha = .9, shape = 21)+
  science_theme +
  scale_fill_gradient2(name = expression(paste("Leaf ", delta^{15}, "N (\u2030)")),
                       low = "#F7FF00", high = "#8C0000", mid = "#DC4600", midpoint = -2.5)+
  lims(x = humus_pyrlim, y = humus_nmrlim) +
  labs(x = paste("Pyr", get_PCA_axislab(pyr_humus_rda)[1]),
       y = paste("NMR", get_PCA_axislab(nmr_humus_rda)[1]))+
  theme(legend.position = "none",
        axis.title = element_text(size = 9),
        axis.text =  element_text(size = 7),
        plot.margin = margin(t = 0, r = 0, b = .5, l = .5, unit = "line")) 

# sp loading for NMR
humus_rda_nmrsp_p <- ggplot(nmr_humus_rda_sp, aes(x = 1, y = RDA1 * 1.8, label = nmr_comp2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = .92) +
  geom_text(aes(x = .92), size = 3, label = "-", hjust = 0) +
  geom_text(hjust = 0, size = 2.5, lineheight = .7, nudge_x = .1) +
  lims(x = c(.87, 2.77), y = humus_nmrlim) +
  labs(x = NULL, y = NULL) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = .5, l = 0, unit = "line"))

# sp loading for Pyrolysis
humus_rda_pyrsp_p <- ggplot(pyr_humus_rda_sp, aes(x = 1, y = RDA1 * 1.1, label = pyr_comp2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = .92) +
  geom_text(aes(x = .92), size = 3, label = "-", hjust = 0, angle = 90) +
  geom_text(angle = 45, hjust = 0, size = 2.5, nudge_x = .1) +
  labs(x = NULL, y = NULL) +
  lims(x = c(.87, 2.17), y = humus_pyrlim) +
  science_theme +
  coord_flip() +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = .5, unit = "line"))
# arrange(pyr_humus_rda_sp, RDA1)

# merge the above plots
humus_rda_nmr_pyr_p <- ggarrange(humus_rda_pyrsp_p, blank_ggplot, humus_rda_p, humus_rda_nmrsp_p,  
                                 ncol = 2, widths = c(2, .9), heights = c(.63, 2))
cairo_pdf(filename = "Output/Figs/RDA_Humus.pdf", width = 3.5, height = 3)
humus_rda_nmr_pyr_p
dev.off()


# . Marge litter and humus ------------------------------------------------
rda_nmr_pyr_p <- ggarrange(litter_rda_pyrsp_p, blank_ggplot, 
                           humus_rda_pyrsp_p, blank_ggplot, 
                           litter_rda_p, litter_rda_nmrsp_p, 
                           humus_rda_p, humus_rda_nmrsp_p, 
                           ncol = 4, widths = c(2, .9, 2, .9), heights = c(.63, 2),
                           labels = c("(a) Litter", "", "(b) Humus", rep("", 5)),
                           label.args = list(gp = grid::gpar(cex = 1), hjust = -.1, vjust = 1.5))
cairo_pdf(filename = "Output/Figs/RDA_Pyr_NMR_v2.pdf", width = 6.5, height = 3)
rda_nmr_pyr_p
dev.off()





# Fig 5. Carbohydrate:Lignin ratios ---------------------------------------

CL_leafd15N_P <- ggplot(pyr_all_raw, aes(x = leaf_d15N, y = CLratio))+
  geom_point(aes(col = layer), size = 3, alpha = .8) +
  geom_smooth(aes(col = layer), method = "lm", se = FALSE, size = 1)+
  scale_color_manual(values = c("black", "red"))+
  labs(x = expression(paste("Leaf ", delta^{15}, "N (\u2030)")),
       y = "Carbohydrate:Lignin ratio")+
  science_theme+
  theme(legend.position = c(.2, .1))

CL_Cmass_P <- ggplot(pyr_all_raw, aes(x = CLratio, y = Cmass))+
  geom_point(aes(col = layer), size = 3, alpha = .8) +
  geom_smooth(aes(col = layer), method = "lm", se = FALSE, size = 1)+
  scale_color_manual(values = c("black", "red"))+
  science_theme+
  labs(x = "Carbohydrate:Lignin ratio",
       y = expression(C~mass~(kg~C~ha^'-1')))+
  theme(legend.position = "NULL")

CL_p <- ggarrange(CL_leafd15N_P, CL_Cmass_P, ncol = 2, labels = c("(a)", "(b)"),
                  label.args = list(gp = grid::gpar(cex = 1), hjust = -.3, vjust = 1.5))
cairo_pdf("Output/Figs/Carbohydrate_Lignin_ratio.pdf", width = 6.5, height = 3)
CL_p
dev.off()



# Fig S. N comp -----------------------------------------------------------

py_Ncomp_prop <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  filter(grp == "N_comp") %>% 
  select(location2, grp, comp, variable, value, layer, variable) %>% 
  group_by(variable) %>% 
  mutate(Nprop = value * 100 / sum(value)) %>% 
  group_by(location2, comp, layer) %>% 
  summarise(value = mean(Nprop)) %>% 
  ungroup()
ggplot(py_Ncomp_prop, aes(x  = location2, y = value, fill = comp))+
  geom_bar(stat = "identity")+
  facet_grid(. ~ layer)+
  labs(y = "Proportion (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 6))

