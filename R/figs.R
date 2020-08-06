
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

save_png600("Output/Figs/Leaf_d15N_distance_all.png", width = 6, height = 3.5)
leafd15N_p
dev.off()


# . Transect --------------------------------------------------------------

plot_cor <- read.table("Data/rosinedal_corner_coords.txt", header = TRUE) %>% 
  filter(Square %in% c(2, 3)) %>% 
  mutate(treatment = mapvalues(Square, 2:3, c("Fertilised", "Control")))

# Define the point of origin
ori <- c(lon = min(transect_cor$lon), lat = min(transect_cor$lat))

# convert GPS to xy cordinate from the origin point
get_xy_cord <- function(p1, p2){
  # get angle and distance between the origin and each point
  distance <- distGeo(p1, p2)
  Rangle   <- bearing(p1, p2) * pi / 180 # degre to radian
  
  # conver to the xy surface
  x_lon <- distance * cos(pi/2 - Rangle)  # bearing = 0 when longitute of two points are the same
  y_lat <- distance * sin(pi/2 - Rangle)
  
  return(data.frame(distance, Rangle, x_lon, y_lat))
}


transect_xy <- adply(transect_cor, 1, function(x) 
  get_xy_cord(p1 = ori, p2 = c(x$lon[1], x$lat[1])))

plot_cor_xy <- adply(plot_cor, 1, function(x)
  get_xy_cord(p1 = ori, p2 = c(x$Long[1], x$Lat[1])))


site_map_p <- ggplot()+
  geom_polygon(data = plot_cor_xy, aes(x = x_lon, y = y_lat, group = treatment),
               fill = NA, colour = "gray50", size = .2)+
  geom_point(data = transect_xy, aes(x = x_lon, y = y_lat), size = 1, shape = 1, stroke = .4)+
  geom_segment(aes(x = 0, xend = 1000, y = 0, yend = 0), col = "black", size = .3) +
  geom_text(aes(x = c(0, 500, 1000), y = 0, label = c("0", "0.5", "1km")), size = 2, nudge_y = 55)+
  geom_point(aes(x = c(0, 500, 1000), y = 0), shape = 108, size = 3)+
  annotate("text", x = c(550, 2200), y = c(900, 400), label = c("Fertilised", "Control"), size = 3)+
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggsavePP("Output/Figs/Transect_map", plot = site_map_p, width = 4, height = 2)


# . Site map --------------------------------------------------------------


map_xlab <- c(expression(12*degree*E),
              expression(18*degree*E),
              expression(24*degree*E))

map_ylab <- c(expression(56*degree*N),
              expression(60*degree*N),
              expression(64*degree*N),
              expression(68*degree*N))

pdf("Output/Figs/Sweden_map.pdf", width = 1.5, height = 2)
# par(omi = c(0, 0, 0, 0))
map("worldHires","Sweden", xlim = c(10, 25), mar = c(1.5, 0, 0, 0))
axis(1, at = seq(12, 24, 6), labels = map_xlab, cex.axis = .5, tck = .02, padj = -4)
axis(2, at = seq(56, 68, 4), labels = map_ylab, las = 2, cex.axis = .5, tck = .02, hadj = .05)
points(mean(plot_cor$Long), mean(plot_cor$Lat), pch = 15)
box()
dev.off()


# Fig 3. Bar grph of composition of NMR and Pyrolysis products ------------

# Pyrolysis
pyr_prop_dd <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  select(category, compound, id, value, layer, location2) %>%
  mutate(layer = factor(layer, levels = c("Litter", "Humus")),
         horizon = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")),
         category = factor(category, levels = pyr_comp_ed)) %>% 
  group_by(category, layer, horizon, location2, id) %>% 
  summarise(value = sum(value) * 100) %>% 
  group_by(category, layer, horizon, location2) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

pyr_prop_p <- ggplot(pyr_prop_dd, aes(x = location2, y = value, fill = category))+
  geom_bar(stat = "identity", col = "black", size = .3) +
  facet_grid(. ~ horizon)+
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
  facet_grid(. ~ horizon)+
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
# nmr_prop_p

# combine
comp_prop_p <- ggarrange(pyr_prop_p, nmr_prop_p, nrow = 2, labels = c("(a)", "(b)"),
                         label.args = list(gp = grid::gpar(cex = 1), hjust = -.3, vjust = 1.5),
                         draw = FALSE)
# comp_prop_p
ggsavePP("Output/Figs/Composition_Pyr_NMR", plot = comp_prop_p,
         width = 6.5, height = 5)


# Fig 4. RDA --------------------------------------------------------------


# . Litter ----------------------------------------------------------------

# Pyrolysis
pyr_litter_rda_site <- pyr_litter %>% 
  mutate(PYR_RDA1 = as.vector(scores(pyr_litter_rda, display = "sites", choices = 1, scaling = 3)),
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, PYR_RDA1)
pyr_litter_rda_sp <-data.frame(scores(pyr_litter_rda, display = "species", choices = 1, scaling = 3)) %>% 
  mutate(pyr_comp = row.names(.)) %>% 
  arrange(RDA1)

# NMR
nmr_litter_rda_site <- nmr_litter %>% 
  mutate(NMR_RDA1 = as.vector(scores(nmr_litter_rda, display = "sites", choices = 1, scaling = 3)),
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, NMR_RDA1)
nmr_litter_rda_sp   <- data.frame(scores(nmr_litter_rda, display = "species", choices = 1, scaling = 3)) %>% 
  mutate(nmr_comp = row.names(.),
         nmr_comp2 = mapvalues(nmr_comp, 
                              c("alkylC" , "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                              c("Alkyl C", "Methoxy/N-\nalkyl C", "O-alkyl C", "Di-O-alkyl C", "Aromatic C", "O-aromatic C", "Carbonyl C"))) %>% 
  arrange(RDA1)

# Combine two
litter_rda_nmr_pyr <- left_join(nmr_litter_rda_site, pyr_litter_rda_site)

# define axis limits
litter_pyr_rda_val <- c(litter_rda_nmr_pyr$PYR_RDA1, pyr_litter_rda_sp$RDA1)
litter_nmr_rda_val <- c(litter_rda_nmr_pyr$NMR_RDA1, nmr_litter_rda_sp$RDA1)

litter_pyrlim <- c(min(litter_pyr_rda_val) * 1.01, max(litter_pyr_rda_val) * 1.01)
litter_nmrlim <- c(min(litter_nmr_rda_val) * 1.01, max(litter_nmr_rda_val) * 1.01)


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
litter_rda_nmrsp_p <- ggplot(nmr_litter_rda_sp, aes(x = 1, y = RDA1, label = nmr_comp2)) +
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
litter_rda_pyrsp_p <- ggplot(pyr_litter_rda_sp, aes(x = 1, y = RDA1, label = pyr_comp)) +
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
                                  ncol = 2, widths = c(2, .9), heights = c(.63, 2),
                                  draw = FALSE)
cairo_pdf(filename = "Output/Figs/RDA_Litter.pdf", width = 3.5, height = 3)
litter_rda_nmr_pyr_p
dev.off()



# . Humus ----------------------------------------------------------------

# Pyrolysis
pyr_humus_rda_site <- pyr_humus %>% 
  mutate(PYR_RDA1 = as.vector(scores(pyr_humus_rda, display = "site", choices = 1, scaling = 3)),
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, PYR_RDA1)
pyr_humus_rda_sp <- data.frame(scores(pyr_humus_rda, display = "species", choices = 1, scaling = 3)) %>% 
  mutate(pyr_comp = row.names(.)) %>% 
  arrange(RDA1)

# NMR
nmr_humus_rda_site <- nmr_humus %>% 
  mutate(NMR_RDA1 = as.vector(scores(nmr_humus_rda, display = "site", choices = 1, scaling = 3)),
         id       = as.character(id)) %>% 
  select(id, location2, leaf_d15N, NMR_RDA1)
nmr_humus_rda_sp   <- data.frame(scores(nmr_humus_rda, display = "species", choices = 1, scaling = 3)) %>% 
  mutate(nmr_comp = row.names(.),
         nmr_comp2 = mapvalues(nmr_comp, 
                               c("alkylC" , "meth_N_alkylC", "O_alkylC", "di_O_alkylC", "aromatic", "O_aromaticC", "carbonylC"),
                               c("Alkyl-C", "Methoxy/N-\nalkyl-C", "O-alkyl-C", "Di-O-alkyl-C", "Aromatic-C", "O-aromatic-C", "Carbonyl-C")))

# Combine two
humus_rda_nmr_pyr <- left_join(nmr_humus_rda_site, pyr_humus_rda_site)

# define axis limits
humus_pyr_rda_val <- c(humus_rda_nmr_pyr$PYR_RDA1, pyr_humus_rda_sp$RDA1)
humus_nmr_rda_val <- c(humus_rda_nmr_pyr$NMR_RDA1, nmr_humus_rda_sp$RDA1)
humus_pyrlim <- c(min(humus_pyr_rda_val) * 1.01, max(humus_pyr_rda_val) * 1.01)
humus_nmrlim <- c(min(humus_nmr_rda_val) * 1.01, max(humus_nmr_rda_val) * 1.01)


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
humus_rda_nmrsp_p <- ggplot(nmr_humus_rda_sp, aes(x = 1, y = RDA1, label = nmr_comp2)) +
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
humus_rda_pyrsp_p <- ggplot(pyr_humus_rda_sp, aes(x = 1, y = RDA1, label = pyr_comp)) +
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
                                 ncol = 2, widths = c(2, .9), heights = c(.63, 2),
                                 draw = FALSE)
cairo_pdf(filename = "Output/Figs/RDA_Humus.pdf", width = 3.5, height = 3)
humus_rda_nmr_pyr_p
dev.off()

save_png600(filename = "Output/Figs/RDA_Humus.png", width = 3.5, height = 3)
humus_rda_nmr_pyr_p
dev.off()

# . Marge litter and humus ------------------------------------------------
rda_nmr_pyr_p <- ggarrange(litter_rda_pyrsp_p, blank_ggplot, 
                           humus_rda_pyrsp_p, blank_ggplot, 
                           litter_rda_p, litter_rda_nmrsp_p, 
                           humus_rda_p, humus_rda_nmrsp_p, 
                           ncol = 4, widths = c(2, .9, 2, .9), heights = c(.63, 2),
                           labels = c("(a) L horizon", "", "(b) F/H horizon", rep("", 5)),
                           label.args = list(gp = grid::gpar(cex = 1), hjust = -.1, vjust = 1.5),
                           draw = FALSE)
cairo_pdf(filename = "Output/Figs/RDA_Pyr_NMR2.pdf", width = 6.5, height = 3)
rda_nmr_pyr_p
dev.off()





# Fig 5. Lignin:carbohydrate ratios ---------------------------------------


lc_newd <- with(pyr_all_raw,
                expand.grid(leaf_d15N = seq(min(leaf_d15N), max(leaf_d15N), length.out = 1000),
                            layer = unique(layer)))
lc_pred <- lc_newd %>% 
  mutate(LCratio = 1/predict(lc_leafd15N_m3, lc_newd, re.form=NA),
         horizon  = mapvalues(layer, c("Litter", "Humus"), paste(c("L", "F/H"), "horizon")))

# Log
LC_leafd15N_P <- ggplot(pyr_all_raw, aes(x = leaf_d15N, y = log(LCratio)))+
  geom_point(aes(col = horizon), size = 2, alpha = .8) +
  geom_line(data = lc_pred, aes(x = leaf_d15N, y = log(LCratio), col = horizon))+
  scale_color_manual(values = c("black", "red"))+
  labs(x = expression(paste("Leaf ", delta^{15}, "N (\u2030)")),
       y = expression(log[e](Lignin:Carbohydrate~ratio)))+
  science_theme+
  theme(legend.key.width = unit(2, "lines"),
        legend.position = c(.2, .6))
cairo_pdf("Output/Figs/Lignin_Carbohydrate_ratio_leafd15N_log.pdf", width = 4, height = 3)
LC_leafd15N_P
dev.off()

cairo_ps("Output/Figs/Lignin_Carbohydrate_ratio_leafd15N_log.eps", width = 4, height = 3)
LC_leafd15N_P
dev.off()

# raw
LC_leafd15N_raw_P <- ggplot(pyr_all_raw, aes(x = leaf_d15N, y = LCratio))+
  geom_point(aes(col = horizon), size = 2, alpha = .8) +
  geom_line(data = lc_pred, aes(x = leaf_d15N, y = LCratio, col = horizon))+
  scale_color_manual(values = c("black", "red"))+
  labs(x = expression(paste("Leaf ", delta^{15}, "N (\u2030)")),
       y = "Lignin:Carbohydrate ratio")+
  science_theme+
  theme(legend.key.width = unit(2, "lines"),
        legend.position = c(.2, .6))
cairo_pdf("Output/Figs/Lignin_Carbohydrate_ratio_leafd15N.pdf", width = 4, height = 3)
LC_leafd15N_raw_P
dev.off()

cairo_ps("Output/Figs/Lignin_Carbohydrate_ratio_leafd15N.eps", width = 4, height = 3)
LC_leafd15N_raw_P
dev.off()


LC_Cmass_P <- ggplot(pyr_all_raw, aes(x = log(LCratio), y = Cmass))+
  geom_point(aes(col = layer), size = 2, alpha = .8) +
  geom_smooth(aes(col = layer), method = "lm", se = FALSE, size = 1)+
  scale_color_manual(values = c("black", "red"))+
  science_theme+
  labs(x = "log(Lignin:Carbohydrate ratio)",
       y = expression(C~mass~(kg~C~ha^'-1')))+
  theme(legend.position = "NULL")

LC_p <- ggarrange(LC_leafd15N_P, LC_Cmass_P, ncol = 2, labels = c("(a)", "(b)"),
                  label.args = list(gp = grid::gpar(cex = 1), hjust = -.3, vjust = 1.5),
                  draw = FALSE)
cairo_pdf("Output/Figs/Lignin_Carbohydrate_ratio.pdf", width = 6.5, height = 3)
LC_p
dev.off()

save_png600("Output/Figs/Lignin_Carbohydrate_ratio.png", width = 6.5, height = 3)
LC_p
dev.off()



# Fig S. N comp -----------------------------------------------------------

py_Ncomp_prop <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  filter(category == "N-compound") %>% 
  select(location2, category, compound, id, value, horizon) %>% 
  group_by(id, horizon) %>% 
  mutate(Nprop = value * 100 / sum(value)) %>% 
  group_by(location2, compound, horizon) %>% 
  summarise(value = mean(Nprop)) %>% 
  ungroup()
ggplot(py_Ncomp_prop, aes(x  = location2, y = value, fill = compound))+
  geom_bar(stat = "identity")+
  facet_grid(. ~ horizon)+
  labs(y = "Proportion (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 6))




# Fig S. Pyrolysates ------------------------------------------------------

pyr_prop_smmry <- rbind.fill(pyr_litter_spec, pyr_humus_spec) %>% 
  group_by(horizon, category, compound, location2) %>% 
  summarise_at(.vars = vars(value), .funs = funs(M = mean, CI = get_ci, N = get_n)) %>% 
  ungroup() %>% 
  mutate(compound = factor(compound, levels = sort(unique(compound), decreasing = TRUE)),
         location2 = mapvalues(location2,
                               c("control", "fertilised:outside", "fertilised:edge","fertilised:inside"),
                               c("Control", "Fertilised:Outside", "Fertilised:Edge","Fertilised:Inside")))

pyr_prop_cat_p <- dlply(pyr_prop_smmry, .(category), function(z){
  gp <- ggplot(z, aes(x = compound, y = M, col = location2))+
    geom_errorbar(aes(ymin = M - CI, ymax = M + CI, width = 0), alpha = .6,
                  position = position_dodge(.3))+
    geom_point(alpha = .6, size = .9, position = position_dodge(.3))+
    facet_grid(. ~ horizon)+
    coord_flip()+
    labs(x = NULL, y = "Proportion (%)")+
    science_theme+
    theme(axis.text.x        = element_text(size = 6.5),
          axis.text.y        = element_text(size = 6.5),
          axis.title.x       = element_text(size = 6.7),
          strip.text.x       = element_text(size = 6.5,
                                            margin = margin(.05,0,.05,0, "cm")),
          legend.title       = element_blank(),
          legend.position    = "right",
          panel.grid.major   = element_line(colour = "grey90", size = .2),
          panel.grid.major.x = element_blank())
  return(gp)
})

# No. of compounds for each group
pyr_fig_length <- pyr_prop_smmry %>% 
  select(category, compound) %>% 
  distinct() %>% 
  group_by(category) %>% 
  summarise(No = length(compound)) %>% 
  mutate(No_r = No / max(No),
         No_r = ifelse(No_r == 1, No_r, No_r + .05))

# edit
pyr_prop_cat_p$`N-compound` <- pyr_prop_cat_p$`N-compound`+
  scale_y_continuous(breaks = c(0, 0.005, 0.01), labels = c(0, 0.005, 0.01))

# save
l_ply(1:8, function(x){
  p <- pyr_prop_cat_p[[x]]
  p_name <- names(pyr_prop_cat_p)[x]
  p_height_r <- pyr_fig_length$No_r[x]
  ggsavePP(filename = paste0("Output/Figs/Pyrolysates/Pyrolysate_", p_name),
           p, width = 6.5, height = 6.5 * p_height_r)
})

