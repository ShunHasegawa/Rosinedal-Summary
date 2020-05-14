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

# map("worldHires","Sweden", col="gray30", fill=TRUE, xlim = c(10, 25)) 
map_xlab <- c(expression(12*degree*E),
              expression(18*degree*E),
              expression(24*degree*E))

map_ylab <- c(expression(56*degree*E),
              expression(60*degree*E),
              expression(64*degree*E),
              expression(68*degree*E))

pdf("Output/Sweden_map.pdf", width = 1.5, height = 2)
# par(omi = c(0, 0, 0, 0))
map("worldHires","Sweden", xlim = c(10, 25), mar = c(1.5, 0, 0, 0))
axis(1, at = seq(12, 24, 6), labels = map_xlab, cex.axis = .5, tck = .02, padj = -4)
axis(2, at = seq(56, 68, 4), labels = map_ylab, las = 2, cex.axis = .5, tck = .02, hadj = -.05)
points(mean(plot_cor$Long), mean(plot_cor$Lat), pch = 15)
box()
dev.off()
