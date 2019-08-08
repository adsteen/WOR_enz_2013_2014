## Creates Fig 2: sum Vmax relative to surface, sum vmax per cell, GRR, and sum Vmax per GRR. It relies on an external script
## Must run 2013-06-24 calculating activities DS version.R first

######
# Merge in cell data
######

sum_kinetics <- kinetics %>%
  left_join(cell_counts, by = "depth") %>%
  group_by(depth, quant.depth, cell.ml, cell.ml.sd) %>%
  summarise(sum.Vmax = sum(Vmax, na.rm = TRUE),
            sum.Vmax.err = sqrt(sum(Vmax.err^2, na.rm = TRUE)),
            cell.norm.Vmax = sum.Vmax / cell.ml[1] * 2.5 * 10^9, # 2.5 g sed per ml converts cell densities (per ml) to enzyme activities (per g sed)
            #cell.norm.Vmax.err = cell.norm.Vmax[1] * sqrt((sum.Vmax.err/sum.Vmax)^2 + (cell.ml.sd[1]/cell.ml[1])^2) * 2.5 * 10^9) # the 10^9 takes us from umol to 
            cell.norm.Vmax.err = cell.norm.Vmax[1] * sqrt((sum.Vmax.err/sum.Vmax)^2 + (cell.ml.sd[1]/cell.ml[1])^2) ) # the 10^9 takes us from umol to 

# Calculate mean cell-specific sum Vmax and trend
mean(sum_kinetics$cell.norm.Vmax)
sd(sum_kinetics$cell.norm.Vmax)
sum_Vmax_lm <- lm(log10(sum.Vmax) ~ quant.depth, data=sum_kinetics)
summary(sum_Vmax_lm)
cell_specific_Vmax_lm <- lm(cell.norm.Vmax ~ quant.depth, data = sum_kinetics)
summary(cell_specific_Vmax_lm)
Km_mod <- lm(log10(Km) ~ quant.depth, data=kinetics)

# Weirdly sum_vmax dissapeared along the way
sum_kinetics <- sum_kinetics %>%
  ungroup() %>%
  mutate(sum.Vmax.rel = sum.Vmax / max(sum.Vmax),
         sum.Vmax.rel.err = sum.Vmax.err / max(sum.Vmax)) 

# Fig 2a: Relative sum vmax vs depth
p_Vmax.rel <- ggplot(sum_kinetics, aes(x=quant.depth)) +
  geom_point(aes(y=sum.Vmax.rel)) + 
  geom_line(aes(y=sum.Vmax.rel)) +
  geom_errorbar(aes(ymin = sum.Vmax.rel-sum.Vmax.rel.err, ymax = sum.Vmax.rel + sum.Vmax.rel.err)) +
  scale_x_reverse() + 
  scale_y_continuous(limits = c(0, 1.2),
                     breaks = seq(from = 0, to = 1, by = 0.25),
                     labels = percent) + 
  xlab("depth, cmbsf") +
  ylab(expression(paste(sum(V[max]), ", relative to surface"))) +
  coord_flip() + 
  theme(axis.text.x=element_text(angle=-45, hjust=0))
if(print.plots) {
  print(p_Vmax.rel)
}

# Per-cell Vmax
p_Vmax_cell_norm <- ggplot(sum_kinetics, aes(x=quant.depth)) + 
  geom_point(aes(y=cell.norm.Vmax)) + 
  geom_line(aes(y=cell.norm.Vmax)) + 
  geom_errorbar(aes(ymin = cell.norm.Vmax - cell.norm.Vmax.err, ymax = cell.norm.Vmax + cell.norm.Vmax.err)) + 
  #xlab("depth, cmbsf") +
  ylab(expression(paste(sum(V[max]), " per cell, amol ", cell^{-1}, " ",hr^{-1}))) +
  scale_x_reverse() + 
  expand_limits(ymin=0) +
  coord_flip() + 
  theme(axis.title.y = element_blank())
if(print.plots) {
  print(p_Vmax_cell_norm)
}



#(kG1 = 0.12 y-1; avg. lifetime = 8 y) accounted for 73.5% of the reactive OC, and the less reactive fraction (kG2 = 0.011 y-1; avg. lifetime = 91 y)
# From Alperin: kG1 = 0.12 y-1; kG2 = 0.011 y-1

GRR <- read.table("data/model results/GRR.mod", header=FALSE, col.names=c("depth", "G")) 
attr(GRR$G, "units") <- "uM wet sediment per day"

# Convert G to units of umol wet sediment per liter per hour
GRR$G.per.hr <- GRR$G * 24
attr(GRR$G, "units") <- "umol wet sediment per liter per hour"

# Figure 2c ##### THE ISSUE IS THAT ANNOTATION LOGTICKS AND COORD_FLIP ARE INCOMPATIBLE
p_GRR <- ggplot(GRR, aes(x=depth, y=G.per.hr)) + 
  geom_line() + 
  scale_y_log10(name=expression(paste("OC oxidation rate, ", mu, "M ", hr^{-1})), breaks=c(10, 100, 1000, 10000)) + 
  scale_x_reverse(name="depth, cmbsf", limits=c(85, 0)) + 
  #annotation_logticks(sides="l") + 
  coord_flip() +
  #xlab(expression(paste("OC oxidation rate, ", mu, "M ", hr^{-1}))) +
  #ylab("depth, cmbsf") + 
  theme(axis.title.y = element_blank())
if(print.plots) {
  print(p_GRR)
}


######
# Plot Vmax relative to G
######

# Interpolate the model output to my precise samplign depths
# Note these are exact, because I am using the same functional form as the model
G_form <- formula(G.per.hr ~ A1 * exp(-k1*depth) + A2 * exp(-k2*depth))
G_mod <- nls(G_form, GRR, start=list(A1=4800, k1=0.1, A2=2400, k2=0.001))
pred_grid <- data.frame(depth=seq(from=4.5, to=100, by=0.1))
G_preds <- data.frame(depth=pred_grid$depth, G.per.hr=predict(G_mod, newdata=pred_grid))

# # Check that my fits worked right. They did.
# ggplot() + 
#   geom_line(data=G_preds, aes(x=depth, y=G.per.hr)) + 
#   geom_point(data=GRR, aes(x=depth, y=G.per.hr))
 # THIS IS ALL F'd UP
Vmax_and_G <- left_join(x=sum_kinetics, y=G_preds, by = c("quant.depth" = "depth"))
#Vmax_and_G <- merge(x=sum_Vmax, y=G_preds, by.x="quantDepth", by.y="depth")
Vmax_and_G <- Vmax_and_G %>%
  mutate(Vmax.rel.G = sum.Vmax / (G.per.hr/1.286/1000), # G per hr needs to be *DIVIDED* by 1000 to go from per liter to per ml
                     Vmax.rel.G.err = sum.Vmax.err / (G.per.hr/1.286/1000)) # Check units
# G per hour is in umol C per liter whole sed whole sed per hour; Vmax is in umol per g wet sed per hour
# The bulk density of the sediment is 1.286 g / ml, so I can divide the G by 1.286*1000 to get a better ratio
attr(Vmax_and_G$G.per.hr, "units") <- "umol C per g wet sed per hour"


# Fig 2d
# Figure 2d
p_Vmax_rel_G <- ggplot(Vmax_and_G, aes(x=quant.depth, y=Vmax.rel.G)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=Vmax.rel.G-Vmax.rel.G.err, ymax=Vmax.rel.G+Vmax.rel.G.err)) + 
  geom_line() + 
  scale_x_reverse(limits=c(85, 0)) + 
  xlab("depth, cmbsf") + 
  ylab(expression(paste(sum(V[max]), ":OC oxidation rate"))) + 
  coord_flip() + 
  theme(axis.title.y = element_blank())
if(print.plots) {
  print(p_Vmax_rel_G)
}

# # Error prop for x = ab/c
# # sigma_x  = x * sqrt((sigma_a / a)^2 + (sigma_b/b)^2 + (sigma_c/c)^2)
# # Here, cell.norm.G  = G / cells. Sigma_g is 0
# # sigma_x = cell.norm.G * sqrt((sigma_cells / cells)^2) = cell.norm.g * sigma_cells/cells
# Vmax_and_G <- mutate(Vmax_and_G, G.per.cell = G.per.hr / cell.ml,
#                      G.per.cell.err = G.per.cell * cell.ml.sd / cell.ml)
# 
# # Figure 2a: respiration rate per cell
# # Keeping this commented because I havent' thought enough about how to make units equivalent: oxidation rate is uM/L (porewater or sediment volume?); cells are (I guess per ml )
# # ggplot(Vmax_and_G, aes(x=quant.depth, y = G.per.cell)) + 
# #   geom_point() + 
# #   geom_line() + 
# #   geom_errorbar(aes(ymin = G.per.cell - G.per.cell.err, ymax = G.per.cell + G.per.cell.err)) + 
# #   scale_x_reverse() + 
# #   scale_y_log10() +
# #   coord_flip()


## How much did G decrease from 4.5 to 82.5?
G_preds$G.rel <- G_preds$G.per.hr / G_preds$G.per.hr[1] # For whatever reason G_preds$G[G_preds$depth==1.5] doesn't work
print(paste("G decreased by", (1-(G_preds$G.rel[G_preds$depth == 82.5] / G_preds$G.rel[G_preds$depth == 4.5]))*100, "% from 4.5 to 82.5 cm"))

# update the graphical theme to work properly 
base_theme <- theme_get()
new_theme <- base_theme +
  theme(plot.margin=grid::unit(c(2,2,2,2), "mm"))
theme_set(new_theme)

if(print.plots) {
  print(grid.arrange(p_Vmax.rel, p_Vmax_cell_norm, p_GRR, p_Vmax_rel_G, nrow=1))
}

if(save.plots) {
  tiff("manuscript/ms_plots/2018_new/Fig_2_Vmax_cells_G_vs_depth.tiff", height = 2.5, width = 6.875, units= "in", res = 300, compression = "lzw", type = "cairo")
  cowplot::plot_grid(p_Vmax.rel, p_Vmax_cell_norm, p_GRR, p_Vmax_rel_G,
                     labels = "auto", label_fontface = "plain", 
                     label_x = c(0.22, 0.12, 0.12, 0.18), 
                     label_y = c(0.99),
                     nrow = 1,  align = 'h', axis = 'b')
  dev.off()
  # tiff("manuscript/ms_plots/2018_new/Fig_2_Vmax_cells_G_vs_depth.tiff", height = 2.5, width = 6.875, units = "in", res = 300, compression = "lzw", type = "cairo")
  # grid.arrange(p_Vmax.rel, p_Vmax_cell_norm, p_GRR, p_Vmax_rel_G, nrow=1)
  # dev.off()
}
# Try this the cowplot way
# Try this the cowplot way
theme_set(base_theme)

