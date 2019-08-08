# Make plot of Km vs depth

set.seed(0)
p_Km_v_depth <- ggplot(kinetics, aes(x=quant.depth, y=Km, colour = substrate)) + 
  geom_smooth(aes(group = 1), method="lm", se=TRUE, colour="black", size=3) +
  geom_pointrange(aes(ymin=Km-Km.err, ymax=Km+Km.err), position = position_jitter(width = 0.5)) + 
  #geom_line() + 
  geom_smooth(method="lm", se=FALSE) + 
  scale_y_log10(name = expression(paste(K[m], ", ", mu, "M")),
                limits=c(10, 1200)
                ) + 
  scale_x_reverse(name = "depth, cmbsf") + 
  scale_color_brewer(palette = "Dark2") + 
  #annotation_logticks(sides = "l") + 
  coord_flip() 
if(print.plots) {
  print(p_Km_v_depth)
}
if(save.plots) {
  ggsave("manuscript/ms_plots/2018_new/Fig_4_Km_downcore.tiff", height = 3, width = 6.875/2, units = "in", dpi = 900, compression = "lzw", type = "cairo")
}

