detachAllPackages <- function() {    
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

detachAllPackages()

library(tidyverse)


rk_phyla <- read_csv("manuscript/ms_data/Baker_16S_to_use.csv") %>%
  gather(key = "depth", value = "reads", 2:5) %>%
  mutate(zone = c("8-12cm" = "SRZ", "36_30" = "SMTZ", "30-32cm" = "SMTZ", "52-54cm" = "MRZ")[depth]) %>%
  filter(!is.na(phylum)) %>%
  group_by(depth) %>% # Note that 'depth' distinguishes between the two SRZ depths
  mutate(rel.reads = reads / sum(reads, na.rm = TRUE)) %>%
  top_n(10, rel.reads) %>%
  ungroup() #%>%
  #mutate(zone = factor(zone, levels = rev(levels(factor(zone))), ordered = TRUE)) 
 


SRZ_reads <- rk_phyla %>% 
  filter(zone == "SRZ") %>%
  arrange(rel.reads)

phyla.in.my.order <- rev(c("Proteobacteria", "Chloroflexi", "Bacteroidetes",
                       "Planctomycetes", "Euryarchaeota", "Bathyarchaeota",
                       "Latescibacteria", "Acidobacteria", "Gemmatimonadetes",
                        "Aminicenantes", "Aenigmarchaeota", "Lokiarchaeota",
                       "Hadesarchaea", "All other bacteria", "All other archaea"))



phyla_to_print <- rk_phyla %>%
  #mutate(ordered.taxon = factor(Taxon, levels = rev(c(SRZ.phyla, "Aenigmarchaeota", "ANME-1", "Hadesarchaea", "all other archaea")), ordered = TRUE)) %>%
  mutate(ordered.phylum = factor(phylum, levels = phyla.in.my.order, ordered = TRUE)) %>%
  group_by(zone, ordered.phylum) %>%
  mutate(n = length(rel.reads),
         mean.rel.reads = mean(rel.reads, na.rm = TRUE),
         min.range.rel.reads = min(range(rel.reads, na.rm = TRUE)),
         max.range.rel.reads = max(range(rel.reads, na.rm = TRUE)))
  

p_phyla <- ggplot(phyla_to_print, aes(x=ordered.phylum, y=mean.rel.reads, colour = zone, shape = zone)) + 
  geom_pointrange(aes(ymin = min.range.rel.reads, ymax = max.range.rel.reads)) + 
  #geom_point() + 
  scale_colour_brewer(name = "redox zone", palette = "Dark2") + 
  scale_shape_discrete(name = "redox zone") + 
  scale_y_continuous(name = "relative abundance", labels = scales::percent) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.title.y = element_blank(),
        legend.position = c(0.83, 0.23),
        legend.background = element_rect(colour = "darkgray"))

print(p_phyla)
detachAllPackages()

