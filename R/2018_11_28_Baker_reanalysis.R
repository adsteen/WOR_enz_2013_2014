########
# Reprocess Brett's data
########

library(readr)
library(plyr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

# read the data & group the data
SRZ <- read_csv("data/Baker_WOR_data_ADS_SRZ.csv") # 22 cols
SRZ$depth <- "SRZ"
SMTZ <- read_csv("data/Baker_WOR_data_ADS_SMTZ.csv")[ , -8]
SMTZ$depth <- "SMTZ"
MRZ <- read_csv("data/Baker_WOR_data_ADS_MRZ.csv") # 22 cols
MRZ$depth <- "MRZ"

ee <- do.call(rbind, list(SRZ=SRZ, SMTZ=SMTZ, MRZ=MRZ))

# How many NA lineages?
sum(is.na(ee$Lineage)) # There's 298 NAs
(1647-1616) + (4665 - 4484) + (2603-2515) #= 300, which is close enough to 298 for my likes
# SO WHERE ARE THE MCGs AND WHATNOT

# Determine which of these are peptidases

# If you have "peptidase" in your name, OR an EC beginning with 3.4, you MAY be a peptidase <\Foxworthy>
ee$peptidase.in.EC <- grepl("EC:3.4", ee$EC_Number)
ee$peptidase.in.name <- grepl("eptidase", ee$Product)

# Strike out non-peptidases
ee_pep <- ee[ee$peptidase.in.EC | ee$peptidase.in.name, ]

# Strike out specifically identified non-peptidases
ee_pep <- subset(ee_pep, Product !="Bacterial protein of unknown function (DUF896)") # lost ~10
ee_pep <- subset(ee_pep, Product !="Chagasin family peptidase inhibitor I42") # lost 56
ee_pep <- subset(ee_pep, Product != "Gag polyprotein, inner coat protein p12") # lost 1

# EC 3.4.11.18 has NAs for product names for some reason - rectify this
ee_pep[ee_pep$EC_Number == "EC:3.4.11.18" & is.na(ee_pep$Product), "Product"] <- "Methionine aminopeptidase"


# # Check this test - how well does it perform? (Note: 8912 rows total)
# sum(ee$peptidase.in.name + ee$peptidase.in.EC == 2) #943
# sum(ee$peptidase.in.name + ee$peptidase.in.EC == 1) #3075
# 
# ee_either <- ee[ee$peptidase.in.name | ee$peptidase.in.EC, ]
# print(ee_either[1:300 , c("Product", "EC_Number", "peptidase.in.name", "peptidase.in.EC")])

# Write output which can be copied into Brett's Excel workbook
# write_csv(ee_pep, "data/genomic_data/Baker_WOR_data_peptidase_only.csv")

### Break out taxonomic IDs

# ee.taxa <- ee$Lineage
# taxa_list <- strsplit(ee.taxa, split=";")

# Split out the various categories of lineage each into their own category
pep_tax <- ddply(ee_pep, c("Lineage"), mutate, 
      domain = str_split(Lineage, ";")[[1]][1],
      phylum = str_split(Lineage, ";")[[1]][2],
      class = str_split(Lineage, ";")[[1]][3],
      order = str_split(Lineage, ";")[[1]][4],
      family = str_split(Lineage, ";")[[1]][5],
      genus = str_split(Lineage, ";")[[1]][6],
      species = str_split(Lineage, ";")[[1]][7],
      strain = str_split(Lineage, ";")[[1]][8]
      )

# Drop out unassigned lineages
pep_tax <- pep_tax[!is.na(pep_tax$Lineage), ]

# Clean up some bad domain assignments - Lokiarchaeota, MBG-D and MCG are assigned as domains
pep_tax$phylum[pep_tax$domain == "Lokiarchaeota"] <- "Lokiarchaeota"
pep_tax$domain[pep_tax$domain == "Lokiarchaeota"] <- "Archaea"
pep_tax$phylum[pep_tax$domain == "MBG-D"] <- "MBG-D"
pep_tax$domain[pep_tax$domain == "MBG-D"] <- "Archaea"
pep_tax$phylum[pep_tax$domain == "MCG-1"] <- "MCG"
pep_tax$domain[pep_tax$domain == "MCG-1"] <- "Archaea"
pep_tax$phylum[pep_tax$domain == "MCG-6"] <- "MCG"
pep_tax$domain[pep_tax$domain == "MCG-6"] <- "Archaea"

# Rename MBG-D as Thorarchaeota
pep_tax$phylum[pep_tax$phylum == "MCG"] <- "Thorarchaeota"

if(save.plots) {
  write_csv(pep_tax, "plots/raw_annotations.csv")
}

# Count the occurrances of each peptidase in each depth zone
pep_summ <- ddply(pep_tax, c("depth", "Product"), summarise, 
                  count=length(species))

# Determine frequency
pep_summ <- ddply(pep_summ, c("depth"), mutate,
                  count.freq = count/sum(count))

# Order by frequency in SRZ (I guess)
   # order in SRZ
SRZ.ord <- order(subset(pep_summ, depth=="SRZ")$count.freq, decreasing = TRUE) # NOTE YOU WILL RUN INTO PROBLEMS IF THERE ARE REPRESENTATIVES THAT AREN'T IN SRZ
SRZ.peps <- pep_summ$Product[pep_summ$depth=="SRZ"][SRZ.ord]
# add the other products, in no particular order
all.peps <- c(SRZ.peps, unique(pep_summ$Product)[!(unique(pep_summ$Product)) %in% unique(SRZ.peps)])
# Order by order of all.peps
pep_summ$Product <- factor(pep_summ$Product, levels = all.peps, ordered = TRUE)



# Note: Could do ddply(pep_tax, c("depth.zone"), count, vars=pep_tax) or something similar
# Check nnumeric contribtuion of most abundant peptidases
pep_summ_print <- pep_summ[with(pep_summ, order(count.freq, decreasing = TRUE)), ]
print(pep_summ_print)

# Calculate summed contribution of peptidases at each depth
ddply(pep_summ[pep_summ$Product == "Peptidase family C25" |
                 pep_summ$Product == "Methionine aminopeptidase" | 
                 pep_summ$Product == "Zinc carboxypeptidase", ],
      c("depth"), function(x) sum(x$count.freq, na.rm = TRUE))

# Reorder depths
pep_summ_1a <- mutate(pep_summ, depth = factor(depth, levels=c("SRZ", "SMTZ", "MRZ" ), ordered=TRUE))


p_pep_ID <- ggplot(pep_summ_1a, aes(x=Product, y=count.freq, colour=depth, shape=depth)) + 
  geom_point() + 
  #scale_colour_manual(values=wes_palette("Zissou", 3)) +
  scale_colour_brewer(palette = "Dark2", guide = FALSE) +
  scale_shape_discrete(guide = FALSE) + 
  scale_y_log10() +
  annotation_logticks(sides = "l") + 
  ylab("frequency") + 
  theme(axis.text.x = element_text(angle=-60, hjust=0, size=6),
        axis.title.x = element_blank())
if(print.plots) {
  suppressWarnings(print(p_pep_ID))
}
if(save.plots) {
  ggsave("plots/Fig_5_peptidase_identities.tiff", p_pep_ID, height=3.5, width=3.5, units="in", dpi=300, compression="lzw", type="cairo")
}
# 

vp6a <- grid::viewport(x = 0.5/2, y = 0.5,  height = 1, width = 0.5)
vp6b <- grid::viewport(x = 0.6 + 0.4/2, y = 0.5, height = 1, width = 0.4)
if(print.plots) {
  print(p_pep_ID, vp = vp6a)
  print(p_phyla, vp = vp6b)
}

if(save.plots) {
  tiff("plots/Fig_6ab_v3.tiff", height = 3, width = 7.5, units = "in", res = 600, compression = "lzw", type = "cairo")
  print(p_pep_ID, vp = vp6a)
  print(p_phyla, vp = vp6b)
  dev.off()
}


#cowplot::plot_grid(p_Fig1a, p_phyla, nrow = 1, rel_heights = c(2, 1))
# Editor's note: late in the process I got a little sketched out by some of the taxonomic assignments, so I dropped them from the analysis.

# ########
# # Calculate frequency by taxonomy
# ########
# # Calculate by domain
# domain_summ <- ddply(pep_tax, c("depth", "domain"), summarise, count=length(species))
# # drop unclassified 
# domain_summ <- domain_summ[!is.na(domain_summ$domain), ]
# domain_summ <- ddply(domain_summ, c("depth"), mutate,
#                      count.freq = count / sum(count))
# 
# domain_summ$quant.depth <- as.numeric(as.factor(domain_summ$depth))
# # Abundance by domain as a function of depth
# p_Fig1b <- ggplot(subset(domain_summ, domain != "Viruses"), aes(x=quant.depth, y=count.freq, colour=domain, shape=domain))+ 
#   geom_point() + 
#   geom_line(aes(group=domain)) + 
#   # scale_colour_manual(values=wes_palette(name="Zissou", n=4), guide=TRUE) +
#   # scale_colour_manual(values=wes_palette(name="Zissou", n=4)) +
#   scale_colour_brewer(palette="Set1") +
#   scale_shape_discrete(guide=FALSE)+
#   scale_x_continuous(breaks=1:3, labels=c("MRZ", "SMTZ", "SRZ")) +
#   ylab("relative abundance") +
#   coord_flip() + 
#   theme(axis.title.y=element_blank())
# if(print.plots) {
#   print(p_Fig1b)
# }
# # print(p_Fig1b)
# 
# # Create a table of all frequencies by domain
# pep_tax$domain[is.na(pep_tax$domain)] <- "unclassified"
# pep_tax$phylum[is.na(pep_tax$phylum) & pep_tax$domain != "unclassified"] <- 
#   paste("unclassified", pep_tax$domain[is.na(pep_tax$phylum) & pep_tax$domain != "unclassified"])
# pep_tax$phylum[pep_tax$domain == "unclassified"] <- "unclassified"
# 
# phy_summ <- plyr::count(pep_tax, vars = c("depth","domain", "phylum")) 
# 
# # Calculate relative frequency for each depth, separately by domains
# phy_bac_arc <- ddply(subset(phy_summ, domain=="Bacteria" | domain=="Archaea"),
#                      c("domain", "depth"), mutate,
#                      rel.freq = freq/sum(freq))
# 
# other.cutoff <- 0.06
# # If all of the relative frequencies for a phylum are less than 0.05, set hte phylum to be "other"
# phy_bac_arc <- ddply(phy_bac_arc, c("domain", "phylum"), mutate,
#                      is.other = max(rel.freq) < other.cutoff)
# 
# # Condense phyla which contribute less than 6% of peptidases, at any depth, into "other"
# phy_bac_arc$condensed.phylum <- phy_bac_arc$phylum
# phy_bac_arc$condensed.phylum[phy_bac_arc$is.other] <- "other"
# 
# # Convert condensed phylum into a factor, so that we can manipulate its order
# phy_bac_arc$condensed.phylum <- as.factor(phy_bac_arc$condensed.phylum)
# 
# # Arrange the data frame so that the bar groupings show up in the right place
# phy_bac_arc <- arrange(phy_bac_arc, condensed.phylum)
# 
# # Split out bacterial data frame, rearrange levels so that the legend assigns gray to other
# phy_bac <- subset(phy_bac_arc, domain=="Bacteria")
# phy_bac$phylum.ord <- relevel(phy_bac$condensed.phylum, ref="other")
# phy_bac$phylum.ord <- factor(phy_bac$phylum.ord, levels=rev(levels(phy_bac$phylum.ord)), ordered=TRUE)
# 
# phy_bac <- rbind(subset(phy_bac, condensed.phylum != "other"), subset(phy_bac, condensed.phylum == "other"))
# 
# # Split out archaeal data frame, rerrange levels so that the legend assigns gray to other
# phy_arc <- subset(phy_bac_arc, domain=="Archaea")
# phy_arc$phylum.ord <- relevel(phy_arc$condensed.phylum, ref="other")
# phy_arc$phylum.ord <- factor(phy_arc$phylum.ord, levels=rev(levels(phy_arc$phylum.ord)), ordered=TRUE)
# 
# ### Calculate sum & unique number of minor phyla
# # Unique number:
# length(unique(phy_bac[phy_bac$phylum.ord == "other" , "phylum"])) #18
# length(unique(phy_arc[phy_arc$phylum.ord == "other", "phylum"]))
# 
# 
# 
# # Figure out what colors are being used here
# named.colors <- brewer_pal(type="qual")(n=4)
# na.gray <- brewer_pal(type="qual")(n=8)[8]
# Archaea.colors <- c(named.colors, na.gray)
# 
# # Plot separately by domain
# p_arc_tax <- ggplot(phy_arc, aes(x=depth, y=rel.freq, fill=condensed.phylum)) + 
#   geom_bar(stat="identity", colour="black") + 
#   #scale_fill_brewer(type="qual", name="phylum") +
#   scale_fill_manual(name="phylum", values=Archaea.colors) +
#   xlab("depth zone") +
#   ylab("relative frequency") + 
#   coord_flip() + 
#   facet_wrap(~domain)
# if(print.plots) {
#   print(p_arc_tax)
# }
# 
# 
# p_bac_tax <- ggplot(phy_bac, aes(x=depth, y=rel.freq, fill=phylum.ord)) +
#   geom_bar(stat="identity", colour="black") +
#   scale_fill_brewer(type="qual", name="phylum") +
#   xlab("depth zone") +
#   ylab("relative frequency") + 
#   coord_flip() + 
#   facet_wrap(~domain)
# if(print.plots) {
#   print(p_bac_tax)
# }
# 
# 
# 
# #### New fig 1: distribution of genes, and domain depth trends
# vp1.frac <- 0.6
# vp1 <- viewport(x=vp1.frac/2, y=0.5, width=vp1.frac, height=1)
# vp2 <- viewport(x=vp1.frac + (1-vp1.frac)/2, y=0.5, width=1-vp1.frac, height=1)
# # # png("manuscript/ms_plots/gene_dist_and_depth_v2.png", height=3, width=7.5, units="in", res=300)
# # print(p_Fig1a, vp=vp1)
# # print(p_Fig1b, vp=vp2)
# # # dev.off()
# 
# 
# #### Giant Fig 1:
# vp1.top.frac <- 0.6
# vp.top.height <- 0.4
# vp1 <- viewport(x=vp1.frac/2, y=(1-vp.top.height) + vp.top.height/2, width=vp1.frac, height=vp.top.height)
# vp2 <- viewport(x=vp1.frac + (1-vp1.frac)/2, y=(1-vp.top.height) + vp.top.height/2, width=1-vp1.frac, height=vp.top.height)
# vp3 <- viewport(x=0.5, y=(1-vp.top.height)/2 + (1-vp.top.height)/4, width=1, height=(1-vp.top.height)/2)
# vp4 <- viewport(x=0.5, y=(1-vp.top.height)/4, width=1, height=(1-vp.top.height)/2)
# 
# if(save.plots) {
#   tiff(file="manuscript/ms_plots/2018_new/Fig_1_xlib.tiff", height=7, width=7.5, units="in", compression="lzw", type="cairo", res=300)
#   print(p_Fig1a, vp=vp1)
#   print(p_Fig1b, vp=vp2)
#   print(p_arc_tax, vp=vp3)
#   print(p_bac_tax, vp=vp4)
#   dev.off()
# }
# 
# 
# ### Calculate how many bacterial phyla are represented
# ddply(phy_bac_arc, c("domain", "depth"), summarise,
#       num.phyla = length(phylum))
# ddply(phy_bac_arc, c("domain", "depth"), summarise, 
#       frac.unclassified = sum(rel.freq[is.other==TRUE]))
# 
# # Which bacterial classes are most abundant
# class_freq <- ddply(pep_tax, c("domain", "depth"), count, "class")
# class_freq <- ddply(class_freq, c("domain", "depth"), mutate,
#                     rel.freq = freq / sum(freq))
# class_freq <- arrange(class_freq, domain, depth, desc(rel.freq))
# if(save.plots) {
#   write_csv(class_freq, "manuscript/ms_data/gene_freq_by_class_and_depth.csv")
# }
# 
# # Identify contributions by phyla
# 
# bac_class <- subset(pep_tax, domain=="Bacteria")
# bac_class <- count(bac_class, "class")
# bac_class$rel.freq <- bac_class$freq / sum(bac_class$freq)
# bac_class <- arrange(bac_class, rel.freq, decreasing=TRUE)
# if(save.plots) {
#   write_csv(bac_class, "manuscript/ms_data/bac_pep_contributions.csv")
# }
# # write.csv("manuscript/ms_data/bac_")
# 
# 
# ### Count how many domains went into 'unclassified'
# # ddply(phy_bac_arc, c("depth"), summarise, 
# #       total.phyla=sum(is.other))
# 
# # Make supplemental tables
# phy_bac_arc <- arrange(phy_bac_arc, desc(domain), desc(depth), desc(phylum))
# 
# if(save.plots) {
#   library(gridExtra)
#   pdf("manuscript/ms_plots/supp_table_pep_by_depth_SRZ.pdf", height=11, width=8.5)
#   grid.table(subset(pep_summ, depth=="SRZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_pep_by_depth_SMTZ.pdf", height=11, width=8.5)
#   grid.table(subset(pep_summ, depth=="SMTZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_pep_by_depth_MRZ.pdf", height=11, width=8.5)
#   grid.table(subset(pep_summ, depth=="MRZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_phyla_by_depth_bac_SRZ.pdf", height=11, width=8.5)
#   grid.table(subset(phy_bac_arc, domain=="Bacteria" & depth=="SRZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_phyla_by_depth_bac_SMTZ.pdf", height=11, width=8.5)
#   grid.table(subset(phy_bac_arc, domain=="Bacteria" & depth=="SMTZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_phyla_by_depth_bac_MRZ.pdf", height=11, width=8.5)
#   grid.table(subset(phy_bac_arc, domain=="Bacteria" & depth=="MRZ"))
#   dev.off()
#   pdf("manuscript/ms_plots/supp_table_phyla_by_depth_arc.pdf", height=11, width=8.5)
#   grid.table(subset(phy_bac_arc, domain=="Archaea"))
#   dev.off()
#   
# }
# #