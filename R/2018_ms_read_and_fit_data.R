# Step 1 of data analysis for WOR peptidases paper:
# Read in data, calculate v0 values, estimate Michaelis Menten parameters, make plots


###########
# Read, plot raw data
###########

# Manually add cell count data
#   (from ~/Dropbox/WORpeptidases/WOR5.13 enzyme-sample,cellcounts.xls)
depth <- 1:5
cell.ml <- c(4.48e8, 3.50e8, 2.52e8, 1.53e8, 7.37e7)
cell.ml.sd <- c(1.21e8, 8.22e7, 1.03e8, 5.58e7, 5.06e7)
cell_counts <- data.frame(depth, cell.ml, cell.ml.sd)

# Read in the data; drop bad columns that Bill Gates likes to add; add quantitative depth 
data <- read_csv("data/2013-06-03 WOR5.13 Peptides Activity.csv") %>%
  subset(select=-c(X8, X9, X10)) %>%
  left_join(quant_depth, by = "depth")

# Relabel the substrates
new_substrate_labels <- c("1: arg-AMC" = "Arg-AMC", 
                          "2: gly-AMC" = "Gly-AMC", 
                          "3: leu-AMC" = "Leu-AMC", 
                          "4: gly-gly-arg-AMC" = "GGR-AMC", 
                          "5: ala-ala-phe-AMC" = "AAF-AMC", 
                          "6: boc-val-pro-arg-AMC" = "BocVPR-AMC")
data$substrate <- new_substrate_labels[data$substrate]

# Calculate elapsed time
data <- data %>%
  group_by(depth, substrate, conc) %>%
  mutate(elapsed = (time - min(time)) * 24)

# Plot the raw data
p_raw <- ggplot(data, aes(x=elapsed, y=fl, colour=as.factor(conc))) + 
  geom_point() + geom_line() +
  facet_grid(depth~substrate, scales="free")
if(print.deep.cuts) {
  print(p_raw)
}

#############
# Scrub bad data points
#############

# The second timepoint for depth 1, subs 5, conc 25 looks bad
# It appears to be offset by exactly 10, somehow. Will just replace it with the average
#    Which is 3.514 +/- 0.00655
# Calculate average elapsed time for AAF-AMC, subsgtrate 1
av_elapsed <- data %>%
  ungroup() %>%
  filter(substrate == "AAF-AMC" &
           depth == 1 &
           elapsed > 0 &
           elapsed < 10) %>%
  summarise(mean.elapsed = mean(elapsed, na.rm=TRUE)) %>%
  as_vector()


data[which(data$substrate=="AAF-AMC" &
             data$depth==1 &
             data$elapsed > 10), "elapsed"] <- av_elapsed

###
# Calculate slopes, and do a bunch of other stuff, using data frame columns
###

data <- data %>%
  group_by(depth, substrate, conc) %>%
  nest() 

# Make a linear model; catch errors 
# In principle possibly() does this but I can't get it to work right


# Pull our slope. Could use broom; but I don't feel like it
get_slope <- function(m) {
  slope <- tryCatch(
    slope <- coefficients(m)[2],
    error = function(x) NA,
    finally = {}
  )
}

# Do the work of making models and making a column for the slopes
data <- data %>%
  mutate(models = map(data, data_lm),
         slopes.raw = map_dbl(models, get_slope))

######
# Calibrate
######

# Function to calculate slopes of calibration curves
calib_lm <- function(x) {
  m <- tryCatch(
    m <- lm(fl ~ conc.AMC, data = x),
    error = function(x) NA,
    finally = {}
  )
  m
}

# DEpth lookup table
depth.codes <- c("Depth 1" = "1", "Depth 2" = "2", "Depth 3" = "3", "Depth 4" = "4", "Depth 5" = "5")

# Read calibration data, calculate slopes
calib <- read_csv("data/2013-06-13 WOR5.13 Calibration Curve Data.csv") %>%
  dplyr::rename(conc.AMC = `Final [AMC]`, fl = Fl) %>%
  #rename(`Final [AMC]` = conc.AMC, Fl = fl)
  filter(Sample != "ASW") %>%
  group_by(Sample) %>%
  nest() %>%
  mutate(models = map(data, calib_lm),
         slopes.raw = map_dbl(models, get_slope), # AMC conc is in units of uM (0 - 0.45), so 
         depth = as.character(depth.codes[Sample]))

# Find mean calibration slope for killed controls - this gets used in the killed control script
mean.calib.slope.raw <- mean(calib$slopes.raw)

# join calibration curves into sample data
data_calibrated <- data %>%
  mutate(depth = as.character(depth)) %>%
  left_join(calib, by = "depth", suffix = c("", ".calib.curve")) %>%
  mutate(v0.unnorm = slopes.raw / slopes.raw.calib.curve) #"This is in units of micromol per liter per hour. Down the chain I divide by sed concentration (g / liter) to get final units  of micromol per g sediment per hour.")

# Merge in the sample masses; normalize slopes to sample mass
sample_masses <- read_csv("data/2013-06-03 WOR5.13 peptidase activity sample mass (1).csv") %>%
  mutate(depth = as.character(depth))
sample_masses$substrate <- new_substrate_labels[sample_masses$substrate]

data_calibrated <- data_calibrated %>%
  left_join(sample_masses, by = c("depth", "substrate", "conc")) %>%
  rename( mass.sample.g = `mass sample g`) %>%
  mutate(sample.conc = mass.sample.g / 0.004) %>% # all samples were filled with 0.004 g sed
  #mutate(mass.norm.slope.uncal = slopes.raw / (`mass filled` - `mass bottle`))
  mutate(v0.mass.norm = v0.unnorm / sample.conc)


##############
# Remove obviously bad slopes
##############

p_dirty <- ggplot(data_calibrated, aes(x=conc, y=v0.mass.norm))  + 
  geom_point() + 
  geom_line() + 
  facet_grid(depth ~ substrate)
if(print.deep.cuts) {
  print(p_dirty)
}

data_clean <- data_calibrated %>%
  filter(!(depth == 1 & substrate == "BocVPR-AMC" & conc == 75)) %>%
  filter(!(depth == 2 & substrate == "Leu-AMC" & conc == 0)) %>%
  filter(!(depth == 3 & substrate == "Arg-AMC" & conc == 25)) %>%
  filter(!(depth == 3 & substrate == "GGR-AMC" & conc==200)) %>%
  filter(!(depth == 4 & substrate == "AAF-AMC" & conc == 75))

p_clean <- ggplot(data_clean, aes(x=conc, y=v0.mass.norm))  + 
  geom_point() + 
  geom_line() + 
  facet_grid(depth ~ substrate, scales = "free")
if(print.plots) {
  print(p_clean)
}



############
# Create Michaelis-Menten fits
############

###Determine Vmax, Km


# Define Mechaelis-Menten formula
mmForm <- formula(I(v0.mass.norm ~ (Vmax * conc)/(Km + conc)))

# Function to perform 
safe_nls <- function(x, KmGuess=200) {
  # Create a unique starting guess for Vmax
  VmaxGuess <- max(x$v0.mass.norm)
  
  # Safely fit the plots by NLS
  mmFit <- tryCatch(
    nls(formula=mmForm, data=x, start=list(Vmax=VmaxGuess, Km=KmGuess)),
    error=function(err) NA)
  
  mmFit
}

# Function to pull out the average v0 for the highest [S] in any dataset

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) { # from 
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
get_v0_at_max_S <- function(df) {
  max.conc.indices <- df$conc == max(df$conc, na.rm = TRUE) | df$conc >= 400 # because I "jittered" the 400 concentrations so I could calculate separate slopes for each rep at 400 uM
  mean.fake.Vmax <- mean(df$v0.mass.norm[max.conc.indices], na.rm = TRUE)
  mean.fake.Vmax
}

# Calculate MM parameters and predictions
kinetics <- data_clean %>%
  ungroup %>%
  select(depth, substrate, conc, data, models, v0.mass.norm) %>%
  group_by(depth, substrate) %>%
  nest(.key = v0.data) %>%
  mutate(mm.model = map(v0.data, safe_nls))


# There are a few samples/depths for which kinetic parameters can be found, but they're not obvious. Need to rewrite the section below to make them work
AAF1 <- kinetics %>%
  filter(substrate == "AAF-AMC" & depth == 1) 
AAF1_mod <- AAF1$v0.data[[1]] %>%
  nls(formula = mmForm, data = ., start = list(Km = 1027, Vmax = 3)) 
kinetics[kinetics$substrate == "AAF-AMC" & kinetics$depth == 1, ]$mm.model[[1]] <- AAF1_mod


AAF2 <- kinetics %>%
  filter(substrate == "AAF-AMC" & depth == 2)
AAF2_mod <- AAF2$v0.data[[1]] %>%
  nls(formula = mmForm, data = ., start = list(Km = 88, Vmax = 0.51))
kinetics[kinetics$substrate == "AAF-AMC" & kinetics$depth == 2, ]$mm.model[[1]] <- AAF2_mod

Arg1 <- kinetics %>%           
  filter(substrate == "Arg-AMC" & depth == 1)
Arg1_mod <- Arg1$v0.data[[1]] %>%
  nls(formula = mmForm, data = ., start = list(Km = 220, Vmax = 0.26))
kinetics[kinetics$substrate == "Arg-AMC" & kinetics$depth == 1, ]$mm.model[[1]] <- Arg1_mod


# from the nls, extract predictions and
kinetics <- kinetics %>%
  mutate(
    preds.df = map(mm.model, predict_mm),
    Km = map_dbl(mm.model, safe_coef, "Km", "Estimate"),
    Vmax = map_dbl(mm.model, safe_coef, "Vmax", "Estimate"),
    Km.err = map_dbl(mm.model, safe_coef, "Km", "Std. Error"),
    Vmax.err = map_dbl(mm.model, safe_coef, "Vmax", "Std. Error")
  ) #%>%

# Substitute v0 at max S for the one case here Vmax wasn't calculated
leu_2_v0_400 <- kinetics %>%
  filter(depth == 2 & substrate == "Leu-AMC") %>%
  select(v0.data) %>%
  unnest() %>%
  filter(conc >= 400) 
leu.2.fake.Vmax <- mean(leu_2_v0_400$v0.mass.norm)
kinetics[kinetics$depth == 2 & kinetics$substrate == "Leu-AMC", "Vmax"] <- leu.2.fake.Vmax

kinetics <- kinetics %>%
  # Fill in v0 values for when kinetic parameters are unavailable
  mutate(Vmax.complete = Vmax) %>%
  mutate_cond(is.na(Vmax.complete),
              Vmax.complete = map_dbl(v0.data, get_v0_at_max_S)) %>%
  # merge in true depths in centimeters
  mutate(depth = as.numeric(depth)) %>% # used depth as a factor above for convenience, because it acts as categorical
  left_join(quant_depth, by = "depth")




## Make plot with predicted MM values
all_kinetic_preds <- kinetics %>%
  unnest(preds.df)

v0_data_long <- kinetics %>%
  unnest(v0.data)

# Calculate the killed controls so I can add them to the plot
source("manuscript/ms_R/2018_killed_controls_rewritten.R") # data frame of killed "activities" is ks


p_all_kinetics <- ggplot(v0_data_long, aes(x=conc, y=v0.mass.norm * 1000))  + # note this is now in nmol g-1 sed hr-1
  geom_line(data = all_kinetic_preds, aes(x=s, y=preds * 1000)) + 
  geom_point() + 
  geom_point(data = ks, aes(x=conc, y=v0.mass.norm * 1000), shape = 2) +
  ylab(expression(paste(v[o], ", ", "nmol ", g^{-1}, " sed ", hr^{-1}))) + 
  xlab(expression(paste("substrate concentration, ", mu, "M"))) + 
  facet_grid(quant.depth ~ substrate, scales = "free") # somehow quant.depth has not been merged in, I think
if(print.plots) {
  print(p_all_kinetics)
}
if(save.plots) {
  ggsave("manuscript/ms_plots/2018_new/Fig1.tiff", height = 4, width = 7.5, units = "in", dpi=300, compression = "lzw", type="cairo")
}

# Create summarise of parameters
kin_summ <- kinetics %>%
  summarise(
    mean.Vmax = mean(Vmax, na.rm = TRUE),
    median.Vmax = median(Vmax, na.rm = TRUE),
    min.Vmax = min(Vmax, na.rm = TRUE),
    IQR.lo.Vmax = quantile(Vmax, 0.25, na.rm = TRUE),
    IQR.hi.Vmax = quantile(Vmax, 0.75, na.rm = TRUE),
    max.Vmax = max(Vmax, na.rm = TRUE),
    n.Vmax = sum(!is.na(Vmax)),
    n.Km = sum(!is.na(Km)),
    
    
    mean.Km = mean(Km, na.rm = TRUE),
    median.Km = median(Km, na.rm = TRUE),
    min.Km = min(Km, na.rm = TRUE),
    IQR.lo.Km = quantile(Km, 0.25, na.rm = TRUE),
    IQR.hi.Km = quantile(Km, 0.75, na.rm = TRUE),
    max.Km = max(Km, na.rm = TRUE)
  ) %>%
  gather(key = "quantity", value = "value", mean.Vmax:max.Km)
if(print.plots) {
  message("Summary of kinetic parameters of the single-cuvette enzyme data:")
  print(kin_summ)
}

# Check which substrates were fastest and slowest. 
# Order substrates in order of decreasing median

Vmax_mod <- aov(log10(Vmax) ~ substrate, data = kinetics)
# plot(Vmax_mod) # could be worse. I'll use it.
summary(Vmax_mod)
Vmax_tukey <- TukeyHSD(Vmax_mod)
print(Vmax_tukey)


p_all_Vmax <- ggplot(kinetics, aes(x=fct_reorder(substrate, Vmax, .fun = median, .desc = TRUE), y=log10(Vmax))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_pointrange(aes(ymin = log10(Vmax - Vmax.err), ymax = log10(Vmax + Vmax.err)), position=position_jitter(width = 0.1, height = 0), size = 0.3) + 
  #ylab(expression(paste(log_[10], "(", V_[max], ")", ", nmol ", g^[-1], " sed ", hr^[-1])) 
  xlab("substrate") + 
  ylab(expression(paste(log[10], "(", V[max], "), nmol ", g^{-1}, " sed ", hr^{-1}))) + 
  annotate("text", x=1:6, y=0.7, label= c("A",  "A", "A", "A,B", "A,B", "B"))
if(print.deep.cuts) {
  print(p_all_Vmax)
}
# Remember, I save this as a 2-panel plot below

Vmax_mod <- aov(log10(Vmax) ~ substrate, data = kinetics)
#plot(Vmax_mod) # could be worse. I'll use it.
summary(Vmax_mod)
Vmax_tukey <- TukeyHSD(Vmax_mod)
print(Vmax_tukey)


Km_mod <- aov(rank(Km) ~ fct_reorder(substrate, Km, mean, .desc = TRUE), data = kinetics)
#plot(Km_mod) # also mediocre but not truly terrible
summary(Km_mod) # Highly significant when expressed as ranks
Km_tukey <- TukeyHSD(Km_mod)
print(Km_tukey)
#Km_tukey <- agricolae::HSD.test(Km_mod, "substrate", group = TRU#plot(Km_tukey)

Km_rank <- kinetics %>%
  mutate(Km.rank = rank(Km))
p_all_Km <- ggplot(Km_rank, aes(x=fct_reorder(substrate, Km.rank, mean, .desc = TRUE), y=Km.rank)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 2) + 
  xlab("substrate") + 
  ylab(expression(paste("rank of ", K[m]))) + 
  annotate("text", x = 1:6, y=35, label = c("A", "A,B", "A,B", "B", "B", "B"))
if(print.deep.cuts) {
  print(p_all_Km)
}
if(save.plots) {
  tiff("manuscript/ms_plots/2018_new/supplemental/Fig_s1_all_Km_and_Vmax.tiff", height = 6, width = 4, units = "in", res = 300, type = "cairo", compression = "lzw")
  gridExtra::grid.arrange(p_all_Vmax, p_all_Km, nrow = 2)
    #ggsave("manuscript/ms_plots/resubmit/supplemental_Vmax.png", height = 3, width = 4, units = "in", dpi = 300)
  dev.off()
  
}


# plot(lm(log10(Vmax) ~ substrate, data=kinetics)) # residuals are badly non-normal with or without log10 transformation
# Vmax.dunn.test <- dunn.test::dunn.test(kinetics$Vmax, as.factor(kinetics$substrate), method = "bh") # sig diff's found, but this is the wrong test - need some kind of repeated measures here
# 
# 
# Vmax_rank_summ <- Vmax_rank %>%
#   group_by(substrate) %>%
#   summarise(mean.rank = mean(Vmax.rank))
# 
# # Check which substrates had highest and lowest Km
# max.Km <- max(kinetics$Km, na.rm = TRUE) # 1312
# kinetics_for_rank <- kinetics %>%
#   mutate_cond(is.na(Km), Km = max.Km) %>%
#   group_by(depth) %>%
#   mutate(Km.rank = rank(Km))
# Km_rank_summ <- kinetics_for_rank %>%
#   group_by(substrate) %>%
#   summarise(mean.Km.rank = mean(Km.rank))
# # ggplot(kinetics_for_rank, aes(x=substrate, y=Km.rank, colour=as.factor(depth))) +
# #   geom_point(position = position_jitter(width = 0.3, height = 0)) + 
# #   geom_point(data=Km_rank_summ, aes(x=substrate, y=mean.Km.rank), colour = "black", size = 4)
# Km.dunn.test <- dunn.test::dunn.test(kinetics$Km, as.factor(kinetics$substrate), method = "bh") # no significant differences in Km among substrates

# Get rid of obsolete variables
rm(AAF1_mod, AAF2_mod, Arg1_mod, av_elapsed, cell.ml, cell.ml.sd, depth, depth.codes, Km.dunn.test, leu.2.v0.400, max.Km, new_substrate_labels, 
   p_all_kinetics, p_dirty, p_raw, Vmax.dunn.test, AAF1, AAF2, Arg1, calib, data, data_cal, data_clean, kinetics_for_rank, Km_rank_summ, 
   quant_depth, Vmax_rank, Vmax_rank_summ, leu_2_v0_400)
message("Drew: look back at why there are a bunch of variables in this list that no longer exist")
