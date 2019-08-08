#######
# Calculate killed controls
#######

# library(plyr)
# library(lubridate)
# library(ggplot2)
# library(lmstats)

##read in the dataframe for our killed controls
#killed<-read.csv("data/2013-06-14 WOR5.13 Killed Controls.csv")
killed <- read_csv("data/2013-06-14 WOR5.13 Killed Controls.csv")
killed<-killed[-25,]

# Elapsed time is in days - convert to hours
# killed <- ddply(killed, c("Sample", "Conc"), mutate,
#                 elapsed = (Time-min(Time))*24)
killed <- killed %>%
  group_by(Sample, Conc) %>%
  mutate(time = (Time - min(Time))* 24) %>% # time is elapsed time in hours; Richard recorded it in days
  rename(fl = Fl) # rename columns for consistency with kinetics

# Pot the raw data. Only two data points per substrate, annoyingly
p_raw <- ggplot(killed, aes(x=time, y=fl, colour=as.factor(Conc))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~Sample)
if(print.plots) {
  print(p_raw)
}
# print(p_raw)
# Relabel the substrates & rename columns
killed_substrate_labels <- c("ASW" = "blank", 
                          "Sub 1" = "Arg-AMC", 
                          "Sub 3" = "Leu-AMC", 
                          "Sub 4" = "GGR-AMC", 
                          "Sub 5" = "AAF-AMC", 
                          "Sub 6" = "BocVPR-AMC")
killed$substrate <- killed_substrate_labels[killed$Sample]

# Calculate slopes
killed_slopes <- killed %>%
  group_by(substrate, Conc) %>%
  nest(time, fl, .key = v0.data) %>%
  mutate(models = map(v0.data, data_lm), # data_lm and get_slope as separate .R files
        slopes.raw = map_dbl(models, get_slope))


######
# Express slopes as per mass and calibrate
######
# Richard did a ton of work to ensure that sample masses were exactly 0.5 g
killed_slopes$sample.conc <- 0.5 / 0.004 # Remember, that the volume was 4 ml
attr(killed_slopes$sample.conc, "units") <- "g L^-1"

# Calibrate
    # This assumes that the WOR peptidase ms script has been run through line 197
#calib.slope.mean <- mean(calib_curves$calib.slope)


killed_slopes$v0.mass.norm <- killed_slopes$slopes.raw / mean.calib.slope.raw / killed_slopes$sample.conc
attr(killed_slopes$v0.mass.norm, "units") <- "mu mol g-1 sed hr-1"

# artificially repeat data frame so that it can be faceted by depth
ks15 <- killed_slopes
ks15$quantDepth <- 1.5
ks45 <- killed_slopes
ks45$quantDepth <- 4.5
ks285 <- killed_slopes
ks285$quantDepth <- 28.5
ks585 <- killed_slopes
ks585$quantDepth <-58.5
ks825 <- killed_slopes
ks825$quantDepth <-82.5
# Put it into a list in order to convert to a data frame using do.call, rename for congruency with 
ksl <- list(ks15, ks45, ks285, ks585, ks825)
ks <- do.call(rbind, ksl) %>%
  rename(conc = Conc) %>%
  filter(substrate != "blank") %>%
  mutate(substrate = as.factor(substrate))

rm(killed, p_raw, killed_substrate_labels, killed_slopes, ks15, ks45, ks285, ks585, ks825, ksl)

