# Producing a figure showing salinity overtime
# I will be using the 100 simulation dataframe

# Setting wd
rm(list=ls())
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/")

# Packages
library(tidyverse)
library(patchwork)

# Importing the 100 simulation df.
df.100 <- read.csv("sim_data_100")
nauplii.100 <- read.csv("nauplii_floquet_100")
cop.100 <- read.csv("cop_flo_100")
chali.100 <- read.csv("chalimi_flo_100")
adult.100 <- read.csv("adult_flo_100")

# I want to plo the abundnace of the sea lice, but only keep the salinity from one sumulation at three seperate pho values

# Pulling a random simulation and manpluating the df
rand.n <- sample(seq(0,100),1)
df.sub.n <- filter(df.100, df.100$Simulation.Number == rand.n &
                     Pho.Value==c(1.8217, 0.6739))
df.sub.n$Pho.Value <- as.factor(df.sub.n$Pho.Value)
colnames(df.sub.n)[1] <- ("Pho Value")
colnames(df.sub.n)[8] <- "Adult Females"

# Building the plots
sal.p <- ggplot(df.sub.n, aes(Time, Salinity)) +
                geom_point(aes(color = `Pho Value`, shape=`Pho Value`)) +
                theme_classic() +
                scale_color_manual(values = c("steelblue", "red")) +
                theme(legend.position = "none")

abund.p <- ggplot(df.sub.n, aes(Time, log(`Adult Females`))) +
                geom_point(aes(color = `Pho Value`)) +
                scale_color_manual(values = c("steelblue", "red")) +
                theme_classic()



# Combing the plot
sal.p + abund.p


# Doesn't look to bad. Still a few edits that need to get done.
 # - I'm not the biggest fan of the color. That might need to change.
 # - The order of the pho values might also need to changed. Have the largest one at the top
 # - The ratio plot also needs to be added to this. More to come.

# Trying to find the max value

f %>% group_by(Gene) %>% summarise(Value = max(Value))
df.sub.n %>%
  group_by(Simulation.Number) %>%
  summarise(Value=max(Nauplii))

lapply(split(df.sub.n, df.sub.n$), function(y) max(y$Value))
