# Further look into the generated results 
rm(list=ls())
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/first_chapter/winter_2020/model_code/model_figures/model_figure_data/")

# Packages needed 
library(tidyverse)
library(patchwork)
library(cowplot)

#### Importing data frames ####
# Bay Data 
baydespoir_prof <- read.csv("../../data_files/baydespoir.nosync/baydespoir_prof.csv")

# Simulation.data.frame for each N simulation
df.sim.100 <- read.csv("sim_data_100")
df.sim.50 <- read.csv("sim_data_50")
df.sim.10.1 <- read.csv("sim_data_10.1")
df.sim.10.2 <- read.csv("sim_data_10.2")
  
# Nauplii life stage
nauplii.100 <- read.csv("nauplii_floquet_100")
nauplii.50 <- read.csv("nauplii_floquet_50")
nauplii.10.1 <- read.csv("nauplii_floquet_10.1")
nauplii.10.2 <- read.csv("nauplii_floquet_10.2")
  
# Copepodid life stage
cop.100 <- read.csv("cop_flo_100")
cop.50 <- read.csv("cop_flo_50")
cop.10.1 <- read.csv("cop_flo_10.1")
cop.10.2 <- read.csv("cop_flo_10.2")

# Chalimis life stage
chalimi.100 <- read.csv("chalimi_flo_100")
chalimi.50  <- read.csv("chalimi_flo_50")
chalimi.10.1 <- read.csv("chalimi_flo_10.1")
chalimi.10.2 <- read.csv("chalimi_flo_10.2")

# Adult life stage 
adult.100 <- read.csv("adult_flo_100")
adult.50 <- read.csv("adult_flo_50")
adult.10.1 <- read.csv("adult_flo_10.1")
adult.10.2 <- read.csv("adult_flo_10.2")

#### Possible Useful Distrubution ####
# Asymmetric Laplace Distrubution ----
dx <- 0.5
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

pho.vector <-  rev(seq(0, 1.913, 1.913/10))
sub.m <- 32.25 # Should stay 32.25
pho <- pho.vector # Should follow the vector above
k <- 1.122308 # Should stay 1.122308

sub.asy.laplace <- function(sub.m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
}

sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, pho[1], k),
                                    cumsum(sub.asy.laplace(sub.m, pho[1], k))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
  
# Salinity Value
sample(sub.sal.df$`Salinity Values`, 1) 

# The initinal Asymmetric Laplace distrubution was fitted in previous R script. T
# The full fitting process can be found within the stage_prediction_matrix code  

# Temperature ----
# Time 
dtime = 1
TT <- 2*365 # Two years
time <-  seq(0, TT, dtime)
time0 <-time[1]

# Temperature Function parameter 
a = 6.000949
b1 =  -4.785375
b2 = -2.518553

# Temperature Functions 
Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
  return(Temp)
}


plot(time, Temp(time), type="l", xlab="Time (days)", ylab="Ocean Temperature (Celsius)")

#### Looking into the frequency that the salinity is being pulled ####
# I'll be looking within the df.sim.100 at pho = 1.9, 0.9565 and 0.
sub.100.1.9 <- subset(df.sim.100, df.sim.100$Pho.Value==1.8217)
sub.100.0.96 <- subset(df.sim.100, df.sim.100$Pho.Value==0.8652)
sub.100.0 <- subset(df.sim.100, df.sim.100$Pho.Value==0.1000)

# Building the histogram 
hist.100.1.9 <- hist(sub.100.1.9$Salinity, sub.edge) 
hist.100.0.96 <- hist(sub.100.0.96$Salinity, sub.edge)
hist.100.0 <- hist(sub.100.0$Salinity, sub.edge)
#They look like what is expected

#### FLoquet Exponent Box Plots ####
simulation.data.frame <- df.sim.100 # What data set you wanna look at 
N = 100 # Number of simulatins

nauplii.inter <- simulation.data.frame %>%
  split(list(.$Pho.Value, .$Simulation.Number), lex.order=F) %>%
  map(~lm(.$Nauplii ~ .$Time, na.action=na.exclude)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1) # Calculates the intercept for each pho values for each simulation. 
# Map dbl tells R what value it wants. In this case it holds the value for the intercept of the regression 
# for life stages ~ time

nauplii.slope <- simulation.data.frame %>%
  split(list(.$Pho.Value, .$Simulation.Number), lex.order = F) %>%
  map(~lm(.$Nauplii ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2) # Calculates the slope/floquet exponent for each pho values for each simulation. 
# Map dbl tells R what value it wants. In this case it holds the value for the intercept of the regression 
# for life stages ~ time

adult.inter <- simulation.data.frame %>%
  split(list(.$Pho.Value, .$Simulation.Number), lex.order = F) %>%
  map(~lm(.$Adult.Females ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1)

adult.slope <- simulation.data.frame %>%
  split(list(.$Pho.Value, .$Simulation.Number), lex.order = F) %>%
  map(~lm(.$Adult.Females ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2)

# Combining all the datasets together to ger the floquet exponets
nauplii.floquet <- as.data.frame(cbind(nauplii.inter, nauplii.slope))
nauplii.names <- data.frame(row.names(nauplii.floquet))
nauplii.names  <- nauplii.names %>%
  separate(row.names.nauplii.floquet., c("pho_ones", "tenths", "simulation_number")
  ) 
nauplii.pho <- paste(nauplii.names$pho_ones, nauplii.names$tenths, sep=".")
nauplii.floquet <- cbind(nauplii.pho, nauplii.names$simulation_number, nauplii.floquet)
names(nauplii.floquet) <- c("Pho_value", "Sim. #", "Intercept", "Floquet")


# # Combining copepodids abundance
# copepodids.floquet <- as.data.frame(cbind(copepodids.inter, copepodids.slope))
# copepodids.names <- data.frame(row.names(copepodids.floquet))
# copepodids.names <- copepodids.names %>%
#   separate(row.names.copepodids.floquet., c("pho_ones", "tenths", "simulation_number")
#   )
# copepodids.pho <- paste(copepodids.names$pho_ones, copepodids.names$tenths, sep=".")
# copepodids.floquet <- cbind(copepodids.pho, copepodids.names$simulation_number, copepodids.floquet)
# names(copepodids.floquet) <- c("Pho_value", "Sim. #", "Intercept", "Floquet")
# 
# # Combing chalimi abundnace
# chalimi.floquet <- as.data.frame(cbind(chalimi.inter, chalimi.slope))
# chalimi.names <- data.frame(row.names(chalimi.floquet))
# chalimi.names <- chalimi.names %>% 
#   separate(row.names.chalimi.floquet., c("pho_ones", "tenths", "simulation_number")
#   )
# chalimi.pho <- paste(chalimi.names$pho_ones, chalimi.names$tenths, sep=".")
# chalimi.floquet <- cbind(chalimi.pho, chalimi.names$simulation_number, chalimi.floquet)
# names(chalimi.floquet) <- c("Pho_value", "Sim. #", "Intercept", "Floquet")

# Combing adult abundance
adult.floquet <- as.data.frame(cbind(adult.inter, adult.slope))
adult.names <- data.frame(row.names(adult.floquet))
adult.names <- adult.names %>%
  separate(row.names.adult.floquet., c("pho_ones", "tenths", "simulation_number")
  )
adult.pho <- paste(adult.names$pho_ones, adult.names$tenths, sep=".")
adult.floquet <- cbind(adult.pho, adult.names$simulation_number, adult.floquet)
names(adult.floquet) <- c("Pho_value", "Sim. #", "Intercept", "Floquet")

# Looking at the head of each data frame
head(nauplii.floquet)
head(copepodids.floquet)
head(chalimi.floquet)
head(adult.floquet)

# Removing the last two pho values from the regression analysis
nauplii.reg <- filter(nauplii.100, nauplii.100$Pho_value > 0.4)
cop.reg <- filter(cop.100, cop.100$Pho_value > 0.4)
chali.reg <- filter(chalimi.100, chalimi.100$Pho_value > 0.4)
adult.reg <- filter(adult.100, adult.100$Pho_value > 0.4)

# Getting slope for the nauplii graph 
coef.naup <- coef(lm(nauplii.reg$Floquet ~ rev(factor(nauplii.reg$Pho_value)), data = nauplii.reg))
ggplot(nauplii.floquet, aes( x= factor(nauplii.floquet$Pho_value), 
                                 y = nauplii.floquet$Floquet)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Salinity Variation", y="Floquet Exponent", 
       title = paste("Nauplii Floquet:",
                     "Intercept =",signif(coef.naup[[1]],2 ),
                     ", Slope =", signif(coef.naup[[2]], 2))) +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.2))  + 
  geom_abline(intercept = coef.naup[1], slope = coef.naup[2], color = "red") + 
  scale_x_discrete(limits = rev(levels(nauplii.floquet$Pho_value)))


coef.cop <- coef(lm(Floquet ~ rev(factor(Pho_value)), data=cop.reg))
ggplot(copepodids.floquet, aes(x=factor(copepodids.floquet$Pho_value), 
                                    y=copepodids.floquet$Floquet)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Salinity Variation", y="Floquet Exponent", 
       title=paste ("Copepodids Floquet:", 
                    "Intercept =", signif(coef.cop[[1]],2 ),
                    ", Slope =", signif(coef.cop[[2]], 2))) +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_abline(intercept = coef.cop[1], slope = coef.cop[2], color="red") +
  scale_x_discrete(limits = rev(levels(copepodids.floquet$Pho_value)))


coef.cham <- coef(lm(Floquet ~ rev(factor(Pho_value)), data=chali.reg))
ggplot(chalimi.floquet, aes(x=factor(Pho_value), y=Floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Salinity Variation", y="Floquet Exponent", 
       title= paste("Chalimis Floquet:", 
                    "Intercept =", signif(coef.cham[[1]], 2), 
                    ", Slope =", signif(coef.cham[[2]], 2))) +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  geom_abline(intercept = coef.cham[1], slope = coef.cham[2], color = "red") +
  scale_x_discrete(limits = rev(levels(chalimi.floquet$Pho_value)))


coef.adult <- coef(lm(adult.reg$Floquet ~ rev(factor(adult.reg$Pho_value)), data=adult.reg))
ggplot(adult.floquet, aes(x=factor(Pho_value), y=Floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Scale Parameter Value", y="Floquet Exponent") +
  ylim(0, 5.4e+05) +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  geom_abline(intercept = coef.adult[1], slope = coef.adult[2], color = "red") +
  scale_x_discrete(limits = rev(levels(adult.floquet$Pho_value))) + 
  annotate("text", x = c(mean.floquet.per.pho$Group.1 ), y = 5.1e+05, 
           label = c(rev(mean.floquet.per.pho$x)), angle=45)



mean.floquet.per.pho <- aggregate(adult.floquet$Floquet,
                                  list(adult.floquet$Pho_value), mean)
mean.floquet.per.pho$x <- format(round(mean.floquet.per.pho$x, 2),nsmall=2)

# Looking into the plots 
plot(adult.floquet$Pho_value, adult.floquet$Floquet) 
plot(factor(adult.reg$Pho_value), adult.reg$Floquet)

# Looking at the abline that  comes from the bacis lm 
reg1 <- lm(adult.floquet$Floquet ~ as.numeric(adult.floquet$Pho_value))
summary(reg1)
plot(adult.floquet$Pho_value, adult.floquet$Floquet) 
abline(reg1, col="red") # Looks fine, time to rest with reg

# Plotting the subset dataset 
reg2 <- lm(adult.reg$Floquet ~ (as.numeric(adult.reg$Pho_value)))

plot(x=(adult.floquet$Pho_value), y=adult.floquet$Floquet)
abline(reg2, col="red")

plot(x=as.numeric(adult.floquet$Pho_value), y=adult.floquet$Floquet)
abline(reg2, col="red")

# Comparing summary
summary(reg1)
summary(reg2)

# Full plot using smaller pho values reg
plot(x=adult.floquet$Pho_value, y=adult.floquet$Floquet)
abline(reg1, col="red")

plot(x=as.numeric(adult.floquet$Pho_value), y=adult.floquet$Floquet)
abline(reg1, col="red")


#### ANOVAs for the life-stages ####
nauplii.anova <- aov(nauplii.100$Floquet ~ as.factor(nauplii.100$Pho_value), data=nauplii.100)
cop.anova <- aov(cop.100$Floquet ~ as.factor(cop.100$Pho_value), data=cop.100)
chali.anova <- aov(chalimi.100$Floquet ~ as.factor(chalimi.100$Pho_value), data=chalimi.100)
adult.anova <- aov(adult.100$Floquet ~ as.factor(adult.100$Pho_value), data=adult.100)

# Summary output 
summary(nauplii.anova)
summary(cop.anova)
summary(chali.anova)
summary(adult.anova)

# Can the different in pho values also be used to decribe salinity and abundance?
sal.anova <- aov(df.sim.100$Salinity ~ as.factor(df.sim.100$Pho.Value), data=df.sim.100)
summary(sal.anova)

# Abundnace 
naup.abund.anova <- aov(df.sim.100$Nauplii ~ as.factor(df.sim.100$Pho.Value), data=df.sim.100)
cop.abund.anova <- aov(df.sim.100$Copepodids ~ as.factor(df.sim.100$Pho.Value), data=df.sim.100)
chali.abund.anova <- aov(df.sim.100$Chalimi.and.pre.adults ~ as.factor(df.sim.100$Pho.Value), data=df.sim.100)
adult.abund.anova <- aov(df.sim.100$Adult.Females ~ as.factor(df.sim.100$Pho.Value), data=df.sim.100)

# Summary output 
summary(naup.abund.anova)
summary(cop.abund.anova)
summary(chali.abund.anova)
summary(adult.abund.anova)

# All the summary output indicate that there is a strong and statistical significant different between the floquet exponets as the pho values change

#### Things to look into ####
# - All nauplii see to have an opposite trend/unexpected trend compared to the rest of the adult stage 
# - Need to update the word document 


# Producing ggplot that shows both the salinity values and the log abundance. 
# Want pho values of 0.6739 and 1.8217. In the future it is important to look at pho values greater then baseline 


# log abundannce plot

# Need to seperate the dataset into two minor ones with pho values 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colourblind palette

q <-1
sample(seq(1,100), 1)

sim.10y <- read.csv("../../../../../../../../Desktop/sim.10y 2.txt")

simulation.data.frame <- sim.10y
simulation.data.frame <- as.data.frame(simulation.data.frame)

adult.pho.18 <- filter(simulation.data.frame, 
                       simulation.data.frame$Pho.Value==1.8217 & 
                         simulation.data.frame$Simulation.Number== 1)
adult.pho.0.6739 <- filter(simulation.data.frame, 
                                  simulation.data.frame$Pho.Value ==0.6739 & 
                                  simulation.data.frame$Simulation.Number== 1)

adult.pho.sub <- rbind(adult.pho.18, adult.pho.0.6739)
adult.pho.sub$Pho.Value <- factor(adult.pho.sub$Pho.Value)
adult.pho.18$Pho.Value <- factor(adult.pho.18$Pho.Value)
adult.pho.0.6739$Pho.Value <- factor(adult.pho.0.6739$Pho.Value)
# Generating the `` 
attach(adult.pho.sub)



adult.abundance <- ggplot(adult.pho.sub, aes(x= Time,  y=log(Adult.Females))) +
  geom_line(aes(colour=factor(Pho.Value)), size=1.5) +
  labs(x="Time (days)", y="Logarithmic Abundance of Adult Females") +
  theme_classic() +  
  scale_color_manual(name = "Scale \nParameter", 
                     breaks = c("1.8217", "0.6739"),
                     values=c("#FF0000", "#00FFFF")) + 
  theme(legend.position="NA") 

adult.salinity <- ggplot(adult.pho.sub, aes(x= Time,  y= Salinity, 
                                            colour=factor(Pho.Value))) +
  geom_line(aes(), size=1.5) +
  labs(x="Time (days)", y="Salinity (psu)") +
  theme_classic() +  
  scale_color_manual(breaks = c("1.8217", "0.6739"),
                     values=c("#FF0000", "#00FFFF")) + 
  theme(legend.position="NA") 


ratio.plot <- ggplot(adult.pho.sub, aes(x= Time,  y=(Nauplii/Copepodids))) +
  geom_line(aes(colour=Pho.Value), size=1) +
  xlim(7*365, 3650) +
  labs(x="Times (days)", y="Ratio of Nauplii to Adult Females") +
  theme_classic() +  
  scale_color_manual(name = bquote("Scale in daily\nsalinity,"~lambda), breaks = c("1.8217", "0.6739"),
                                                  values=c("#FF0000", "#00FFFF"))

plot_grid(adult.salinity, adult.abundance, ratio.plot, labels = c("(a)", "(b)", "(c)"), 
          nrow = 1)
# Now to put them together using the patchwork package

max(sub.100.0$Time)


# Want to calculat the mean floquet exponents for all simulations 

aggregate(adult.floquet[, 3:4], list(adult.floquet$Pho_value), mean)


# Correlation?
cor(adult.100$Pho_value, adult.100$Floquet)
