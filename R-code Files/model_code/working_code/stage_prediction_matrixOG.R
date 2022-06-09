# Sample Code to create and predict data that will be recorded into a dataset. 
rm(list=ls())

# Setting wd and loading additional Rdata files
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/first_chapter/winter_2020/model_code/data_files/amy_regional/")
load("ParamsFuncs.RData")
load("farmNL.RData")

# Needed packages ---------------------------------------------------------------
# install.packages(c("tidyverse", "stringr", "dplyr")) # If packages need to be installed
x <-c("tidyverse", "stringr", "dplyr", "base", 
        "ggplot2", "patchwork", "bbmle")
lapply(x, require, character.only = TRUE)


# Loading Bay D'Esspoir dataset
baydespoir_prof <- read.csv("../baydespoir.nosync/baydespoir_prof.csv") 

# Setting data frame that will be used in the loop. 
simulation.data.frame=NULL

# Formatting it into joulain time ----
# The date format that is needed for the chron data 
monthly.dates <- dates(paste(baydespoir_prof$OBS_MONTH, baydespoir_prof$OBS_DAY, 
                             baydespoir_prof$OBS_YEAR,
                             sep="/")) 

# The date format that ggplot used 
baydespoir_prof$ggplot.data <- as.Date(paste(baydespoir_prof$OBS_YEAR, baydespoir_prof$OBS_MONTH, 
                                        baydespoir_prof$OBS_DAY,
                                        sep="-")) 

baydespoir_prof$chron_date <-chron(dates=monthly.dates, 
                              origin. = c(month = 01,day = 01,
                                          year = 1956)) # setting into chron
baydespoir_prof <- arrange(baydespoir_prof, chron_date) # ordering dates 

o <- as.Date("1956-01-01") # The origin or first date
baydespoir_prof$days_since_origin <- (as.numeric(as.Date(baydespoir_prof$chron_date) - o)) # Producing a seqfrom smallest to largest

# Formating the data (only want data that has recorded temperature)
profile.df <- filter(profile.df, !is.na(TEMP)) 
profile.df.sub <- filter(profile.df, LATITUDE...N. > 46.1 & LATITUDE...N. < 48.04)
profile.df.sub <- filter(profile.df.sub, LONGITUDE...E. < -54 & LONGITUDE...E. > -59.5)


# Plot of Bay d'Espoir ---------
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  geom_point(data= baydespoir_prof, aes(x=LONGITUDE...E., y=LATITUDE...N.)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(4, "in"), pad_y = unit(4, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -56.07, y = 47.40, label = "Bay d'Espoir",
           fontface = "bold", color = "black", size = 4) +
  coord_sf(xlim = c(-58, -54), ylim = c(46, 49), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
  ) +
  labs(x="", y="")


# Plotting the empirical data ---------
# Plotting the temperature data 
baydespoir_prof.sub <- subset(baydespoir_prof, baydespoir_prof$OBS_YEAR > 1999)

ggplot(baydespoir_prof.sub, aes(x=ggplot.data, y=TEMP)) +
  geom_line(aes(colour=factor(OBS_YEAR)))  + 
  theme_classic() +
  scale_x_date(date_labels = "%Y %b %d", date_breaks = "2 year") +
  xlab("") +
  ylab("Temperature (celsius)") + 
  labs(colour = "Observation Year") +
  theme(legend.position="top", 
        axis.text.x=element_text(size= 10), 
        axis.text.y=element_text(size= 10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size = 10), 
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) 

# Plotting the salinity data 
  ggplot(baydespoir_prof, aes(x=OBS_MONTH, y=PSAL)) +
  geom_jitter(aes(colour=factor(OBS_YEAR))) +
  theme_classic() +
  ylim(21, 33) +
  xlab("") +
  ylab("Salinity (psu)") + 
  labs(colour = "Observation Year") +
  theme(legend.position="right",
        axis.text.x=element_text(size= 10), 
        axis.text.y=element_text(size= 10), 
        axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size = 10), 
        legend.text = element_text(size=10), 
        legend.title = element_text(size = 10)) 

vertical.plot <- plot_grid(record.dfo.temp, record.dfo.sal, 
                           labels = c("(b)", ("(c)")), ncol=1)

plot_grid(bay.observation, vertical.plot, labels= c("(a)", ""))

# Loop parameters ---------------------------------------------------------------
# k=10 # Salinity variance for j loop
N=10 # Number of simulation that the loop will rerun for each variacne for v loop
TT= 10*365 # Number of days the model will run

# Time 
dtime = 0.5
time <-  seq(0, TT, dtime)
time0 <-time[1]

# Temperature Functions 
Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
  return(Temp)
}

Temp1 = function(t){
  Temp = 6.4294 + 3.32*sin(2 * pi * t/365) + 4.61*cos(2 * pi * t/365)
  return(Temp)
}

df <- as.data.frame(cbind(seq(0, 365*2), Temp(seq(0, 365*2))))
ggplot(df, aes(x=V1, y=V2)) +
  geom_line() + 
  theme(axis.text.x=element_text(size= 14), 
        axis.text.y=element_text(size= 14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          labs(x = "Time (days)",
               y = "Temperature (celsius)")


Temp1(time.seq)
time.seq <- seq(0, 365, 0.5)
Temp(time)

# time.temp <- as.data.frame(cbind(t, Temp(t)))
# ggplot(time.temp, aes(x=t, y=Temp(t))) +
#   geom_line()+
#   geom_point() +
#   theme_classic() +
#   labs(x="Time (days)", y="Ocean Temperature (Celsius)")
# 
# plot(t, Temp(t), type="l", xlab="Time (days)",
#      ylab="Ocean Temperature (Celsius)", bty="L", lwd = 5,
#      cex.axis=1.25, cex.lab=1.25)
# axis(1:2)

# Life stage initial conditions
P0 <- 1
C0 <- 1
A0 <- 1 
I0 <- 1
P = P0
C = C0
A = A0
I = I0

# Paramater conditions 
iota <- iota1 # Attachment rate of lice
f <- 500000 # Number of fish within each pen

#  Asymetric Laplace distrubution -----------------------------------------------
# Using the Asymmetric Laplace distrubution to produce a cummulative distrubution 
# Changing the bin and x-vectors
dx <- 0.5
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Subsetting to only take values that are greater then 15
sub.bay <- subset(baydespoir_prof, PSAL >= 15)
sub.bay <- na.omit(sub.bay$PSAL)

# Plotting the subset distrubution 
sub.hist <- hist(sub.bay, sub.edge) # Plotting the histogram when subsetted from 15-40
sub.freq <- sub.hist$counts
sub.freq <- sub.freq/sum(sub.freq) # Converting the count data into frequency
sum(sub.freq) # Should sum to 1

# Checking to make sure that they're the same length 
length(sub.freq)
length(sub.x)

# Find the m value 
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2 
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1

# Need to update the asymmetric laplace distrubution so that it pulls the right 
# x value

sub.asy.laplace <- function(sub.m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
} # Same distrubution but x was changed to sub.x. Sub.x is the x-vector for the subset data

sub.negLL = function(pho, k){
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
} # Only want to fit to pho and k so only include them within the code

sub.mle.asy.laplace.Ex <- mle2(sub.negLL, 
                               start=list(pho=1, k=3), ) # Similary, only want to fit to pho and k

a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]

# fitted asymmetric distrubution
print(a.est)

# Initial fitted values 
pho0 <- a.est[1]
k0 <- a.est[2]

# # The probability distrubution is described by
plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", ylim=c(0,1.1)) # Probability 
#lines(sub.x, cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2])), col="red") # Cumlative
points(sub.x, sub.freq, col="blue", pch=20)

# Combing the distrubutions into one distrubution 
sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]),
                                  cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2]))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
head(sub.sal.df)

# Pulling a probability
simulation.data.frame<- as.data.frame(matrix(NA, nrow=0, ncol=9))
# Function is based off Hurford_etal_no_delay function --------------------------
for(pho in c(0.4826, 0.6739,
             1.8217, 4)){  # Salinity values (rounded up)
  
  # Initial value
  sub.m <- 32.25 # Should stay 32.25
  pho <- pho # Should follow the vector above
  k <- 1.122308 # Should stay 1.122308

  sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, pho, k),
                                    cumsum(sub.asy.laplace(sub.m, pho, k))))
  names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
  
  Temp = function(t){
    Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
    return(Temp)
  }
  
  Temp0 <- Temp(1)

  # Salinity Value
  salinity0 = 32.25 # Need to define salinity0 within the loop
  
  for(v in seq(1, N)){ # The number of simulations that will run 
    
    # Setting initial conditions again
    S <- sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF)
    P = P0
    C = C0
    A = A0
    I = I0
    
    result <- as.data.frame(matrix(c(pho, v, time0, Temp0, S, P, C, A, I),ncol=9, nrow=1))
    simulation.data.frame <- bind_rows(simulation.data.frame, result)
    
    for(t in time[2:length(time)]){ # Populates the matrix for each time step (each day a new estimate)
      
      # Random variable where the loop reruns v 10 times 
      S = sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF)
      
      # Stage
      P = max((P + (eta(Temp(t))*epsilon(Temp(t))*A*V(Temp(t), S) - mup(S)*P - gammaP(Temp(t))*P)*dtime), 0)
      I = max((I + (gammaP(Temp(t))*P -iota*f*I - mui(S)*I)*dtime), 0)
      C = max((C + (iota*f*I - gammaC(Temp(t))*C - muc(S)*C)*dtime), 0)
      A = max((A + (gammaC(Temp(t))*C - mua(S)*A)*dtime), 0)
      
      # Row binding everything into a single data frame
      result2 <- as.data.frame(matrix(c(pho, v, t, Temp(t), S, P, C, A, I), ncol=9, nrow=1))
      simulation.data.frame <- bind_rows(simulation.data.frame, result2)
      
    } # End of the time loop
  }  # End of the v loop
} # End of the j salinity loop

simulation.data.frame <- as.data.frame(simulation.data.frame)

names(simulation.data.frame) <- c( "Pho Value", "Simulation Number", "Time", "Salinity", 
                                   "Nauplii", "Chalimi and pre-adults", "Adult Females", 
                                   "Copepodids")

simulation.data.frame$`Pho Value` <- factor(simulation.data.frame$`Pho Value`)
str(simulation.data.frame)

sim.10y <- simulation.data.frame 

# The code to generate the Floquet exponents is actually in another code. 
# It's within a R code called auto_floq_code
sim.10y <- simulation.data.frame
sim.10y$
write.csv(sim.10y, "../../model_figures/model_figure_data/sim.10y.txt")

sim.data <- read.csv("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/first_chapter/winter_2020/model_code/model_figures/model_figure_data/sim_data_100")

# Code to calculate the Floquet Exponent
# Calculating the Floquet Exponent when the scale parameter varies 
sim.data$Pho.Value <- as.factor(sim.data$Pho.Value)
#simulation.data.frame$Pho.Value <- as.factor(simulation.data.frame$`Pho Value`)

comb.df <- rbind(sim.data, simulation.data.frame)

#sim.data # The simulated predicted abundance of salmon lice made by the model
A0 <- 1 # The starting abundance of adult females within the system at t=0
C0 <- 1
I0 <- 1
P0 <- 1

Floq.data <- NULL # The dataframe that will hold the calculated Floquet exponents
names(sim.data) <- c( "Pho Value", "Simulation Number", "Time", "Salinity", 
                      "Nauplii", "Chalimi and pre-adults", "Adult Females", 
                      "Copepodids")
simulation.data.frame <- sim.data

for(p in c(0.1, 0.2913, 0.4826, 0.6739,
           0.8652, 1.0565, 1.2478, 1.4391,
           1.6304, 1.8217)){
  test.pho <- subset(simulation.data.frame, simulation.data.frame$`Pho Value`==p) # separates sim.data matrix into smaller systems
  
  for (z in seq(1,100)){
    test.sub <- subset(test.pho, test.pho$`Simulation Number`==z) # seperates and calculates the
    # Floquet theory for each of the n=100 simulations 
    
    SpectralBond = function(lambda, g){
      index <- NULL
      pastvals <- matrix(c(rep(A0, TT+1)),
                         ncol=1, nrow=(TT+1)) # Produces a vector of A0 of length = TotalTime (TT)
      
      for(i in seq(1,100)){
        # Calculates 
        index[i] = max(c(max(test.sub$`Adult Females`)))# Finds the max abundance of adult females
        
        A0 <- tail(test.sub$`Adult Females`, n=1)/index[i]
        
        pastvals=matrix(c(test.sub$`Adult Females`/index[i],
                          ncol = 1, nrow=(TT+1)
        ))
        
        # if(i>1 && (abs(index[i]-index[i-1])<1e-5)){
        #   break
        # }
        
        
        r=index[i]-1
        return(r)
      }
    }
    # Now calculating the Floquet
    FEst = log(SpectralBond(1,1)+1)/365
    Floq.data <-rbind(Floq.data, c(FEst, z, p)) # combines the calculated Floquet exponents, 
    # simulation number, and the scale parameter
  }
}
Floq.data <- as.data.frame(Floq.data)
names(Floq.data) <- c("floquet", "sim", "pho")
head(Floq.data)

# Building a regression group of data. Removes the last two groups
Floq.data$pho <- factor(Floq.data$pho)
adult.flo.reg <- filter(Floq.data, Floq.data$pho != c("0.1000", "0.2913"))

floq.mean <- aggregate(Floq.data$floquet, list(Floq.data$pho), mean)[,2]
floq.mean <- format(round(floq.mean, 5), nsmall = 2)
floq.mean

# Now looking into the Floquet for the skewness 
sim.data <- read.csv("../../model_figures/model_figure_data/skewness_100")

A0 <- 1 # The starting abundance of adult females within the system at t=0
C0 <- 1
I0 <- 1
P0 <- 1
simulation.data.frame <- sim.data
TT.time <- 365 + 1
Floq.data.skew <- NULL
for(p in c(0.1, 1.2, 2.3, 3.4, 4.5, 5.6,
           6.7, 7.8, 8.9, 10)){
  test.pho <- subset(simulation.data.frame, simulation.data.frame$Skewness ==p)
  
  for (z in seq(1, 100)){ 
    # z loop controls the simulation number
    test.sub <- subset(test.pho, test.pho$Simulation.Number ==z) # Pulling a subset data frame that resets each
    # interaction
    
    SpectralBond = function(lambda, g){
      index <- NULL
      pastvals <- matrix(rep(A0, TT.time), ncol=1, nrow=(TT+1))
      
      for(i in seq(1,10)){
        index[i] = (max(test.sub$Adult.Females))
        
        A0 <- tail(test.sub$Adult.Females, n=1)/index[i]
        
        pastvals=matrix(c(test.sub$Adult.Females/index[i]),
                        ncol = 1, nrow=(TT+1))
        r=index[i]-1
        return(r)
      }
    }
    # Now calculating the Floquet 
    FEst = log(SpectralBond(1,1)+1)/TT.time
    Floq.data.skew <-rbind(Floq.data.skew, c(FEst, z, p))
  }
}
Floq.data.skew <- as.data.frame(Floq.data.skew)
names(Floq.data.skew) <- c("floquet", "sim", "skewness")
print(Floq.data.skew)

floq.skew <- aggregate(Floq.data.skew$floquet, list(Floq.data.skew$skewness), mean)[,2]
floq.skew <- format(round(floq.skew, 5), nsmall = 2)

ggplot(Floq.data.skew, aes(x=factor(skewness), y=floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Groups", y="Response Value") +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_x_discrete(limits = rev(levels(Floq.data.skew$skewness))) + 
  annotate("text", x=1:10, y = 0.07, label = floq.skew, angle = 45, size=4) +
  geom_hline(yintercept = 0.05681907, color = "red") + 
  labs(x= expression("Skewness in daily salinity," ~italic(k)),
       y= expression("Floquet exponent"~ (`phi`)))



# Building a regression group of data. Removes the last two groups
Floq.data$pho <- factor(Floq.data$pho)
adult.flo.reg <- filter(Floq.data, Floq.data$pho != c("0.1000", "0.2913"))

floq.mean <- aggregate(Floq.data$floquet, list(Floq.data$pho), mean)[,2]
floq.mean <- format(round(floq.mean, 5), nsmall = 2)
floq.mean



# Now producing the Floquet graph
coef.adult <- coef(lm(floquet ~ rev(factor(pho)), data=adult.flo.reg))
coef.full.adult <- coef(lm(floquet ~ rev(factor(pho)), data=Floq.data))

# Producing figure 5 
# Some of the code to produce the effects of the various scale parameters 

# Find the m value
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1
pho0 <- 1.9137
k0 <- 1.223080

sub.asy.laplace <- function(sub.m, pho, k){
  d <- 0.45*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
} # Same distrubution but x was changed to sub.x. Sub.x is the x-vector for the subset data

sub.dx <-
  # Produce the asymmetric laplce with the dx
  sub.asy.laplace.dx <- function(sub.m, pho, k){
    d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
  }

sub.negLL = function(pho, k){
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
} # Only want to fit to pho and k so only include them within the code

sub.mle.asy.laplace.Ex <- mle2(sub.negLL,
                               start=list(pho=1, k=3), ) # Similary, only want to fit to pho and k

a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]
print(a.est)

# The probability distrubution is described by
# Combing the distrubutions into one distrubution
sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]),
                                  cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2]))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
head(sub.sal.df)


hist.sal.plot <- as.data.frame(cbind(sub.x, sub.bay/sum(sub.bay), sub.bay, 
                                     sub.edge,sub.asy.laplace(sub.m, 1.91, 1.22), rep(1.8217, 50)))
names(hist.sal.plot) <- c("Salinity", "Frequency", "Recorded Sal", "Bin Edge", "Probability", "Scale Parameter")

# Building a plot that highlights the difference between the different pho values and salintiy
# Creating a new data frame for the ggplot
var.salinity.dist.df <- as.data.frame(
  cbind(rep(sub.x, 11),
        
        c(sub.asy.laplace(sub.m, 4, k0), sub.asy.laplace(sub.m, 3.5, k0), sub.asy.laplace(sub.m, 3, k0),
          sub.asy.laplace(sub.m, 2.5, k0), sub.asy.laplace(sub.m, 2, k0), sub.asy.laplace(sub.m, 1.8217, k0),
          sub.asy.laplace(sub.m, 1.4391, k0), sub.asy.laplace(sub.m, 1.0565, k0), sub.asy.laplace(sub.m, 0.6739, k0),
          sub.asy.laplace(sub.m, 0.4826, k0), sub.asy.laplace(sub.m, 0.2913, k0)),
        
        c(rep(4, 50), rep(3.5, 50), rep(3, 50), rep(2.5, 50), rep(2, 50), rep(1.8217, 50), rep(1.4391, 50),
          rep(1.0565, 50), rep(0.6739, 50), rep(0.4826, 50), rep(0.2913, 50))
  ))

names(var.salinity.dist.df) <- c("Salinity", "Probability", "Scale Parameters")
var.salinity.dist.df$`Scale Parameters` <- as.factor(var.salinity.dist.df$`Scale Parameters`)

revalue(var.salinity.dist.df$`Scale Parameters`, c("0.2913" = "0.29", 
                                             "0.4826" = "0.48", 
                                             "0.6739" = "0.67",
                                             "1.0565" = "1.06", 
                                             "1.4391" = "1.44",
                                             "1.8217" = "1.82", 
                                             "2" = "2.00", 
                                             "2.5" = "2.50", 
                                             "3" = "3.00", 
                                             "3.5"= "3.50", 
                                             "4"= "4.00")) -> var.salinity.dist.df$`Scale Parameters`

var.sal.plot<- ggplot(var.salinity.dist.df, aes(x=Salinity, y=Probability, 
                                   colour=factor(`Scale Parameters`))) +
  geom_line() + 
  theme_classic() +
  theme(axis.text.x=element_text(size= 14), 
        axis.text.y=element_text(size= 14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size = 14))+
  scale_colour_manual(values = rev(c("black", "brown", "blue", "purple", "pink",
                                     "red", "orange", "yellow", "cyan", "green", 
                                     "darkgreen")), 
                      name = bquote('Scale in\n daily salinity '~lambda))  +
  xlim(26, 35) +
  ylim(0, 1.1) +
  labs(x = "Salinity (psu)",
       y = "Probability") 

# Code for the random salinity 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colourblind palette

q <-1
sample(seq(1,100), 1)

sim.data <- read.csv("../../../../../../../sim.10y 2.txt")
simulation.data.frame <- sim.data
simulation.data.frame <- as.data.frame(simulation.data.frame)
# names(simulation.data.frame) <- c( "Pho Value", "Simulation Number", "Time", "Temp", "Salinity", 
#                                    "Nauplii", "Chalimi and pre-adults", "Adult Females", 
#                                    "Copepodids")
simulation.data.frame$Pho.Value <- as.factor(simulation.data.frame$Pho.Value)

revalue(simulation.data.frame$Pho.Value, c("0.1"="0.10", 
                                             "0.2913" = "0.29", 
                                             "0.4826" = "0.48", 
                                             "0.6739" = "0.67", 
                                             "0.8652" = "0.87", 
                                             "1.0565" = "1.06", 
                                             "1.2478" = "1.25", 
                                             "1.4391" = "1.44", 
                                             "1.6304" = "1.63", 
                                             "1.8217" = "1.82")) -> simulation.data.frame$Pho.Value
revalue(Floq.data$pho, c("0.1"="0.10", 
                                             "0.2913" = "0.29", 
                                             "0.4826" = "0.48", 
                                             "0.6739" = "0.67", 
                                             "0.8652" = "0.87", 
                                             "1.0565" = "1.06", 
                                             "1.2478" = "1.25", 
                                             "1.4391" = "1.44", 
                                             "1.6304" = "1.63", 
                                             "1.8217" = "1.82")) -> Floq.data$pho


adult.pho.18 <- filter(simulation.data.frame, simulation.data.frame$Pho.Value == 1.82 &
  simulation.data.frame$Simulation.Number== 1)
adult.pho.0.6739 <- filter(simulation.data.frame, simulation.data.frame$Pho.Value == 0.67 &
                             simulation.data.frame$Simulation.Number == 1)

adult.pho.sub <- rbind(adult.pho.18, adult.pho.0.6739)
adult.pho.sub$Pho.Value <- factor(adult.pho.sub$Pho.Value)
adult.pho.18$Pho.Value <- factor(adult.pho.18$Pho.Value)
adult.pho.0.6739$Pho.Value <- factor(adult.pho.0.6739$Pho.Value)

adult.salinity <- ggplot(adult.pho.sub, aes(x= Time,  y= Salinity)) +
  geom_line(aes(colour=Pho.Value), size=1.5) +
  labs(x="Time (days)", y="Salinity (psu)") +
  ylim(0, 43) +
  theme_classic() +  
  scale_color_manual(breaks = c("1.82", "0.67"),
                     values=c("#FF0000", "#00FFFF"), 
                     name=bquote('Scale in\n daily salinity '~lambda)) + 
  theme(legend.position="right") 


floquet.scale.parameter <- 
  ggplot(Floq.data, aes(x=factor(pho), y=floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Groups", y="Response Value") +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle=45, hjust=1))+
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_x_discrete(limits = rev(levels(Floq.data$pho))) + 
  annotate("text", x=1:10, y = 0.0650, label = rev(floq.mean), angle = 45, size=4) +
  geom_hline(yintercept = 0.05681907, color = "red") + 
  labs(x= expression("Scale in daily salinity,"~ `lambda`),
       y= expression("Floquet exponent"~ (`phi`)))
  

vertical.plots <- plot_grid(var.sal.plot, adult.salinity, labels = c("(a)", "(c)"),
                            ncol = 1, nrow=2)
plot_grid(vertical.plots, floquet.scale.parameter, labels= c("", "(b)"))





TukeyHSD((aov(floquet ~ pho, data=Floq.data))) # Cleary the phos are very different
anova(L0(floquet ~ pho, data=Floq.data))



# Comparing the difference between a constant salinity and the first pho
# First need to calculate the Floquet exponents
A0 <- 1

const.sal.flo <- NULL
for (z in seq(1, N)){ 
  # z loop controls the simulation number
  test.sub <- subset(simulation.data.frame, 
                     simulation.data.frame$`Simulation Number` ==z) # Pulls each simulation number from the df
  # interaction
  
  SpectralBond = function(lambda, g){
    index <- NULL
    pastvals <- matrix(rep(A0, TT), ncol=1, nrow=TT)
    for(i in seq(1,10)){
      index[i] = (max(test.sub$`Adult Females`))
      A0 <- tail(test.sub$`Adult Females`, n=1)/index[i]
      pastvals=matrix(c(test.sub$`Adult Females`/index[i]),
                      ncol = 1, nrow=TT)
      r=index[i]-1
      return(r)
    }
  }
  # Now calculating the Floquet 
  FEst = log(SpectralBond(1,1)+1)/365
  const.sal.flo <-rbind(const.sal.flo, c(FEst, z, 1))
}

const.sal.flo <- as.data.frame(const.sal.flo)
names(const.sal.flo) <- c("floquet", "sim", "pho")
const.sal.flo$pho <- factor(const.sal.flo$pho)

# Only want the first flow. 

floquet.pho.1.8 <- filter(Floq.data, Floq.data$pho == 1.8217)

const.sal.pho.1.8 <- rbind(const.sal.flo, floquet.pho.1.8)

summary(aov(floquet ~ pho, const.sal.pho.1.8))

const.sum <- sum(const.sal.flo$floquet)
pho.sum <- sum(floquet.pho.1.8$floquet)
pho.sum/const.sum
const.sum/pho.sum  
# Finding the mean of floquet based on mean 
floq.mean <- aggregate(Floq.data$floquet, list(Floq.data$pho), mean)[,2]
floq.mean <- format(round(floq.mean, 5), nsmall = 2)
floq.mean
attach(Floq.data)



