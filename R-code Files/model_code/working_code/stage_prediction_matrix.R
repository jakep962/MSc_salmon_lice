# Sample Code to create and predict data that will be recorded into a dataset. 
rm(list=ls())

# Setting wd and loading additional Rdata files
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/data_files/amy_regional/")

# Needed packages ---------------------------------------------------------------
# install.packages(c("tidyverse", "stringr", "dplyr")) # If packages need to be installed
x <-c("tidyverse", "stringr", "dplyr", "base", 
        "ggplot2", "patchwork", "bbmle")
lapply(x, require, character.only = TRUE)


# Loading Bay D'Esspoir and data
baydespoir_prof <- read.csv("../baydespoir/baydespoir_prof.csv") 
load("ParamsFuncs.RData")
load("farmNL.RData")

# Setting data frame that will be used in the loop. 
simulation.data.frame=NULL

# Loop parameters 
N=100 # Number of simulation that the loop will rerun for each variacne for v loop
TT= 2*365 # Number of days the model will run

# Time 
dtime = 1
time <-  seq(0, TT, dtime)
time0 <-time[1]

# Temperature Functions 
Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
  return(Temp)
}

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

#  Asymetric Laplace distrubution 
# Produce a cummulative distrubution 
# Changing the bin and x-vectors
dx <- 0.5 # how large the timesteps/bins are
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Subsetting to only take values that are greater then 15
sub.bay <- subset(baydespoir_prof, PSAL >= 15, na.omit=T)
sub.bay <- na.omit(sub.bay$PSAL)

# Plotting the subset distrubution 
sub.hist <- hist(sub.bay$PSAL, sub.edge) # Plotting the histogram when subsetted from 15-40
sub.freq <- sub.hist$counts
sub.freq <- sub.freq/sum(sub.freq) # Converting the count data into frequency
sum(sub.freq) # Should sum to 1

# Find the m value 
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2 
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1 




# Need to update the asymmetric laplace distrubution so that it pulls the right 
# x value

sub.asy.laplace <- function(sub.m, pho, k){
  
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
  
}

sub.negLL = function(pho, k){
  
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
  
} 

sub.mle.asy.laplace.Ex <- mle2(sub.negLL, start=list(pho=1, k=3)) # Similary, only want to fit to pho and k
a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]

# Initial fitted values 
pho0 <- a.est[1]
k0 <- a.est[2]

# Combing the distrubutions into one distrubution 
sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]),
                                  cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2]))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")

# Pulling a probability
sample(sub.sal.df$`Salinity Values`, 1, replace = T, prob=sub.sal.df$PDF)

# Function is based off Hurford_etal_no_delay function --------------------------
for(pho in rev(seq(0.1, 1.913, 1.913/10))){  # Salinity values (rounded up)
  
  # Initial value 
  sub.m <- 32.25 # Should stay 32.25
  pho <- pho # Should follow the vector above
  k <- 1.122308 # Should stay 1.122308
  
  sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, pho, k),
                                    cumsum(sub.asy.laplace(sub.m, pho, k))))
  names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
  
  # Salinity Value
  salinity0 = sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF) # Need to define salinity0 within the loop
  
  for(v in seq(1, N)){ # The number of simulations that will run 
    
    # Setting initial conditions again
    S <- salinity0[1]
    P = P0
    C = C0
    A = A0
    I = I0
    
    simulation.data.frame <- rbind(simulation.data.frame,c(pho, v, time0, S, P, C, A, I))
    
    for(t in time[2:length(time)]){ # Populates the matrix for each time step (each day a new estimate)
      
      # Random variable where the loop reruns v 10 times 
      S = sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF)
      
      # Stage
      P = max((P + (eta(Temp(t))*epsilon(Temp(t))*A*V(Temp(t), S) - mup(S)*P - gammaP(Temp(t))*P)*dtime), 0)
      I = max((I + (gammaP(Temp(t))*P -iota*f*I - mui(S)*I)*dtime), 0)
      C = max((C + (iota*f*I - gammaC(Temp(t))*C - muc(S)*C)*dtime), 0)
      A = max((A + (gammaC(Temp(t))*C - mua(S)*A)*dtime), 0)
      
      # Row binding everything into a single data frame
      simulation.data.frame <- rbind(simulation.data.frame,c(pho, v, t, S, P, C, A, I))
      
    } # End of the time loop
  }  # End of the v loop
} # End of the j salinity loop

simulation.data.frame <- as.data.frame(simulation.data.frame) 
names(simulation.data.frame) <- c( "Pho Value", "Simulation Number", "Time", "Salinity", 
                                   "Nauplii", "Chalimi and pre-adults", "Copepodids", 
                                   "Adult Females")
simulation.data.frame$`Pho Value` <- factor(simulation.data.frame$`Pho Value`)

test <- (aov(`Adult Females` ~ `Simulation Number`, data=simulation.data.frame))
summary(test)

# Comparing now the different between non-variation and future ones 
bind <- filter()



# Model prediction at a constant salinity, will be using sub.m
  for(v in seq(1, 1)){ # The number of simulations that will run 
    constant.sal <- NULL
    # Salinity Value
    salinity0 = sub.m 
    # Setting initial conditions again
    S <- salinity0[1]
    P = P0
    C = C0
    A = A0
    I = I0
    
    constant.sal <- rbind(constant.sal,c(pho, v, time0, S, P, C, A, I))
    
    for(t in time[2:length(time)]){ # Populates the matrix for each time step (each day a new estimate)
      
      # Random variable where the loop reruns v 10 times 
      S = sub.m
      
      # Stage
      P = max((P + (eta(Temp(t))*epsilon(Temp(t))*A*V(Temp(t), S) - mup(S)*P - gammaP(Temp(t))*P)*dtime), 0)
      I = max((I + (gammaP(Temp(t))*P -iota*f*I - mui(S)*I)*dtime), 0)
      C = max((C + (iota*f*I - gammaC(Temp(t))*C - muc(S)*C)*dtime), 0)
      A = max((A + (gammaC(Temp(t))*C - mua(S)*A)*dtime), 0)
      
      # Row binding everything into a single data frame
     constant.sal <- rbind(constant.sal,c(pho, v, t, S, P, C, A, I))
      
    } # End of the time loop
  }  # End of the v loop
constant.sal <- as.data.frame(constant.sal)
names(constant.sal) <- c( "Pho Value", "Simulation Number", "Time", "Salinity", 
                                   "Nauplii", "Chalimi and pre-adults", "Copepodids", 
                                   "Adult Females")
plot(time, log(constant.sal$`Adult Females`))

# Checking to see if the pho values was calculated correctly

box.pho.plot <- ggplot(simulation.data.frame, aes(x=simulation.data.frame$`Pho Value`, 
                                                  y=simulation.data.frame$Salinity, 
                                                  group=simulation.data.frame$`Pho Value`)) +
                geom_boxplot()
box.pho.plot # This is what expected. As the pho value decreases, the spread in the salinity frequencies inscreases

# Subsetting simulation.data.frame by sumulation number -------------------
sep.variance <- split(simulation.data.frame, f=simulation.data.frame$`Pho Value`)

naup.df <- as.data.frame(cbind(sep.variance$`1.8217`$Time, sep.variance$`1`$Nauplii,
                               sep.variance$`1`$`Simulation Number`))


# is this section still needed?

# Plotting the log abundnace of Nauplii life stage against time  ----------
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(sep.variance$`1.8217`$Time, log(sep.variance$`1.8217`$Nauplii), 
     xlab="Time", ylab="Log Abundance")
lines(time, Temp(time), col="red")
axis(side = 4, at = pretty(range(Temp(time)))) 
mtext("Temperature", side=4, line=3, col="red")


# Plotting all the stages on the same graph -------------------------------
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(simulation.data.frame$Time, log(simulation.data.frame$Nauplii),
     type="p", 
     ylab = "log(Abundance per stage)", xlab="Time (day)", main = "Life stages for all var and sim number")
axis(side = 4, at = pretty(range(Temp(time)))) 
mtext("Temperature", side=4, line=3, col="purple")
points(simulation.data.frame$Time, log(simulation.data.frame$Copepodids), 
       type="p", col="blue")
points(simulation.data.frame$Time, log(simulation.data.frame$`Chalimi and pre-adults`), 
       type="p", col="red")
points(simulation.data.frame$Time, log(simulation.data.frame$`Adult Females`), 
       type="p", col="green")
lines(time, Temp(time), col="purple")
legend(1, 20, legend=c("Nauplii", "Copepodids", "Chalimi and pre-adults", "Adult Females", "Temperature"), col=c("black", "blue", "red", "green", "purple"),
       pch=19, cex=0.8)

# Floquet exponent code ---------------------------------------------------
nauplii.inter <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order=F) %>%
  map(~lm(.$Nauplii ~ .$Time, na.action=na.exclude)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1) # Calculates the floquent exponent for each pho values for each simulation. V1 column is the floquent expoent and left column is pho.simulation
# The map_dbl responds to the coefficent that is taken.

nauplii.slope <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$Nauplii ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2)

copepodids.inter <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$Copepodids ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1)

copepodids.slope <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$Copepodids ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2)

chalimi.inter <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$`Chalimi and pre-adults` ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1)

chalimi.slope <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$`Chalimi and pre-adults` ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2)

adult.inter <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$`Adult Females` ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(1)

adult.slope <- simulation.data.frame %>%
  split(list(.$`Pho Value`, .$`Simulation Number`), lex.order = F) %>%
  map(~lm(.$`Adult Females` ~ .$Time)) %>%
  map(summary) %>%
  map(coef) %>%
  map_dbl(2)

# Combining all the datasets together to ger the floquet exponets
# Combing adult abundance
adult.floquet <- as.data.frame(cbind(adult.inter, adult.slope))
adult.names <- data.frame(row.names(adult.floquet))
adult.names <- adult.names %>%
  separate(row.names.adult.floquet., c("pho_ones", "tenths", "simulation_number")
           )
adult.names$simulation_number[seq(11,N*11, 11)] <- seq(1, N)
adult.names$tenths[seq(11, N*11, 11)] <- 0
adult.pho <- paste(adult.names$pho_ones, adult.names$tenths, sep=".")
adult.floquet <- cbind(adult.pho, adult.names$simulation_number, adult.floquet)
names(adult.floquet) <- c("Pho_value", "Sim. #", "Intercept", "Floquet")

# Looking at the head of each data frame
head(adult.floquet)

# Plotting Floquet Exponet ------------------------------------------------------
# Getting slope for the nauplii graph 

coef.adult <- coef(lm(Floquet ~ rev(factor(Pho_value)), data=adult.reg))
ggplot(adult.floquet, aes(x=factor(adult.floquet$Pho_value), y=adult.floquet$Floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Salinity Variation", y="Floquet Exponent", 
       title= paste("Adult Floquet:", 
                    "Intercept =", signif(coef.adult[[1]], 2), 
                    ", Slope =", signif(coef.adult[[2]], 2))) +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  geom_abline(intercept = coef.adult[1], slope = coef.adult[2], color = "red") +
  scale_x_discrete(limits = rev(levels(adult.floquet$Pho_value)))



# Exporting df to csv !!!!!!! Needs to be changed each time!!!!!!!!!
write.csv(simulation.data.frame,"/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/sim_data_100", row.names = F)
write.csv(nauplii.floquet, "/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/nauplii_floquet_100", row.names = F)
write.csv(copepodids.floquet, "/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/cop_flo_100", row.names = F)
write.csv(chalimi.floquet, "/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/chalimi_flo_100", row.names = F)
write.csv(adult.floquet, "/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/model_figures/model_figure_data/adult_flo_100", row.names = F)


# Running some analysis 
nauplii.floquet$Pho_value <- as.factor(nauplii.floquet$Pho_value)
copepodids.floquet$Pho_value <- as.factor(copepodids.floquet$Pho_value)
chalimi.floquet$Pho_value <- as.factor(chalimi.floquet$Pho_value)
adult.floquet$Pho_value <- as.factor(adult.floquet$Pho_value)

# Mean of floquet for each life stage 
nauplii.flo.mean <- with(nauplii.floquet, tapply(nauplii.floquet$Floquet, nauplii.floquet$Pho_value, mean)) # Prints the mean of each pho value 
cop.flo.mean <- with(copepodids.floquet, tapply(copepodids.floquet$Floquet, copepodids.floquet$Pho_value, mean))
chalimi.flo.mean <- with(chalimi.floquet, tapply(chalimi.floquet$Floquet, chalimi.floquet$Pho_value, mean))
adult.flo.mean <- with(adult.floquet, tapply(adult.floquet$Floquet, adult.floquet$Pho_value, mean))

# Print the means 
nauplii.flo.mean
cop.flo.mean
chalimi.flo.mean
adult.flo.mean











