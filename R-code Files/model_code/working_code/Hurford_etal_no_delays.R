# Modified from Figure3_Dynamics.R archived on figshare. The goal is to reformulate
# this model to be solved using Euler's method and without any maturation delays.

rm(list=ls())

# Set the path to your directory

setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/winter_2020/model_code/data_files /amy_regional//")

# Download the files below from:
# https://figshare.com/articles/Regional_climate_affects_salmon_lice_dynamics_stage_structure_and_management/7824152

load("ParamsFuncs.RData")
load("farmNL.RData")

library(dplyr)

#plotting plot(saltime, salNL, ylim = c(0,35)) will recover figure 1(g) from Hurford et al. 2019

plot(saltime, salNL, ylim = c(0,35))
abline(h=mean(salmean), col="red", lwd=2)

# The file ParamsFuncs.RData is the output of Figure_S1_S2_ParamFunc.R and all the data underlying
# the parameterizations is described in more detail there.

# iota1 is the non-seasonal estimate of iota and for the BC site iota is parameterized from
# it appears that salinity is non-seasonal. Set iota = iota1.

iota <- iota1

# Set number of fish on farm equal to

f <- 500000

# Under this reformulation, the assumptions around the distribution of maturation times at any fixed temperature
# have changed, however, the data, at a fixed temperature, are insufficient to know if maturation times more
# closely follow a geometric distribution (as assumed here) or a delay formulation (as per the previous formulation).

# Note that gammaP(12) = 0.38 per day. This means that the distribution of times for maturation from nauplii to
# copepodid is geometrically distributed with the mean days to maturation equal to 1/0.38 = 2.6 days at 12C.


# Length of simulation

T = 15*365 # How many days there will be within the model
t.v = seq(1, 15*365, 1) # A vector with all values of T

# Temperture as a function of time

Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
  return(Temp)
}

# Plot showing temperature through time 

plot(t.v, Temp(t.v))

# Preallocate
S = rep(0,T)
P = rep(0,T+1)
C = P
A = P
I = P
# Initial values
A[1]=1

for(t in seq(1,T)){
  # Salinity - if salinity is changed to be a function, be careful that the same value of salinity is used
  # for the same t. Note that S[t] occurs in each of P, I, C, and A, and you don't want different realizations
  # of this value for the same day.
  
  S[t] = rnorm(1,30,5)
  
  # Stages
  
  P[t+1] =  P[t] + eta(Temp(t))*epsilon(Temp(t))*A[t]*V(Temp(t),S[t]) - mup(S[t])*P[t] - gammaP(Temp(t))*P[t]
  I[t+1] =  I[t] + gammaP(Temp(t))*P[t] - iota*f*I[t] - mui(S[t])*I[t]
  C[t+1] = C[t] + iota*f*I[t] - gammaC(Temp(t))*C[t] - muc(S[t])*C[t]
  A[t+1] = A[t] + gammaC(Temp(t))*C[t] - mua(S[t])*A[t]
  
}

# Building the dataframe 

lice.fish <- data.frame((seq(1,T+1)/365), A/f) # combining data into a data frame
names(lice.fish) <- c("Time_(years)", "A_per_fish")
year.data <- NULL
year.data <- lice.fish[ lice.fish$`Time_(years)` %in% c(1:15),] # Fish count at the start of each year


par(mfrow=c(2,1))
plot(seq(1,T+1)/365,A/f, typ="l", xlab = "time in years", ylab = "A per fish")
plot(seq(1,T)/365, S, xlab = "time in years", ylab = "salinity", pch = 20) # See need to look at more seasonal patterns
# in the salinity data

#Finding the Floquet exponent,(could be another name)

par(mfrow=c(1,1))
plot(year.data$`Time_(years)`, log(year.data$A_per_fish))

#Build linear model to fit to the points

flow.mew <- lm(log(year.data$A_per_fish) ~ year.data$`Time_(years)`, data = year.data) #Fitted line to data, seems to
flow.mew1 <- lm(year.data$A_per_fish ~ year.data$`Time_(years)`, data = year.data) #Fitted line to data, seems to

abline(flow.mew, col="red")
summary(flow.mew)
summary(flow.mew1)

# For the Floquet exponent would a simple version of it just be the slope of the model? The more times that you reran the 
# model and averaged the slope, would it be more accurate? 


