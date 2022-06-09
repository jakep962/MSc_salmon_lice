#### Code for the asymetric Laplace distrubution mle2 ####
rm(list=ls())

# Packages 
library(bbmle)
library(ggplot2)

# Fitting a Laplace distrubution to fake data 
dx <- 0.1
xedge <- seq(0, 40, dx)
xmid <- seq(0+dx/2, 40-dx/2, dx)
x <- xmid

asy.laplace <- function(m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((x-m)*pho*(sign(x-m))*k^(sign(x-m))))
}


plot(x, asy.laplace(33, 1, 3), type="l", ylim=c(0,0.04))
sum(asy.laplace(33, 1,3)) # Want the function to sum to 1. 

y.obs <- asy.laplace(33,1,3)*runif(length(x)) # Produce simmulated data
y.obs <- y.obs/sum(y.obs) # This changes y.obs from a number to a frequency
sum(y.obs) # Should now sum to 1

# Plotting the unfitted Asy.Laplace ditsrubution populated with runif points
plot(x, asy.laplace(33, 1, 3), col="black", type="l", ylim=c(0, 0.04))
points(x, y.obs, col="red") # Simulated data plotted against unfitted distrubution

# Calculating the location parmater
ran.df <- as.data.frame(cbind(x, y.obs))
ran.max <- ran.df[which.max(ran.df$y.obs),]
ran.max # Finding the m parameter

# Fitting the data distrubution
m <- ran.max[1,1]
sigma <-1 
negLL = function(pho, k){
  
  -sum(dnorm(y.obs, mean=asy.laplace(m, pho, k), sd=sigma, log=T))

}

mle.asym.laplace.Ex <- mle2(negLL, 
                            start=list(pho=1, k=3))

a.est <- unname(coef(mle.asym.laplace.Ex))
lines(x, asy.laplace(m, a.est[1], a.est[2]), col = "blue")
print(a.est) # The parameter values for pho and k. The m parameter is the max frequency and sd is set to 1.


#### Importing the Bay d'Espoir data and trying to fit to that ####
# Setting wd and loading additional Rdata files
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/data_files/amy_regional/")
load("ParamsFuncs.RData")
load("farmNL.RData")

# X vectors for binning the data
dx <- 0.1
xedge <- seq(0, 40, dx)
xmid <- seq(0+dx/2, 40-dx/2, dx)
x <- xmid

# Loading Bay D'Esspoir dataset
baydespoir_prof <- read.csv("../baydespoir/baydespoir_prof.csv") 
y.bay <- baydespoir_prof$PSAL
y.bay <- na.omit(y.bay)
bay.hist <- hist(y.bay, xedge) # Uses xedge to bin the data into catergories. This is the whole data
y.freq <- bay.hist$counts
y.freq <- y.freq/sum(y.freq) # Turns it from a number to a freq that sums to 1
sum(y.freq) # Both the freq and unfitted Asy.Laplace should sum to 1

# The asymetric function
asy.laplace <- function(m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((x-m)*pho*(sign(x-m))*k^(sign(x-m))))
}

# Plotting the distrubution
y.val <- asy.laplace(33, 1, 3)
plot(xmid, y.freq, col="red", pch=20)

# Fitting the data 
negLL = function(m, pho, k, sigma){
  -sum(dnorm(y.freq, mean=asy.laplace(m, pho, k), sd=sigma, log=TRUE))
}

mle.asym.laplace.Ex <- mle2(negLL, start=list(m=30, pho=1, k=1, sigma=1))

a.est = unname(coef(mle.asym.laplace.Ex))[-4]
lines(x, asy.laplace(a.est[1], a.est[2], a.est[3]), col = "blue")
print(a.est) # Produces negative numbers but this is an issue. 

# This maybe a result of the low salinity values at the left tail. 

# Subset the dataframe to focus more on the abundance of points
# Only going to be taking values that are greater then 15

#### Subsetting and refitting the data ####

# Changing the bin and x-vectors
dx <- 0.5
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Plot the unchanged distrubution
plot(bay.hist)

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
plot(sub.x, sub.freq) # still a break left tail, could bring it in even closer (20-40)?

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

plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="blue") # Plotting the 
points(sub.x, sub.freq, col="red", pch=20, ylim=c(0,0.2)) # Plotting the empircal data

# fitted asymmetric distrubution. 
print(a.est)

# Coefficents are pho: 3.7387943  and k: 0.8711824

# How/where the parameters are being fitted 
# m: m is the location parameter and thus where the peak in the distrubution is. Count 
# data for the frequency is calculated and the peak is taken as the parameter. This 
# can be found on line 119-122.
# The pho and k parameters are fitted by the mle2 package. 
# The sd is set to 1. 


#### Changing dx to 0.5 ####

# Changing the bin and x-vectors
dx <- 0.5
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Plot the unchanged distrubution
plot(bay.hist)

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
plot(sub.x, sub.freq) # still a break left tail, could bring it in even closer (20-40)?

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

plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="blue", main="dx=0.5") # Plotting the 
points(sub.x, sub.freq, col="red", pch=20, ylim=c(0,0.2)) # Plotting the empircal data

# fitted asymmetric distrubution. 
print(a.est)

# Coefficents are pho: 1.912716  and k: 1.223080

# How/where the parameters are being fitted 
# m: m is the location parameter and thus where the peak in the distrubution is. Count 
# data for the frequency is calculated and the peak is taken as the parameter. This 
# can be found on line 119-122.
# The pho and k parameters are fitted by the mle2 package. 
# The sd is set to 1. 

#### Changing dx to 1 ####


# Changing the bin and x-vectors
dx <- 1
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Plot the unchanged distrubution
plot(bay.hist)

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
plot(sub.x, sub.freq) # still a break left tail, could bring it in even closer (20-40)?

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

plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="blue", main="dx=1") # Plotting the 
points(sub.x, sub.freq, col="red", pch=20, ylim=c(0,0.2)) # Plotting the empircal data

# fitted asymmetric distrubution. 
print(a.est)

# Coefficents are pho: 1.912716  and k: 1.223080

# How/where the parameters are being fitted 
# m: m is the location parameter and thus where the peak in the distrubution is. Count 
# data for the frequency is calculated and the peak is taken as the parameter. This 
# can be found on line 119-122.
# The pho and k parameters are fitted by the mle2 package. 
# The sd is set to 1. 


#### Using dx=0.5 and the subsetted data to see how chaning pho effects the distrubution ####

## Re running code for dx set to 0.1 and subset --------------------------
dx <- 0.5
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Plot the unchanged distrubution
plot(bay.hist)

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
plot(sub.x, sub.freq) # still a break left tail, could bring it in even closer (20-40)?

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

plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="blue") # Plotting the 
points(sub.x, sub.freq, col="red", pch=20, ylim=c(0,0.2)) # Plotting the empircal data

# fitted asymmetric distrubution. 
print(a.est)

## Changing the values of pho as a proxy for variance 
plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="red") # pho is set to 1.91


# Creating a ggplot of all the different pho values 
pho.values <- rev(seq(0, a.est[1], a.est[1]/10)) # A vector of pho value that spread

# Plotti
plot(sub.x, sub.asy.laplace(sub.m, pho.values[1], a.est[2]), type="l", col="red")
points(sub.x, sub.asy.laplace(sub.m, pho.values[2], a.est[2]), type="l", col="blue")
points(sub.x, sub.asy.laplace(sub.m, pho.values[3], a.est[2]), type="l", col="green")
points(sub.x, sub.asy.laplace(sub.m, pho.values[4], a.est[2]), type="l", col="yellow")
points(sub.x, sub.asy.laplace(sub.m, pho.values[5], a.est[2]), type="l", col="pink")
points(sub.x, sub.asy.laplace(sub.m, pho.values[6], a.est[2]), type="l", col="purple")
points(sub.x, sub.asy.laplace(sub.m, pho.values[7], a.est[2]), type="l", col="black")
points(sub.x, sub.asy.laplace(sub.m, pho.values[8], a.est[2]), type="l", col="orange")
points(sub.x, sub.asy.laplace(sub.m, pho.values[9], a.est[2]), type="l", col="gray")
points(sub.x, sub.asy.laplace(sub.m, pho.values[10], a.est[2]), type="l", col="brown")
legend(15, 0.43, legend=c("Pho Values", "1.913", "1.7", "1.5", "1.3", "1.14", "0.9", "0.76", "0.57", "038", 
                          "0.191"),
       col=c("white", "red", "blue", "green", "yellow", "pink", "purple", "black", "orange", "gray", "brown"), 
       lty=1:2, cex=0.8)


