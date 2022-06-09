# Building the Bay D'Espoir plot
rm(list=ls())

# Packages
library(bbmle)
library(ggplot2)

# Setting wd and loading additional Rdata files
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/first_chapter/winter_2020/model_code/data_files/baydespoir.nosync//")

# Loading Bay D'Esspoir dataset
baydespoir_prof <- read.csv("baydespoir_prof.csv")

# Using a lower dx will produce a smoother line for the curve. When not using to generate the plot, it should
# be set to 0.1

#  Asymetric Laplace distrubution -----------------------------------------------
# Using the Asymmetric Laplace distrubution to produce a cummulative distrubution
# Changing the bin and x-vectors
dx <- 0.5
sub.edge <- seq(15, 40, 0.5)
sub.mid <- seq(15+dx/2, 40-dx/2, 0.5)
sub.x <- sub.mid

# Subsetting to only take values that are greater then 15
sub.bay <- subset(baydespoir_prof, PSAL >= 15)
sub.bay <- na.omit(sub.bay$PSAL)

# Plotting the subset distrubution
sub.hist <- hist(sub.bay, sub.edge) # Plotting the histogram when subsetted from 15-40
sub.freq <- sub.hist$counts
sub.freq <- sub.freq/sum(sub.freq) # Converting the count data into frequency


# Find the m value
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1

pho0 <- 1.9137
k0 <- 1.223080
# Need to update the asymmetric laplace distrubution so that it pulls the right
# x value

sub.asy.laplace <- function(sub.m, pho, k){
  d <- (pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
} # Same distrubution but x was changed to sub.x. Sub.x is the x-vector for the subset data

sub.dx <-
# Produce the asymmetric laplce with the dx
sub.asy.laplace.dx <- function(sub.m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
}

sum(sub.asy.laplace.dx(sub.m, pho0, k0))


sub.negLL = function(pho, k){
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
} # Only want to fit to pho and k so only include them within the code

sub.mle.asy.laplace.Ex <- mle2(sub.negLL,
                               start=list(pho=1, k=3), ) # Similary, only want to fit to pho and k

a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]

plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", col="blue", main="dx=0.5") # Plotting the
points(sub.x, sub.freq, col="red", pch=20, ylim=c(0,0.2)) # Plotting the empircal data

# fitted asymmetric distrubution
print(a.est)

# Initial fitted values


# The probability distrubution is described by
plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", ylim=c(0,1.1)) # Probability
lines(sub.x, cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2])), col="red") # Cumlative
points(sub.x, sub.freq, col="blue", pch=20)

# Combing the distrubutions into one distrubution
sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]),
                                  cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2]))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
head(sub.sal.df)


hist.sal.plot <- as.data.frame(cbind(sub.x, sub.bay/sum(sub.bay), sub.bay,
                                     sub.edge,sub.asy.laplace(sub.m, 1.91, 1.22)))



bay.sal <-
  ggplot(hist.sal.plot, aes(x=sub.bay)) +
    geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_line(aes(x=sub.x, y=V5), colour="red") +
    xlim(26,36) +
    labs(x="Salinity (PSU)", y="Probability") +
    theme_classic()


# # Building the plot -------------
# hist(sub.bay, sub.edge, prob=T, xlim=c(25,37), ylim=c(0,1), xlab="Salinity (PSU)", ylab="Probability", main="")
# lines(sub.x, sub.asy.laplace(sub.m, 1.91, 1.22), col="red")
# lines(sub.x, sub.asy.laplace.dx(sub.m, 1.91, 1.22), col="blue")

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

 
var.sal.plot <- ggplot(var.salinity.dist.df, aes(x=Salinity, y=Probability, colour=factor(`Scale Parameters`))) +
                    geom_line() + 
                    scale_colour_manual(values = rev(c("black", "brown", "blue", "purple", "pink",
                                                   "red", "orange", "yellow", "cyan", "green", 
                                                   "darkgreen")), 
                                        name = "Scale Parameter") +
                    xlim(26, 35) + 
                    labs(x = "Salinity (PSU)",
                         y = "Probability") +
                    theme_classic()
