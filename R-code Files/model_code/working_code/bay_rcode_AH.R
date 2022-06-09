# Code for the asymetric Laplace distrubution mle2 ------------------------
rm(list=ls())

# Packages 
library(bbmle)

# Fitting a Laplace distrubution 
# AH: Define an x-vector to bin data across later 
dx <- 0.1
# AH: Edge points and mid-points for the vector to bin on
xedge <- seq(0,40,dx)
xedge2 <- seq(5, 40, dx)
xmid <- seq(5+dx/2,40-dx/2,dx)
x<-xmid

asy.laplace <- function(m, pho, k){
  # AH: added the dx multiplier
  d <- dx*(pho/(k+(1/k)))*exp(-((x-m)*pho*(sign(x-m))*k^(sign(x-m))))
}

plot(x, asy.laplace(33, 1, 3), typ = "l", ylim=c(0, 0.040))
sum(asy.laplace(33, 1, 3)) # Does sum to 1
# AH: I want the data and the function to both sum to 1.
# AH: Here the data is changed from numbers to frequency
y.obs <- asy.laplace(33, 1, 3)*runif(length(x))
y.obs <- y.obs/sum(y.obs)
points(x,y.obs, col = "red")
lines(x, asy.laplace(33, 1, 3))
# Above, I'm plotting the un-fitted laplace with the data.


negLL = function(m, pho, k, sigma){
  -sum(dnorm(y.obs, mean=asy.laplace(m, pho, k), sd=sigma, log=T))
}

mle.asym.laplace.Ex <- mle2(negLL, 
                            start=list(m=33, pho=3, k=1, sigma=1))

a.est <- unname(coef(mle.asym.laplace.Ex)[-4])
# AH: The fit here looks much better
lines(x, asy.laplace(a.est[1], a.est[2], a.est[3]), col = "blue", lwd=2)
# Above we add the fitted laplace to the previos plot, where the fitted Laplace
# is the thick black line.

print(a.est)
# plot(asy.laplace(a.est))
# asy.laplace(a.est)
# 
# asy.laplace(a.est[1], a.est[2], a.est[3])
#  

# plot(x, asy.laplace(a.est[1], a.est[2], a.est[3]))


# Importing the Bay d'Espoir data and trying to fit to that --------------
# Setting wd and loading additional Rdata files
#setwd("~/Desktop/Work/Students/MSc/Prosser/Code/") ------
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/data_files/amy_regional/")
load("ParamsFuncs.RData")
load("farmNL.RData")

# Loading Bay D'Espoir dataset
#baydespoir_prof <- read.csv("/Users/amyhurford/Desktop/Work/Students/MSc/Prosser/Code/baydespoir_prof.csv") 
baydespoir_prof <- read.csv("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/winter_2020/model_code/data_files/baydespoir/baydespoir_prof.csv")
hist(baydespoir_prof$PSAL, prob=T, breaks= 18, xlim=c(0,40), xlab="Salinity", 
     ylab="Probability", main="")

y.bay <- na.omit(baydespoir_prof$PSAL)
y.bay2 <- subset(baydespoir_prof, baydespoir_prof$PSAL > 5)
y.bay2 <- na.omit(y.bay2$PSAL)


# AH: I use the hist() function to get the frequencies associated with the
# salinity values which go on the x-axis.

y.hist <- hist(y.bay, breaks=xedge)
y.hist2 <- hist(y.bay2, breaks = xedge2)
y.freq <- y.hist$counts
y.freq <- y.freq/sum(y.freq)
y.freq2 <- y.hist2$counts
y.freq2 <- y.freq2/sum(y.freq2)
sum(y.freq2)

# AH: again, both the data and the unfitted Laplace values should sum to 1.
# plot the unfitted laplace and the data
y.val <- asy.laplace(33,3,1)
plot(xmid, y.freq2, col = "red", pch = 20)
lines(xmid,y.val,pch=20, ylim = c(0, max(y.freq,y.val)))

# Check everything sums to 1
sum(y.freq2)
sum(y.val)

negLL = function(m, pho, k, sigma){
  
  -sum(dnorm(y.freq2, mean=asy.laplace(m, pho, k), sd=sigma, log=TRUE))
  
}

mle.asym.laplace.Ex <- mle2(negLL, start=list(m=30, pho=3, k=1, sigma=1))

a.est = unname(coef(mle.asym.laplace.Ex)) 

lines(x, asy.laplace(a.est[1], a.est[2], a.est[3]), col = "blue")

# Fitted curve plotted by itself
plot(xmid,y.freq, col = "red", pch = 20)
lines(x, asy.laplace(a.est[1], a.est[2], a.est[3]), type="l")

# AH: The first fit that I did of this, the estimated values of m and pho were
# negative. I don't think these values should be negative, so please check this with
# the fitting. If maybe necessasry to write the functions so that it is impossible
# to get negative parameter estiamtes back.

# The current fit, I doubt it is optimal. Note the poor fit to the highest frequency
# salinities. The problem is both fitting those peak values, while also fitting the long
# tail on the low salinity values. These low salinity values are relatively rare. I
# think for use, we maybe should forget about those values and focus on the region
# with the higher probabilities. I think we should probability define a min salinity
# cutoff and only fit to data higher than that value.
 

