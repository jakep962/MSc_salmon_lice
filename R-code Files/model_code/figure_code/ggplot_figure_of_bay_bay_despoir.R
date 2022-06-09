# Producing a map of Bay d'Espoir
rm(list=ls())

setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/first_chapter/winter_2020/model_code/data_files/baydespoir.nosync/")

#install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel",
                   # "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

# Loading Bay d'Espoir and data ----
baydespoir_prof <- read.csv("baydespoir_prof.csv")



# Packages needed
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library("ggspatial")
library(patchwork)
library(cowplot)
library(bbmle)
bay.sub <- subset(baydespoir_prof, baydespoir_prof$PSAL>=0)

# Plot of Bay d'Espoir
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

bay.map <- 
  ggplot(data = world) +
        geom_sf() +
        annotation_scale(location = "bl", width_hint = 0.5) +
        annotation_north_arrow(location = "bl", which_north = "true",
                               pad_x = unit(3.5, "in"), pad_y = unit(4, "in"),
                               style = north_arrow_fancy_orienteering) +
        annotate(geom = "text", x = -56.07, y = 47.64, label = "Bay d'Espoir",
                 fontface = "bold", color = "black", size = 4) +
        coord_sf(xlim = c(-58, -54), ylim = c(46, 49), expand = FALSE) +
        theme(
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)
          ) +
        labs(x="", y="")
bay.map

# Changing the time 

# The date format that is needed for the chron data 
monthly.dates <- dates(paste(bay.sub$OBS_MONTH, bay.sub$OBS_DAY, 
                             bay.sub$OBS_YEAR,
                             sep="/")) 

# The date format that ggplot used 
bay.sub$ggplot.data <- as.Date(paste(bay.sub$OBS_YEAR, bay.sub$OBS_MONTH, 
                                     bay.sub$OBS_DAY,
                                        sep="-")) 

bay.sub$chron_date <-chron(dates=monthly.dates, 
                              origin. = c(month = 01,day = 01,
                                          year = 1956)) # setting into chron
bay.sub <- arrange(bay.sub, chron_date) # ordering dates 

o <- as.Date("1956-01-01") # The origin or first date
bay.sub$days_since_origin <- (as.numeric(as.Date(bay.sub$chron_date) - o)) #

bay.sub.test <- subset(bay.sub, bay.sub$OBS_YEAR >= 1994 & bay.sub$OBS_YEAR <=2003)
# Plotting the salinity data 
  ggplot(bay.sub, aes(x=OBS_MONTH, y=PSAL)) +
  geom_jitter(aes(colour=factor(OBS_YEAR)))  + 
  theme_classic() +
  xlab("") +
    scale_x_discrete(name ="Month Sampled", 
                     limits=c("1","2","3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12")) +
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


# Information needed for salinity plots
# Loading Bay D'Esspoir dataset
baydespoir_prof <- read.csv("baydespoir_prof.csv")

# Using a lower dx will produce a smoother line for the curve. When not using to generate the plot, it should
# be set to 0.1
#  Asymetric Laplace distrubution -----------------------------------------------
# Using the Asymmetric Laplace distrubution to produce a cummulative distrubution
# Changing the bin and x-vectors
dx <- 0.25
sub.edge <- seq(15, 40, dx)
sub.mid <- seq(15+dx/2, 40-dx/2, dx)
sub.x <- sub.mid

# Subsetting to only take values that are greater then 15
sub.bay <- subset(baydespoir_prof, PSAL >= 15)
sub.bay <- na.omit(sub.bay$PSAL)

# Plotting the subset distrubution
sub.hist <- hist(sub.bay, sub.edge, xlab="Salinity Bins (0.5 psu)", 
                 main="", cex.lab= 1.5, cex.axis=1.5) # Plotting the histogram when subsetted from 15-40
sub.freq <- sub.hist$counts
sub.freq <- sub.freq/sum(sub.freq)# Converting the count data into frequency

# Find the m value
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1
pho0 <- 1.9137
k0 <- 1.223080

sub.asy.laplace <- function(sub.m, pho, k){
  d <- (pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
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
        
        c(rep(4.00, 50), rep(3.50, 50), rep(3.00, 50), rep(2.50, 50), rep(2.00, 50), rep(1.82, 50), rep(1.44, 50),
          rep(1.06, 50), rep(0.67, 50), rep(0.48, 50), rep(0.29, 50))
  ))

names(var.salinity.dist.df) <- c("Salinity", "Probability", "Scale Parameters")
var.salinity.dist.df$`Scale Parameters`<- round(var.salinity.dist.df$`Scale Parameters`, digits=2)

var.sal.plot <-ggplot(var.salinity.dist.df, aes()) +
    geom_line(aes(x=Salinity, y=Probability, 
                  colour=factor(`Scale Parameters`))) + 
    scale_colour_manual(values = rev(c("black", "brown", "blue", "purple", "pink",
                                       "red", "orange", "yellow", "cyan", "green", 
                                       "darkgreen")), 
                        name = bquote('Scale in\n daily salinity '~lambda)) +
    xlim(28, 36) + 
    ylim(0, 1) +
    labs(x = "Salinity (PSU)",
         y = "Probability") +
  theme_classic() +
  theme(axis.text.x=element_text(size= 14), 
        axis.text.y=element_text(size= 14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size = 14)) 


# Producing the single histogram 
hist.sal.plot <- as.data.frame(cbind(sub.x, sub.bay, 
                                     sub.edge, sub.asy.laplace(sub.m, 1.91, 1.22), 
                                     rep(1.8217, 50)))

bay.sal <- 
ggplot(hist.sal.plot, aes(x=sub.bay), colour=V5) +
  geom_histogram(aes(y=..density..), colour="black", fill="gray") +
  geom_line(aes(x=sub.x, y=V4), col="red", show.legend = TRUE) +
  xlim(26,36) +
  theme_classic() +
  theme(axis.text.x=element_text(size= 15), 
        axis.text.y=element_text(size= 15), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size = 15)) + 
  labs(x="Salinity (psu)", y="Probability") +
  scale_colour_manual(values = "red", 
                      name = "Scale Parameter")

vertical.plots <- plot_grid(bay.sal, var.sal.plot, labels = c("(b)", "(c)"),
                            ncol = 1, nrow=2,  align = 'v')
plot_grid(bay.map, vertical.plots, labels= c("(a)", ""))




ggplot(baydespoir_prof, aes(x=OBS_MONTH, y=PSAL)) +
  geom_jitter() +
  labs(x = "Month",
       y = "Salinity (PSU)") +
  theme_classic() +
  theme(axis.text.x=element_text(size= 14), 
        axis.text.y=element_text(size= 14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size = 14)) 


