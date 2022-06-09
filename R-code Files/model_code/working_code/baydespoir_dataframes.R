# R code to open and view the Bay D'espoir files. "baydespoir" file contains a zip file with all three files compressed, 
# cruise, prof, and curf .csv files. 

rm(list=ls())

#Columns:
# CR_NUMBER: mission identifier
# DATA_TYPE
# BA: Data from expendable bathythermograph (XBT) formatted for real-time transmission on the World Meteorological Organization Global Telecommunication System
# *BO: Data from bottle
# *CD: Data from CTD (downcast)
# *MB: Data from mechanical bathythermograph *TE: Data formatted for real-time transmission on the World Meteorological Organization Global Telecommunication System
# *CD: Data from CTD (Towed mode)
# *TK: Data from surface underway thermosalinograph formatted for real-time transmission on the World Meteorological Organization Global Telecommunication System
# -D_P_CODE: indicates whether value reported in DEPTH_PRESS is Depth (m) or Pressure (dbar)
# -Q_* fields: quality flag qualifying preceding value: 1-good, 3-suspicious, 4-bad
# -OBS_YEAR, OBS_MONTH, OBS_DAY, OBS_TIME: observation time in UTC
# -SOURCE_ID: organization who provided us the data
# *BIO = Bedford Institute of Oceanography, Dartmouth NS (DFO)
# *NAFC = Northwest Atlantic Fisheries Centre, St. John's NL (DFO)
# *MEMU = Memorial University, St. John's NL
# *CANA = (Royal) Canadian Navy
# *CWOW = Data inserted in the World Meteorological Organization Global Telecommunication System by Canada
# -TEMP: water temperature at depth (ºC)
# -PSAL: practical salinity at depth
# -USAL: salinity at depth (uknown method)
# -SSAL : salinity at depth (g/kg)
# -PLAT: platform type (ship onboard which data were collected; see ICES list: http://vocab.ices.dk/?sortby=description&ref=315)
# -PRO1: Project name
# -CSC1: Chief scientist name
# -other 4 character fields: see https://www.nodc.noaa.gov/GTSPP/document/codetbls/gtsppcode.html for definition and units

# Set wd 

setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/winter_2020/model_code/data_files /baydespoir/")

# Import the data 

cruise <- read.csv("baydespoir_cruise.csv", header = T)
prof <- read.csv("baydespoir_prof.csv", header = T) # Has the most inforamtion on salinity
surf <- read.csv("baydespoir_surf.csv", header = T)

# Looking for patterns within the salinity data 
plot(prof$OBS_MONTH, prof$PSAL) # Pretty even, some drops in salinity around 10 (October)
plot(prof$OBS_YEAR, prof$PSAL) # Data isn;t taken every year and has some gaps. Lower salinity
# levels taken ~1997
plot(prof$OBS_TIME, prof$PSAL) # Missing a lot of inforamtion before 10 am. After this, the data
# becomes a more homogenous

# Seeing if there is a pattern with when the observation was taken (time in 24H)
plot(prof$OBS_YEAR, prof$OBS_TIME) 
plot(prof$OBS_MONTH, prof$OBS_TIME)


# Subset data; focus from 10 am on forward 
prof_subset <- subset(prof, PSAL >0 & OBS_TIME > 1000, select= OBS_YEAR:PSAL) # Removing all

plot(prof_subset$OBS_MONTH, prof_subset$PSAL)
plot(prof_subset$OBS_YEAR, prof_subset$PSAL)
hist(prof_subset$OBS_MONTH)  # Not equal distrubution of observations across months 
hist(prof_subset$OBS_YEAR)   # 90% of the data comes from 1996

# Looks like most of the data is comig from 1997. I thused subset the data to only contain 1997. The following
# plots aim to look at the distrubution of the data. 

prof_1997 <- subset(prof, OBS_YEAR == 1997, select=OBS_MONTH:PSAL)
hist(prof_1997$PSAL) # Most of the data is between 30:35 PSU 
hist(prof_1997$OBS_MONTH) # Data is only recorded at four months 

#Looking at how the mean, mode, and medium change between the raw data and the subset

mean(prof$PSAL, na.rm = T) # 32.02758
mean(prof_1997$PSAL, na.rm = T) # 31.80257

median(prof$PSAL, na.rm = T) # 32.177
median(prof_1997$PSAL, na.rm = T) # 32.117

# There is little change between the raw and subet data, likely due to large weight of the data at the upper 
# tail. 






