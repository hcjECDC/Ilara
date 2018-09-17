rm(list=ls())
dir <- "P:/Measles/Greek data/"
setwd(dir)
library(plyr)
library(ggmap)
library(devtools)
#install_github("ropengov/eurostat")
library(eurostat)
library(sf)
library(sp)
library(tmap)
library(geosphere)
library(lubridate)
library(reshape2)
library(stringr)
library(kSamples)
library(tidyverse)
library(purrrlyr)
library(dplyr) # this needs to be loaded last due to conflicts with plyr
options(digits=8)

## In matrices, i (row) is infected by j (column)

# Create custom colour scale
myColours <- c("26 107 133", "241 214 118", "168 45 23")
ECDCcol <- sapply(strsplit(myColours, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))  # convert to hexadecimal

# Measles parameters
mean1 <- 11.7 # mean of normal distribution for serial interval (days) ##Vink et al. 2014
sd1 <- 0.92 # standard deviation of normal distribution for serial interval (days)  ##Vink et al. 2014

# Read in line list
data0 <- read_csv("workingGreekData.csv", col_names=TRUE)
MMRData <- read_csv("MMR2.csv", col_names=TRUE, col_types = cols(resPostcode = col_character()))  # need to specify residencePostcode as character here because two entries include letters, not just numbers. Otherwise parser would assume it was a column of integers.
popData <- read_csv("pop2015.csv", col_names=TRUE) 

latLong <- read_csv("LatLong.csv", col_names=TRUE) # With 'workingGreekData.csv the home town co-ordinates are included

data0$onsetDate <- as.Date(data0$symptomOnset, "%d/%m/%Y")  # convert to date format
data0$notifyDate <- as.Date(data0$dateNotification, "%d/%m/%Y")  # convert to date format
data0$hospitalDate <- as.Date(data0$dateHospitalisation, "%d/%m/%Y")  # convert to date format
data0$hosp1Lat <- as.numeric(data0$hosp1Lat)
data0$hosp1Lon <- as.numeric(data0$hosp1Lon)
data0$exp1Lat <- as.numeric(data0$exp1Lat)
data0$exp1Lon <- as.numeric(data0$exp1Lon)
data0$exp2Lat <- as.numeric(data0$exp2Lat)
data0$exp2Lon <- as.numeric(data0$exp2Lon)
data0$HCWHospLat <- as.numeric(data0$HCWHospLat)
data0$HCWHospLon <- as.numeric(data0$HCWHospLon)
data0$deathDate <- as.Date(data0$dateDeath, "%d/%m/%Y")  # convert to date format
data0$dischargeDate <- as.Date(data0$dateDischarge, "%d/%m/%Y")  # convert to date format
data0$symptomOnset <- NULL  # remove old variable
data0$dateNotification <- NULL  # remove old variable
data0$dateHospitalisation <- NULL  # remove old variable
data0$dateDeath <- NULL  # remove old variable
data0$dateDischarge <- NULL  # remove old variable
data0 <- data0[order(data0$onsetDate),]

write.csv(data0, file = "180509greekData.csv", row.names = FALSE)			

##################
## MMR coverage ##
##################

# Rearrange prescribing data from long to wide format (group all vaccine codes together (all MMR))
data2 <- MMRData %>% 
 mutate(purchaseDate = substr(purchaseDate, 1, 10))%>% # drop times
 select(-initialOrder, -vaccineCode) %>% # Need to drop this otherwise record number coerces R into keeping all rows
 arrange(id, purchaseDate) %>% # Ensure that dates are in correct order before inferrinstr(g dose number
 group_by(id) %>% 
 dplyr::mutate(dose = row_number()) %>% # Dose number for spreading data
 dplyr::mutate(doses = max(row_number())) %>% # Total number of doses for record
 unite("datePlace", c(purchaseDate, pharmPostcode, pharmDistrict)) %>%
 spread(dose, datePlace, fill="")
	
# Separate out date and place fields. Label columns
data3 <- data2 %>% select(id, dateOfBirth, gender, doses, resPostcode, resDistrict, resNUTS3)

for(i in (ncol(data3)+1):(ncol(data2))){   # NB. Start counting after constant columns
	name <- names(data2[,i])
	temp <- separate(data=data2, col=name, c(paste(name,"Date", sep=""),paste(name,"Pharm", sep=""),paste(name,"District", sep="")), "_", extra = "drop", fill = "right")
	data3 <- bind_cols(data3,temp[,i], temp[,i+1], temp[,i+2])
	}

#write.csv(data3, file = "180815MMRoutput.csv", row.names = FALSE)	

# Number of doses by NUTS-3 
oneDose <- data3 %>% filter(resNUTS3 != "#N/A" & doses==1) %>%
					arrange(resNUTS3) %>%
					group_by(resNUTS3) %>%
					summarise(Total = n()) %>%
					rename(NUTS3Code = resNUTS3) 
					
twoPlusDoses <- data3 %>% filter(resNUTS3 != "#N/A" & doses>1) %>%
						  arrange(resNUTS3) %>%
						  group_by(resNUTS3) %>%
					      summarise(Total = n())%>%
						  rename(NUTS3Code = resNUTS3) 
						  
pop2015 <- popData %>% arrange(NUTS3Code) %>%
					   group_by(NUTS3Code) %>%
					   summarise(Total = sum(popCohort15))
					   
MMRcoverage <- full_join(pop2015, oneDose, by="NUTS3Code")
MMRcoverage <- full_join(MMRcoverage, twoPlusDoses, by="NUTS3Code") %>%
			   rename(population = Total.x) %>%
			   rename(oneDose = Total.y) %>%
			   rename(twoPlusDoses = Total) %>%
			   mutate(coverage1 = oneDose/population) %>%
			   mutate(coverage2plus = twoPlusDoses/population) %>%
			   mutate(overallCoverage = coverage1 + coverage2plus) 
			   
# 2010 NUTS-3 code (Attiki) was subdivided into 7 codes in 2013 but we need to update population data. Average for now.
athensCodes <- c(paste("EL30", seq(1,7), sep=""))			   
expandAthens <- MMRcoverage %>% select(-NUTS3Code) %>%
								slice(rep(1, each=7)) %>%	          # copy Attiki entry 7 times
								mutate(NUTS3Code = athensCodes) %>%	  #	assign list of 7 Attiki codes 
								select(names(MMRcoverage))
								
										   
MMRcoverage <- 	bind_rows(MMRcoverage[-1,], expandAthens)	%>%	   
			    mutate(geo = NUTS3Code)   # to join with Eurostat shapefile 
							
							
							
			   
# Classes for coverage		
# 1 dose
lower <- 0.6
upper <- 1
sep = "-"
by <- 0.05
 
labs1 <- c("0-0.59", paste(seq(lower, upper - 0.05, by = by),
                 seq(lower + by - 0.01, upper, by = by),
                 sep = sep))
		   
MMRcoverage$classCoverage1 <- cut(MMRcoverage$overallCoverage, breaks = c(0, seq(lower,upper,by=0.05)), labels = labs1, include.lowest=FALSE, right=FALSE)
   
# 2+ doses
lower <- 0.3
upper <- 1
sep = "-"
by <- 0.1
 
labs2 <- c("0-0.29", paste(seq(lower, upper - 0.1, by = by),
                 seq(lower + by - 0.01, upper, by = by),
                 sep = sep))
		   
MMRcoverage$classCoverage2 <- cut(MMRcoverage$overallCoverage, breaks = c(0, seq(lower,upper,by=0.1)), labels = labs2, include.lowest=FALSE, right=FALSE)


#mapCoverage <- MMRcoverage %>% select(geo, classCoverage)
			 
###################
## Plotting maps ##			   
###################

# Plot Greece
mapGreece <- get_eurostat_geospatial(output_class = "sf", resolution = "10", nuts_level = 3)
mapGreece <- mapGreece %>% filter(CNTR_CODE == "EL")
mapGreece %>% select(geo) %>%
				plot()

# Plot MMR coverage				  
mapCoverage <- left_join(mapGreece, MMRcoverage, by="geo") %>%
			   filter(geo != "NA")


# Reading in shapefiles
#localDir <- paste(dir,"Shapefiles/LocalCommunities/", sep="")	
#setwd(localDir)		
#localMap <- st_read("local.shp")
#localMap %>% select(geometry) %>%
#				plot()
#setwd(dir)

muncDir <- paste(dir,"Shapefiles/Kallikratikoi/", sep="")	
setwd(muncDir)		
mapMunc <- st_read("Kallikratikoi.shp")
mapMunc <- st_transform(mapMunc, crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
mapMunc %>% select(geometry) %>%
				plot()
setwd(dir)

# Convert names from Greek to Latin script
#latinNames <- stri_trans_general(mapMunc$LEKTIKO, "Greek-Latin")
#write.csv(latinNames, file = "muncNames.csv", row.names = FALSE)	

# Format data for mapping
thisData <- data0 %>% filter(is.na(homeLat)==FALSE) # exclude cases with missing latitude/longitude
data1 <- thisData %>% group_by(week=ceiling_date(onsetDate, "week")) # group by week for plotting (this is 'week ending...')
hospData <- data0 %>% filter(is.na(hosp1Lat)==FALSE) # exclude cases with missing latitude/longitude
mapHome <- st_as_sf(data1, coords = c('homeLon', 'homeLat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
mapHosp <- st_as_sf(hospData, coords = c('hosp1Lon', 'hosp1Lat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
#mapSchool <- st_as_sf(data1, coords = c('schoolLon', 'schoolLat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")

# Compare with LatLong matches reported region (safety check as Lat and Long derived from place name)
#matrix1 <- sf::st_intersects(mapHome, mapGreece, sparse = FALSE) # matrix showing which NUTS-3 region point belongs to.
matrix1 <-  sf::st_intersects(mapHome, mapMunc, sparse = FALSE) # matrix showing which municipality point belongs to.

# N.B. Does it matter that these values are not planar? Could convert (see end of script)
#colnames(matrix1) <- mapGreece$geo # assign names of NUTS-3 regions
colnames(matrix1) <- mapMunc$KALCODE # assign names of NUTS-3 regions
rownames(matrix1) <- mapHome$ourID
thisData <- thisData %>% mutate(findMunc = colnames(matrix1)[apply(matrix1,1,which.max)]) # Find name of NUTS-3 region corresponding to Lat. and Long.

checkThese <- thisData %>% filter(findRegion != placeOfResidence) %>%
			select(ourID, placeOfResidence, homeTown, homeLat, homeLon, findRegion, placeOfNotification) #%>%
			#print(n=1000)
#write.csv(checkThese, file = "180509checkPlaces.csv", row.names = FALSE)			   

# Plot MMR uptake with cases
map2 <- tm_shape(mapGreece) + tm_fill("lightgrey") +
tm_shape(mapCoverage, is.master = TRUE) +
tm_polygons("classCoverage1", title = "MMR uptake",
palette = "Greens", border.col = "white") +
#tm_text("geo", just = "center") +
#tm_scale_bar()# +
tm_format_Europe(legend.position = c("right", "top"), attr.outside = TRUE,
				 main.title = "Uptake of at least 1 dose of MMR vaccine (2015 birth cohort)", main.title.size = 1)+
 tm_shape(mapHome) + tm_dots(col="darkred", size=0.05)#+
# tm_shape(mapHosp) + tm_dots(col="darkblue", size=0.075) 
 
 map2
	
plotCases <- qtm(mapGreece) + tm_shape(mapHome) + tm_dots(col="darkred", size=0.05) +
							  tm_shape(mapHosp) + tm_dots(col="darkblue", size=0.075) 
print(plotCases)

plot_anim <- qtm(mapGreece) + tm_shape(mapHome) + tm_dots(col="darkred", size=0.05) + tm_facets(by = "week", free.coords = FALSE, nrow = 6, ncol = 10) # plot cases by week
# Require ImageMagick for animation
#tmap_animation(plot_anim, filename = "plot_anim.gif", delay = 25)


###########################
## Calculate shortest    ##
## distance between each ##
## pair of cases         ##
###########################

# List of locations for transmitting (include hospitals where patient)
locVars <- c("home", "school", "hospLoc1", "hospLoc2", "hospLoc3", "expLoc1", "expLoc2", "HCWLoc")
distanceData <- thisData %>% select(ID, homeTown, homeLat, homeLon, hosp1, hosp1Lat, hosp1Lon, hosp2, hosp2Lat, hosp2Lon, hosp3, hosp3Lat, hosp3Lon, schoolCode, schoolLat, schoolLon, exp1Code, exp1Lat, exp1Lon, exp2Code, exp2Lat, exp2Lon, HCWHospCode, HCWHospLat, HCWHospLon, onsetDate, hospitalDate, hospitalisationLag) %>%
							 unite(home, c(homeTown, homeLat, homeLon), sep = "*", remove = TRUE) %>%
							 unite(school, c(schoolCode, schoolLat, schoolLon), sep = "*", remove = TRUE) %>%
							 # add notification place?
							 unite(hospLoc1, c(hosp1, hosp1Lat, hosp1Lon), sep = "*", remove = TRUE) %>%
							 unite(hospLoc2, c(hosp2, hosp2Lat, hosp2Lon), sep = "*", remove = TRUE) %>%
							 unite(hospLoc3, c(hosp3, hosp3Lat, hosp3Lon), sep = "*", remove = TRUE) %>%
							 unite(expLoc1, c(exp1Code, exp1Lat, exp1Lon), sep = "*", remove = TRUE) %>%
							 unite(expLoc2, c(exp2Code, exp2Lat, exp2Lon), sep = "*", remove = TRUE) %>%
							 unite(HCWLoc, c(HCWHospCode, HCWHospLat, HCWHospLon), sep = "*", remove = TRUE) %>%
							 gather(key = location, value = locLat, locVars) %>%   # gather so that there is one entry for each reference to a hospital
						     arrange(ID) %>%
							 separate(locLat, c("place", "lat", "lon"), sep="[*]") %>%
							 filter(lat != "NA") 

distanceData <- distanceData %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon)) # no. of decimal places
splitDistance <- split(distanceData, as.numeric(distanceData$ID))
caseIndices <- distanceData %>% select(ID) %>% distinct(ID) 							 

# List of locations for being infected
infVars <- c("home", "school", "expLoc1", "expLoc2", "HCWLoc")
infLocData <- thisData %>% select(ID, homeTown, homeLat, homeLon, schoolCode, schoolLat, schoolLon, exp1Code, exp1Lat, exp1Lon, exp2Code, exp2Lat, exp2Lon, HCWHospCode, HCWHospLat, HCWHospLon, onsetDate, hospitalDate, hospitalisationLag) %>%
							 unite(home, c(homeTown, homeLat, homeLon), sep = "*", remove = TRUE) %>%
							 unite(school, c(schoolCode, schoolLat, schoolLon), sep = "*", remove = TRUE) %>%
							 # add notification place?
							 unite(expLoc1, c(exp1Code, exp1Lat, exp1Lon), sep = "*", remove = TRUE) %>%
							 unite(expLoc2, c(exp2Code, exp2Lat, exp2Lon), sep = "*", remove = TRUE) %>%
							 unite(HCWLoc, c(HCWHospCode, HCWHospLat, HCWHospLon), sep = "*", remove = TRUE) %>%
							 gather(key = location, value = locLat, infVars) %>%   # gather so that there is one entry for each reference to a hospital
						     arrange(ID) %>%
							 separate(locLat, c("place", "lat", "lon"), sep="[*]") %>%
							 filter(lat != "NA") 

infLocData <- infLocData %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon)) # no. of decimal places
splitInfDistance <- split(infLocData, as.numeric(infLocData$ID))

## Find shortest distance between all locations of transmitter and pre-admission locations of receiver
distIndexCase <- function(indexCase){
  # Co-ordinates for thisCase
  thisCase <- splitDistance[[indexCase]] %>% select(lon,lat)
  
  # Reduce to cases for which distance has not yet been calculated
  otherCases <- caseIndices %>% filter(caseIndices > indexCase)
  otherCases <- as.character(otherCases$ID)  
   
  distanceLocation <- function(oneLocation){
    oneLocationM <- as.matrix(oneLocation)     
    
    # Compute distance between this location for thisCase and each location of otherCase
      distanceLocationCase <- function(oneLocCase2){
      coordsCase2 <- as.matrix(splitInfDistance[[oneLocCase2]] %>% select(lon,lat)) 
	  distLocCase2 <- distVincentyEllipsoid(coordsCase2, oneLocationM)  # calculate in metres using Vincenty distance
      return(min(distLocCase2)) # shortest distance between unPoint and country which index is unIndicePays2
   }
   
   # Find minimum distance to each other case
   distLocCase2 <- lapply(otherCases, distanceLocationCase)
   minDistances <- unlist(distLocCase2)
   return(minDistances)
  }
  
  # Create table with shortest distance to each other case (use 'by_row' from purrrlyr)
  distanceCases <- by_row(thisCase, distanceLocation, .collate="cols")
  distanceCases <- select(distanceCases, -c(lat,lon)) # remove lat and lon. of index case
  
  
  # Shortest distances between oneLocation and every other case
  if(count(distanceCases)==1){
    # For the last case on the list
    shortestDistance <- distanceCases
  }else{
    shortestDistance <- as_tibble(apply(distanceCases, 2, min))
  }
  
  result <- bind_cols(case1 = rep(indexCase, nrow(shortestDistance)),case2 = otherCases, dist = shortestDistance)
  return(result)
}


# Matrix of distances between each pair of cases
distanceMatrix <- lapply(caseIndices$ID[-length(caseIndices$ID)], distIndexCase)
distanceMatrix <- bind_rows(distanceMatrix)
distanceMatrix <- distanceMatrix %>% spread(case2,value)

####################################################

#############################
## Likelihood calculation: ##
## Wallinga-Teunis plus    ##
## gravity model           ##
#############################


# Likelihood based on date of symptom onset alone
dates1 <- thisData$onsetDate; names(dates1) <- thisData$ourID

lags <- outer(dates1,dates1, FUN = "-")
lags <- as.numeric(lags, units="days")
dim(lags) <- c(length(dates1),length(dates1))
like1 <- dnorm(lags, mean1, sd1)
like1[col(like1)==row(like1)] <- 0  # set diagonal elements to zero as not included in denom. of normalisation

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))


minRow <- apply(like1, 1, FUN=min)
nLike1 <- like1/pmax(minRow,rowSums(like1)) # normalise and adjust for potential zero sum here
nLike1[is.nan(nLike1)] <- 0
nLike1[1,] <- 0  # set likelihood of first case being infected by future cases to zero 
thisData$RnTime <- colSums(nLike1) # calculation of effective reproduction number


# Likelihood based on distance alone (symmetrical, no arrow of time)
rho <- 2 # Sensitivity analysis. Alter to 0.5, 1 and 5
like2 <- 1/(distanceMatrix)^rho # power law
like2[is.infinite(like2)] <- 1  # replace infinite values with 1 (pole)
# N.B. see Meyer & Held for alternatives to the basic power law which causes such a pole at x=0

minRow2 <- apply(like2, 1, FUN=min)
nLike2 <- like2/pmax(minRow2,rowSums(like2)) # normalise and adjust for potential zero sum here
nLike2[is.nan(nLike2)] <- 0
thisData$RnPlace <- colSums(nLike2) # calculation of effective reproduction number based on distance

# Combined likelihood based on distance and time
like3 <- nLike1*nLike2
minRow3 <- apply(like3, 1, FUN=min)
nLike3 <- like3/pmax(minRow3,rowSums(like3)) # normalise and adjust for potential zero sum here
nLike3[is.nan(nLike3)] <- 0

thisData$RnCombined <- colSums(nLike3) # calculation of effective reproduction number

# Secondary infections before / after hospital admission. **Do not use this!**
preAd <- outer(dates1, thisData$hospitalDate, FUN="-")  # Time lag between date of hospitalisation and symptom onset in secondary case
preAd <- as.numeric(preAd, units="days")
dim(preAd) <- c(length(dates1),length(thisData$hospitalDate))
preAd <- as_tibble(preAd)
postAd <- preAd
preAd[preAd <= -11] <- 1 # select for secondary infections before hospital admission (assume 14 days before rash appears. Term this symptom onset)
preAd[preAd > -11] <- 0 
preAd <- preAd * nLike3

postAd[postAd > -11] <- 1 # select for secondary infections following hospital admission 
postAd[postAd <= -11] <- 0 
postAd <- postAd * nLike3

thisData <- thisData %>% mutate(RnPreAd = colSums(preAd)) # effective reproduction number for cases caused prior to hospital admission
thisData <- thisData %>% mutate(RnPostAd = colSums(postAd)) # effective reproduction number for cases caused post hospital admission

thisData %>% filter(is.na(hospitalDate) =="FALSE") %>% 
			select(RnCombined, RnPreAd, RnPostAd)

# Plot effective reproduction number over time
ggplot( data = thisData, aes(onsetDate, RnCombined)) +
				geom_point(col=ECDCcol[1]) +
				#geom_smooth(span=0.15) +
				geom_hline(yintercept = 1, lty = 2, lwd=1.2) +
				labs(
					x = "Date of onset of symptoms",
					y = "Effective reproduction number of case j (Rj)") +
				theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=18)) 

# Plot histogram of Rn 
ggplot(data = thisData, aes(x=RnCombined)) +
  geom_histogram(binwidth = 0.1, fill=ECDCcol[1]) +
  labs(	x = "Effective reproduction number", y = "Count") +
  theme(legend.position="none",  text = element_text(size=18))

median(thisData$RnCombined)
var(thisData$RnCombined)	


##################
## Matrix plots ##
##################	
	
## Group data by age
thisData$age1 <- ifelse(is.na(thisData$monthAge), thisData$age, thisData$monthAge/12) # Updated age variable to account for age in months for 0-1 year olds
thisData$age1[thisData$age1==0] <- 0.00001 # workaround for babies being only classified as 0 

lower <- 0
upper <- 60
sep = "-"
above.char = "+"
by <- 5
 
labs <- c("0-1", "2-4", paste(seq(5, upper - by, by = by),
                 seq(5 + by - 1, upper - 1, by = by),
                 sep = sep), "60-69", "70-79",
           paste(80, above.char, sep = ""))
		   
thisData$ageGroup <- cut(thisData$age1, breaks = c(0,2,seq(5,upper,by=5),70,80,Inf), labels = labs, include.lowest=FALSE, right=FALSE)

	
## Code age groups
thisData <- thisData %>%      
  mutate(popGroup = case_when(isRoma==1    ~ 1,    	## Roma
							  isMigrant==1 ~ 2,		## Migrant
							  isHCW==1     ~ 3,		## Healthcare worker
							  TRUE         ~ 4))	## General population

# Create matrix of age groups
popMatrix <- data.matrix(nLike3)
popMatrix <- cbind(thisData$popGroup, nLike3) # Define matrix with population group in first row and column and likelihood of i being infected by j in body
popMatrix <- rbind(c(0,thisData$popGroup), popMatrix) 


# Sum the number of inferred cases in each row for each population group
nPopGroups <- length(unique(thisData$popGroup))   # No. of age groups. Check: will this work if one group is missing?
a <- rep(-1, (nPopGroups*nrow(nLike3)))
dim(a) <- c(nrow(nLike3),nPopGroups)

for (i in 1:nPopGroups){
a[,i] <- rowSums(popMatrix[-1,which(popMatrix[,1]==i)]) # Group population groups of 'causes'
}


WAIFWPop <- rep(-1, nPopGroups^2)   ## Normalise
dim(WAIFWPop) <- c(nPopGroups,nPopGroups)

b <- cbind(thisData$popGroup,a)

for(j in 1:nPopGroups){
d <- aggregate(b[,(j+1)], by=list(Category=b[,1]), FUN=sum)
WAIFWPop[,j] <- d[c(1:nPopGroups),2]
}

WAIFWPop[is.na(WAIFWPop)] <- 0

groupPopCases <- rowSums(WAIFWPop)  # total no. of cases in each age group
groupPopCauses <- colSums(WAIFWPop)  # total no. of cases caused by each age group

## Plot WAIFW matrix by age
melted_WAIFWPop <- melt(WAIFWPop)
melted_WAIFWPop <- as.data.frame(melted_WAIFWPop)

# Rename population groups
popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")
melted_WAIFWPop$Var1 <- factor(melted_WAIFWPop$Var1,
levels = c(1:length(popLabs)),
labels = popLabs)
melted_WAIFWPop$Var2 <- factor(melted_WAIFWPop$Var2,
levels = c(1:length(popLabs)),
labels = popLabs)

head(melted_WAIFWPop)

#create new variable from counts
WAIFWPopFactor <- cut(melted_WAIFWPop$value, 
				  breaks = c(0, 1, 10, 50, 100, 500, 1000, 10000),
				  labels = c("0", "1-9", "10-49", "50-99", "100-499", "500-999", "1000+"), include.lowest=TRUE)
				   
melted_WAIFWPop$value <- factor(WAIFWPopFactor)	
colours1 <- colorRampPalette(ECDCcol)(length(unique(melted_WAIFWPop$value)))
N <- nlevels(melted_WAIFWPop$value)
	
ggplot(data = melted_WAIFWPop, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "No. of cases", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Number of cases attributable to each population group", 
    x = "Population group of infected cases",
    y = "Population group of infecting cases") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(angle=60, hjust=1.1),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
  
  
	
## Plot normalised WAIFW matrix	
WAIFWPop1 <- sweep(WAIFWPop, 1, groupPopCases, `/`) # normalise by number of cases in each age group 
melted_WAIFWPop1 <- melt(WAIFWPop1)
melted_WAIFWPop1 <- as.data.frame(melted_WAIFWPop1)

# Rename age groups
melted_WAIFWPop1$Var1 <- factor(melted_WAIFWPop1$Var1,
levels = c(1:length(popLabs)),
labels = popLabs)
melted_WAIFWPop1$Var2 <- factor(melted_WAIFWPop1$Var2,
levels = c(1:length(popLabs)),
labels = popLabs)

head(melted_WAIFWPop1)

#create new variable from counts
lower <- 0
upper <- 1
sep = "-"
above.char = "+"
byThis <- 0.05
 
labels2 <- c(paste(seq(lower, upper - byThis, by = byThis),
                 seq(lower + byThis, upper, by = byThis),
                 sep = sep))
		
WAIFWPop1Factor <- cut(melted_WAIFWPop1$value, 
				   breaks = seq(lower,upper, by=byThis),
				   labels = labels2,
				   include.lowest=TRUE, dig.lab = 3)
				   
melted_WAIFWPop1$value <- factor(WAIFWPop1Factor)	
colours1 <- colorRampPalette(ECDCcol)(length(unique(melted_WAIFWPop1$value)))
N <- nlevels(melted_WAIFWPop1$value)
	
ggplot(data = melted_WAIFWPop1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "Prop. of cases", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Proportion of cases attributable to each population group", 
    x = "Population group of infected cases",
    y = "Population group of infecting cases") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
  
# Group hospitalisation lag as proxy for severity
thisData <- thisData %>% 
  mutate(lagGroup = case_when(hospitalisationLag < 0  ~ 1,		## admitted before symptom onset
							  hospitalisationLag == 0 ~ 2,		## admitted at symptom onset
							  hospitalisationLag == 1 | hospitalisationLag == 2 ~ 3,	## 1-2 days after symptom onset
							  hospitalisationLag >= 3 & hospitalisationLag <= 7 ~ 4,	## 3-7 days after symptom onset
							  hospitalisationLag >= 8 | hospitalisationLag == 14 ~ 5,	## 8-14 days after symptom onset
							  hospitalisationLag >= 15 ~ 6))	## more than two weeks after symptom onset

thisData %>% filter(hospitalisationLag != "NA") %>% select(hospitalisationLag, lagGroup)
  
# Effective reproduction numbers by group
thisData <- thisData %>% mutate(isHospitalised = factor(isHospitalised, levels=c("Y","N"))) # assign factors
thisData <- thisData %>% mutate(isRoma = factor(isRoma, levels=c(0,1))) # assign factors
thisData <- thisData %>% mutate(isHCW = factor(isHCW, levels=c(0,1))) # assign factors
thisData <- thisData %>% mutate(isMigrant = factor(isMigrant, levels=c(0,1))) # assign factors
thisData <- thisData %>% mutate(lagGroup = factor(lagGroup, levels=c(1,2,3,4,5,6))) # assign factors

# *Repeat for each group*
thisData %>% filter(isHospitalised == "N") %>%
			 group_by(popGroup) %>%
			 #group_by(lagGroup) %>%
		     summarise(meanRn = mean(RnCombined), medianRn = median(RnCombined), sdRn = sd(RnCombined))

# Density plots of distribution of Rn		
# ** Repeat for each risk factor **  
ggplot(thisData %>% filter(isRoma == 1), aes(x=RnCombined, fill=(isHospitalised))) +
  geom_density(alpha=0.4) + 
  labs(	x = "Effective reproduction number",
		y = "Probability density") +
  theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=16)) +
  scale_x_continuous(limits=c(0, 8), breaks=seq(0,upper,by=1)) 


# Plot effective reproduction number by age group and population group
popAge <- thisData %>% filter(isHospitalised != "NA") %>% group_by(popGroup, ageGroup) %>%
		  summarise(meanRn = mean(RnCombined), medianRn = median(RnCombined), sdRn = sd(RnCombined))

popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")
popAge$popGroup <- factor(popAge$popGroup,
					levels = c(1:length(popLabs)),
					labels = popLabs)
	

lower <- 0
upper <- 3
sep = "-"
above.char = "+"
byThis <- 0.25
 
popAgeLabs <- c(paste(seq(lower, upper - byThis, by = byThis),
                      seq(lower + byThis, upper, by = byThis),
                      sep = sep))
		   
popAge$meanRn <- cut(popAge$meanRn, breaks = seq(lower,upper,by=byThis), labels = popAgeLabs, include.lowest=TRUE, right=FALSE)
				   

colours1 <- colorRampPalette(ECDCcol)(length(unique(popAge$meanRn)))
N <- nlevels(popAge$meanRn)
	
ggplot(data = popAge, aes(x=ageGroup, y=popGroup, fill=meanRn)) + 
  geom_tile() +
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "Rn", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Effective reproduction number", x = "Age group", y = "Population group") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(angle=60, hjust=1.1),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") +
   scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
  
  


# Anderson-Darling tests to compare distributions
roma <- thisData %>% filter(isRoma == 1) %>% select(RnCombined)
nonRoma <-  thisData %>% filter(isRoma == 0) %>% select(RnCombined)
ad.test(roma$RnCombined, nonRoma$RnCombined)

HCW <- thisData %>% filter(isHCW == 1) %>% select(RnCombined)
nonHCW <-  thisData %>% filter(isHCW == 0) %>% select(RnCombined)
ad.test(HCW$RnCombined, nonHCW$RnCombined)

migrant <- thisData %>% filter(isMigrant == 1) %>% select(RnCombined)
nonMigrant <-  thisData %>% filter(isMigrant == 0) %>% select(RnCombined)
ad.test(migrant$RnCombined, nonMigrant$RnCombined)

hospitalised <- thisData %>% filter(isHospitalised == "Y") %>% select(RnCombined)
nonHospitalised <-  thisData %>% filter(isHospitalised == "N") %>% select(RnCombined)
ad.test(hospitalised$RnCombined, nonHospitalised$RnCombined)

lag1 <- thisData %>% filter(lagGroup == 1) %>% select(RnCombined)
lag2 <- thisData %>% filter(lagGroup == 2) %>% select(RnCombined)
lag3 <- thisData %>% filter(lagGroup == 3) %>% select(RnCombined)
lag4 <- thisData %>% filter(lagGroup == 4) %>% select(RnCombined)
lag5 <- thisData %>% filter(lagGroup == 5) %>% select(RnCombined)
lag6 <- thisData %>% filter(lagGroup == 6) %>% select(RnCombined)
ad.test(lag1$RnCombined, lag2$RnCombined, lag3$RnCombined, lag4$RnCombined, lag5$RnCombined)

# Plot histogram of age of cases
ggplot(data = thisData, aes(age1)) +
  geom_histogram(binwidth = 1, fill=ECDCcol[1]) +
  labs(	x = "Age of case",
		y = "Count") +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=18))  +
  scale_x_continuous(limits=c(0, 80), breaks=seq(0,80,by=5)) 
  
# Create matrix of age groups
ageMatrix <- data.matrix(nLike3)
ageMatrix <- cbind(as.numeric(thisData$ageGroup), nLike3) # Define matrix with age group in first row and column and likelihood of i being infected by j in body
ageMatrix <- rbind(c(0,as.numeric(thisData$ageGroup)), ageMatrix) 


# Sum the number of inferred cases in each row for each age group
nGroups <- length(unique(thisData$ageGroup))   # No. of age groups. Check: will this work if one group is missing?
a <- rep(-1, (nGroups*nrow(nLike3)))
dim(a) <- c(nrow(nLike3),nGroups)

for (i in 1:nGroups){
a[,i] <- rowSums(ageMatrix[-1,which(ageMatrix[,1]==i)]) # Group ages of 'causes'
}


WAIFW <- rep(-1, nGroups^2)   ## Normalise
dim(WAIFW) <- c(nGroups,nGroups)

b <- cbind(as.numeric(thisData$ageGroup),a)

for(j in 1:nGroups){
d <- aggregate(b[,(j+1)], by=list(Category=b[,1]), FUN=sum)
WAIFW[,j] <- d[c(1:nGroups),2]
}

WAIFW[is.na(WAIFW)] <- 0

groupCases <- rowSums(WAIFW)  # total no. of cases in each age group
groupCauses <- colSums(WAIFW)  # total no. of cases caused by each age group
nAgeGp <- count(as.numeric(thisData$ageGroup)) # alternative count of number of cases in each age group (check)
groupRn <- groupCauses/nAgeGp$freq[c(1:nGroups)] # average effective reproduction number by age group


library(reshape2)
library(stringr)

## Plot WAIFW matrix by age
melted_WAIFW <- melt(WAIFW)
melted_WAIFW <- as.data.frame(melted_WAIFW)

# Rename age groups
melted_WAIFW$Var1 <- factor(melted_WAIFW$Var1,
levels = c(1:length(labs)),
labels = labs)
melted_WAIFW$Var2 <- factor(melted_WAIFW$Var2,
levels = c(1:length(labs)),
labels = labs)

head(melted_WAIFW)

#create new variable from counts
WAIFWFactor <- cut(melted_WAIFW$value, 
				   breaks = c(0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130),
				   labels = c("0", "1-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100-109", "110-119", "120+"), include.lowest=TRUE)
				   
melted_WAIFW$value <- factor(WAIFWFactor)	
colours1 <- colorRampPalette(ECDCcol)(length(unique(melted_WAIFW$value)))
N <- nlevels(melted_WAIFW$value)
	
ggplot(data = melted_WAIFW, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "No. of cases", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Number of cases attributable to each age group", 
    x = "Age group of infected cases",
    y = "Age group of infecting cases") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(angle=60, hjust=1.1),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
  
  
	
## Plot normalised WAIFW matrix	
WAIFW1 <- sweep(WAIFW, 1, groupCases, `/`) # normalise by number of cases in each age group 
melted_WAIFW1 <- melt(WAIFW1)
melted_WAIFW1 <- as.data.frame(melted_WAIFW1)

# Rename age groups
melted_WAIFW1$Var1 <- factor(melted_WAIFW1$Var1,
levels = c(1:length(labs)),
labels = labs)
melted_WAIFW1$Var2 <- factor(melted_WAIFW1$Var2,
levels = c(1:length(labs)),
labels = labs)

head(melted_WAIFW1)

#create new variable from counts
lower <- 0
upper <- 1
sep = "-"
above.char = "+"
byThis <- 0.05
 
labels2 <- c(paste(seq(lower, upper - byThis, by = byThis),
                 seq(lower + byThis, upper, by = byThis),
                 sep = sep))
		
WAIFW1Factor <- cut(melted_WAIFW1$value, 
				   breaks = seq(lower,upper, by=byThis),
				   labels = labels2,
				   include.lowest=TRUE, dig.lab = 3)
				   
melted_WAIFW1$value <- factor(WAIFW1Factor)	
colours1 <- colorRampPalette(ECDCcol)(length(unique(melted_WAIFW1$value)))
N <- nlevels(melted_WAIFW1$value)
	
ggplot(data = melted_WAIFW1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "Prop. of cases", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Proportion of cases attributable to each age group", 
    x = "Age group of infected cases",
    y = "Age group of infecting cases") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(angle=60, hjust=1.1),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") 
  




i <- 1419
plot(nLike3[i,])
points(nLike1[i,], col="red")
points(nLike2[i,], col="blue")





#####################################
## Code for finding lat. and long. ##
#####################################

distinctPlace <- data0 %>% distinct(homeTown) %>% top_n()

# Obtain lat. and long. from Google. It is possible to obtain further information on boundaries of the location etc.
#distinctPlace <- distinctPlace %>% mutate(geocode(distinctPlace$homeTown))
a <- as_tibble(geocode(distinctPlace$homeTown, output = "all") %>% 
			   select(lat, lon, postal_code)) %>%
			   mutate(postal_code = as.character(postal_code))
distinctPlace2 <- bind_cols(distinctPlace,a)

# Manual check for missing values (iterate since often due to Google API query limit)
distinctPlace1 <- distinctPlace2 %>% filter(is.na(distinctPlace2$lat))
b <- as_tibble(geocode(distinctPlace1$homeTown, output = "more") %>% 
			   select(lat, lon, postal_code)) 
distinctPlace1 <- bind_cols(distinctPlace1 %>% select(homeTown),b)

distinctPlace2 <- left_join(distinctPlace2, distinctPlace1, by="homeTown") %>% 
				   transmute(homeTown, lat = ifelse(is.na(lat.y), lat.x, lat.y),
									   lon = ifelse(is.na(lon.y), lon.x, lon.y),
									   postal_code = ifelse(is.na(postal_code.y), postal_code.x, postal_code.y)) 

write.csv(distinctPlace2, file = "postcodes.csv", row.names = FALSE)			   

# Postcodes from lat. and lon.
postcodes <- revgeo(latLong$lon, latLong$lat, output="frame") %>% select(zip)

# Date checking
preOnset <- thisData %>% filter(hospitalDate < onsetDate)
write.csv(preOnset, file = "preOnsetHospitalisation.csv", row.names = FALSE)	

thisData %>% filter(notifyDate < onsetDate)
thisData %>% filter(dischargeDate < hospitalDate)
thisData %>% filter(deathDate < onsetDate)



planarData <- st_transform(mapData, 2163)
planarGreece <- st_transform(mapGreece, 2163)

matrix2 <- st_intersects(planarData, planarGreece, sparse = FALSE)


localMap <- st_read("local.shp")


#######################
## Hospital analysis ##
#######################
hospVars <- c("hosp1", "hosp2", "hosp3", "exp1Code", "exp2Code", "HCWHospCode")
hospData <- thisData %>% gather(key = anyHospital, value = code, hospVars) %>%   # gather so that there is one entry for each reference to a hospital
						drop_na(code) %>% 
						arrange(ourID)
						  
hospData %>% group_by(code) %>%
            
hospitals <- hospData %>% group_by(code) %>%
						  summarise(nCases = n())
						  
  
## Histogram of cases per hospital
ggplot(data = hospitals, aes(x=nCases)) +
  geom_bar() +
  labs(	x = "", y = "Count") +
  
  theme(legend.position="none")# +
  
  
 