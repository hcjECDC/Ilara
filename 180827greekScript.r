#rm(list=ls())
dir <- "P:/Measles/Greek data/"
setwd(dir)
library(plyr)
library(ggmap) #for LatLong API (Google Maps)
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
options(pillar.sigfig = 6) # set number of significant figures to be displayed as standard
api <- "AIzaSyBn2K2--zdOLXfhkBqSXr6XAl7T0t4Pc-M" # Key for using Google Maps API to obtain geo. data
register_google(key = api)

## In matrices, i (row) is infected by j (column)

## Create custom colour scale
myColours <- c("26 107 133", "241 214 118", "168 45 23")
ECDCcol <- sapply(strsplit(myColours, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))  # convert to hexadecimal

## Measles parameters
mean1 <- 11.7 # mean of normal distribution for serial interval (days) ##Vink et al. 2014
sd1 <- 2 # standard deviation of normal distribution for serial interval (days)  ##Vink et al. 2014

## Read in line list
data0 <- read_csv("workingGreekData.csv", col_names=TRUE)
data0$onsetDate <- as.Date(data0$onsetDate, "%d/%m/%Y")  # convert to date format
data0$notifyDate <- as.Date(data0$notifyDate, "%d/%m/%Y")  # convert to date format
data0$hospitalDate <- as.Date(data0$hospitalDate, "%d/%m/%Y")  # convert to date format
data0$hosp1Lat <- as.numeric(data0$hosp1Lat)
data0$hosp1Lon <- as.numeric(data0$hosp1Lon)
data0$exp1Lat <- as.numeric(data0$exp1Lat)
data0$exp1Lon <- as.numeric(data0$exp1Lon)
data0$exp2Lat <- as.numeric(data0$exp2Lat)
data0$exp2Lon <- as.numeric(data0$exp2Lon)
data0$HCWHospLat <- as.numeric(data0$HCWHospLat)
data0$HCWHospLon <- as.numeric(data0$HCWHospLon)
data0$deathDate <- as.Date(data0$deathDate, "%d/%m/%Y")  # convert to date format
data0$dischargeDate <- as.Date(data0$dischargeDate, "%d/%m/%Y")  # convert to date format

data0 <- data0[order(data0$onsetDate),] # order by date of symptom onset 

write.csv(data0, file = "180509greekData.csv", row.names = FALSE)			

MMRData <- read_csv("MMR2.csv", col_names=TRUE, col_types = cols(resPostcode = col_character()))  # need to specify residencePostcode as character here because two entries include letters, not just numbers. Otherwise parser would assume it was a column of integers.
popData <- read_csv("pop2015.csv", col_names=TRUE) 
muncPopData <- read_csv("popMunc.csv", col_names=TRUE) 
muncPopData$muncCode[c(1:41)] <- sprintf("%04d", muncPopData$muncCode[c(1:41)]) # re-format to match mapMunc
#latLong <- read_csv("LatLong.csv", col_names=TRUE) # With '180509greekData.csv' co-ordinates and postcodes are included

##################
## MMR coverage ##
##################

## Rearrange prescribing data from long to wide format (group all vaccine codes together (all MMR))
data2 <- MMRData %>% 
 mutate(purchaseDate = substr(purchaseDate, 1, 10))%>% # drop times
 select(-initialOrder, -vaccineCode) %>% # Need to drop this otherwise record number coerces R into keeping all rows
 arrange(id, purchaseDate) %>% # Ensure that dates are in correct order before inferring dose number
 group_by(id) %>% 
 dplyr::mutate(dose = row_number()) %>% # Dose number for spreading data
 dplyr::mutate(doses = max(row_number())) %>% # Total number of doses for record
 unite("datePlace", c(purchaseDate, pharmPostcode, pharmDistrict)) %>%
 spread(dose, datePlace, fill="")
	
## Separate out date and place fields. Label columns
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
						  
# Number of doses by postcode
oneDosePC <- data3 %>% filter(resPostcode != "#N/A" & doses==1) %>%
					arrange(resPostcode) %>%
					group_by(resPostcode) %>%
					summarise(oneDose = n()) %>%
					rename(postcode = resPostcode) 
					
twoPlusDosesPC <- data3 %>% filter(resPostcode != "#N/A" & doses>1) %>%
						  arrange(resPostcode) %>%
						  group_by(resPostcode) %>%
					      summarise(twoPlusDoses = n())%>%
						  rename(postcode = resPostcode) 
						  

# Incidence by postcode
postcodeRatio <- thisData %>% filter(age <= 3) %>%
							    group_by(postcode) %>% 
								mutate(nCases=n()) %>% 
								distinct(postcode, .keep_all=TRUE) %>% # this way we keep all other variables, as opposed to summarise where we lose them
								select(postcode, nCases, muncCode)

postcodeRatio <- left_join(postcodeRatio, oneDosePC, by="postcode")
postcodeRatio <- left_join(postcodeRatio, twoPlusDosesPC, by="postcode") %>%
							mutate(ratio1 = nCases/(oneDose+twoPlusDoses),
								   ratio2 = nCases/twoPlusDoses)

propPopGroup <- thisData %>% filter(age <= 3) %>%
							  group_by(muncCode, popGroup) %>% 
							  dplyr::summarise(n = n()) %>%
							  mutate(freq = n / sum(n)) %>%
							  filter(popGroup == 1) # select proportion of Roma
							  
 

muncRatio <- postcodeRatio %>% group_by(muncCode) %>%
								summarise(nCases = sum(nCases),
										  oneDose = sum(oneDose),
										  twoPlusDoses = sum(twoPlusDoses)) %>%
								mutate(ratio1 = nCases/(oneDose+twoPlusDoses),
									   ratio2 = nCases/twoPlusDoses,
									   logRatio1 = log(ratio1),
									   logRatio2 = log(ratio2))

muncRatio <- left_join(muncRatio, propPopGroup)	
							   
muncRatio$muncCode[c(1:41)] <- sprintf("%04d", muncRatio$muncCode[c(1:41)]) # re-format to match mapMunc						
				
						  
# Population by postcode
popPostcode <- popData %>% filter(is.na(postcode)==FALSE) %>%
							arrange(postcode) %>%
							group_by(postcode) %>%
							summarise(total = sum(popCohort15)) %>%
							mutate(relSize = total/max(total))
write.csv(popPostcode, file = "popPostcode.csv", row.names = FALSE)								
					   
pop2015 <- popData %>% arrange(NUTS3Code) %>%
					   group_by(NUTS3Code) %>%
					   summarise(Total = sum(popCohort15))

# MMR coverage at NUTS-3 level					   
MMRcoverage <- full_join(pop2015, oneDose, by="NUTS3Code")
MMRcoverage <- full_join(MMRcoverage, twoPlusDoses, by="NUTS3Code") %>%
			   rename(population = Total.x) %>%
			   rename(oneDose = Total.y) %>%
			   rename(twoPlusDoses = Total) %>%
			   mutate(coverage1 = oneDose/population) %>%
			   mutate(coverage2plus = twoPlusDoses/population) %>%
			   mutate(overallCoverage = coverage1 + coverage2plus) %>%
			   mutate(RnVaccine = 15 * (1-overallCoverage)) %>%
			   mutate(unvaccinated = population * (1-overallCoverage))


## 2010 NUTS-3 code (Attiki) was subdivided into 7 codes in 2013 but we need to update population data. Average for now.
athensCodes <- c(paste("EL30", seq(1,7), sep=""))	
#athensCodes1 <- athensCodes[c(1:6)] # omit Piraeus
		   
expandAthens <- MMRcoverage %>% select(-NUTS3Code) %>%
								slice(rep(1, each=7)) %>%	          # copy Attiki entry 7 times
								mutate(NUTS3Code = athensCodes) %>%	  #	assign list of 7 Attiki codes 
								select(names(MMRcoverage))
								
										   
MMRcoverage <- 	bind_rows(MMRcoverage[-1,], expandAthens)	%>%	   
			    rename(geo = NUTS3Code)   # to join with Eurostat shapefile 
				
	
infantCasesPop <- thisData %>% filter(age == 2) %>%
							  group_by(popGroup, placeOfResidence) %>% 
							  #group_by(postcode) %>% 
							  summarise(nCases = n()) %>%
							  mutate(geo = placeOfResidence) %>%
							  select(geo, popGroup, nCases)
							  #select(postcode, nCases)
							  
MMRcoverage <- full_join(MMRcoverage, infantCasesPop, by="geo") %>%
					mutate(proportion = 100 * (nCases / unvaccinated)) 

MMRcoverage$popGroup <- factor(MMRcoverage$popGroup,
							levels = c(1:length(popLabs)),
							labels = popLabs)

					
# Plot proportion of susceptibles who are infected. Caution! Could be mismatch with population denominator due to migration or age cohort effects (2015 birth vs. 2 years old)
ggplot(MMRcoverage, aes(x=geo, y=proportion)) +
     geom_point() +
	 labs(x = "NUTS-3 region",
	   y = "nCases",
	   fill = "") +
	 theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14)) 

					

# MMR coverage at postcode resolution
MMRpopCoverage <- full_join(popPostcode, oneDosePC, by="postcode")
MMRpopCoverage <- full_join(MMRpopCoverage, twoPlusDosesPC, by="postcode") %>%
					rename(population = popPostcode) %>%
					mutate(coverage1 = oneDose/population) %>%
					mutate(coverage2plus = twoPlusDoses/population) %>%  
					rowwise() %>% 
					mutate(overallCoverage = sum(coverage1,coverage2plus, na.rm=TRUE)) %>%
					mutate(overallCoverage = min(overallCoverage,1)) %>% # discrepancy with population denominator means a couple of postcodes have coverage > 1. Recode.
					mutate(unvaccinated = population * (1-overallCoverage))
						
									   
## Define classes for plotting coverage		
## 1 dose
lower <- 0.6
upper <- 1
sep = "-"
by <- 0.05
 
labs1 <- c("0-0.59", paste(seq(lower, upper - 0.05, by = by),
                 seq(lower + by - 0.01, upper, by = by),
                 sep = sep))
		   
MMRcoverage$classCoverage1 <- cut(MMRcoverage$overallCoverage, breaks = c(0, seq(lower,upper,by=0.1)), labels = labs2, include.lowest=FALSE, right=FALSE)
   
## 2+ doses
lower <- 0.3
upper <- 1
sep = "-"
by <- 0.1
 
labs2 <- c("0-0.29", paste(seq(lower, upper - 0.1, by = by),
                 seq(lower + by - 0.01, upper, by = by),
                 sep = sep))
		   
MMRcoverage$classCoverage2 <- cut(MMRcoverage$coverage2plus, breaks = c(0, seq(lower,upper,by=0.1)), labels = labs2, include.lowest=FALSE, right=FALSE)


## RnVaccine
lower <- 0
upper <- 7
sep = "-"
by1 <- 1
 
labs1 <- c(paste(seq(lower, (upper - by1), by = by1),
                 seq((lower + by1), upper, by = by1),
                 sep = sep))
		   
MMRcoverage$classRnVaccine <- cut(MMRcoverage$RnVaccine, breaks = c(seq(lower,upper,by=by1)), labels = labs1, include.lowest=FALSE, right=FALSE)
							   
## Define classes for municipality population	
breaks3 <- c(100, 1000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 1000000)
labs3 <- c("100-999", "1000-9999", "10000-19999", "20000-29999", "30000-39999", "40000-49999", "50000-59999", "60000-69999", "70000-79999", "80000-89999", "90000-99999", "100000-1000000")
		   
muncPopData$classPop <- cut(muncPopData$muncPop, breaks = breaks3, labels = labs3, include.lowest=FALSE, right=FALSE)
   
								   
## Define classes for log (case / vaccine ratio)	
lower <- -8
upper <- 0
sep = " to "
by1 <- 1
 
labs1 <- c(paste(seq(lower, (upper - by1), by = by1),
                 seq((lower + by1), upper, by = by1),
                 sep = sep))
		   
muncRatio$classLogRatio1 <- cut(muncRatio$logRatio1, breaks = c(seq(lower,upper,by=by1)), labels = labs1, include.lowest=FALSE, right=FALSE)
muncRatio$classLogRatio2 <- cut(muncRatio$logRatio2, breaks = c(seq(lower,upper,by=by1)), labels = labs1, include.lowest=FALSE, right=FALSE)
      
###################
## Plotting maps ##			   
###################

## Plot Greece
mapGreece <- get_eurostat_geospatial(output_class = "sf", resolution = "60", nuts_level = 3) # Bug in EUROSTAT package (16/11) that prevents 
mapGreece <- mapGreece %>% filter(CNTR_CODE == "EL")
mapGreece %>% select(geo) %>%
				plot()

## Plot Athens
mapAthens <- mapGreece %>% filter(NUTS_ID %in% athensCodes) 
mapAthens %>% select(geo) %>%
				plot()
				
## Plot MMR coverage				  
mapCoverage <- left_join(mapGreece, MMRcoverage, by="geo") %>%   # combine geographic data with coverage
			   filter(geo != "NA")

mapCoverageAthens <- left_join(mapAthens, MMRcoverage, by="geo") %>%   # combine geographic data with coverage
			   filter(geo != "NA")

## Reading in shapefiles (different spatial resolutions)
## Local communities
#localDir <- paste(dir,"Shapefiles/LocalCommunities/", sep="")	
#setwd(localDir)		
#localMap <- st_read("local.shp")
#localMap %>% select(geometry) %>%
#				plot()
#setwd(dir)

## Kallikratikoi
muncDir <- paste(dir,"Shapefiles/Kallikratikoi/", sep="")	
setwd(muncDir)		
mapMunc <- st_read("Kallikratikoi.shp")
mapMunc <- st_transform(mapMunc, crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
mapMunc %>% select(geometry) %>%
				plot()
setwd(dir)


# Format data for mapping
thisData <- data0 %>% filter(is.na(homeLat)==FALSE) %>% # exclude cases with missing latitude/longitude
                      arrange(ourID)

## Code and factorise population groups
thisData <- thisData %>%      
  mutate(popGroup = case_when(isRoma==1    ~ 1,    	## Roma
							  isMigrant==1 ~ 2,		## Migrant
							  isHCW==1     ~ 3,		## Healthcare worker
							  TRUE         ~ 4))	## General population
popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")

	
data1 <- thisData %>% group_by(postcode, popGroup, week=ceiling_date(onsetDate, "week")) %>% # group by week for plotting (this is 'week ending...')
						mutate(nCases=n()) %>% 
						distinct(postcode, .keep_all=TRUE) %>% # this way we keep all other variables, as opposed to summarise where we lose them
						select(postcode, popGroup, nCases, week, homeLat, homeLon, placeOfResidence, age, muncCode)
						
data1$popGroup <- factor(data1$popGroup,
					levels = c(1:length(popLabs)),
					labels = popLabs)
					
hospData <- data0 %>% filter(is.na(hosp1Lat)==FALSE) # exclude cases with missing latitude/longitude
mapHome <- st_as_sf(data1, coords = c('homeLon', 'homeLat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
mapHosp <- st_as_sf(hospData, coords = c('hosp1Lon', 'hosp1Lat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
#mapSchool <- st_as_sf(data1, coords = c('schoolLon', 'schoolLat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")

# N.B. Does it matter that these values are not planar? Could convert (see end of script)
# Find region of data point from lat. and long.
# Compare with LatLong matches reported region (safety check as Lat and Long derived from place name)
matrix1 <- sf::st_intersects(mapHome, mapGreece, sparse = FALSE) # matrix showing which NUTS-3 region point belongs to.
colnames(matrix1) <- mapGreece$geo # assign names of NUTS-3 regions
rownames(matrix1) <- mapHome$ourID
thisData <- thisData %>% mutate(findRegion = colnames(matrix1)[apply(matrix1,1,which.max)]) # Find name of NUTS-3 region corresponding to Lat. and Long.

map1 <- st_as_sf(distanceData, coords = c('lon', 'lat'), crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
st_crs(map1) <- 4326
matrix1 <- sf::st_intersects(map1, mapGreece, sparse = FALSE) # matrix showing which NUTS-3 region point belongs to.
colnames(matrix1) <- mapGreece$geo # assign names of NUTS-3 regions
distanceData <- distanceData %>% mutate(geo = colnames(matrix1)[apply(matrix1,1,which.max)]) # Find name of NUTS-3 region corresponding to Lat. and Long.


## Check whether inferred lat and long match with region
#checkThese <- thisData %>% filter(findRegion != placeOfResidence) %>%
			#select(ourID, placeOfResidence, homeTown, homeLat, homeLon, findRegion, placeOfNotification)# %>%
			#print(n=1000)
#write.csv(checkThese, file = "180509checkPlaces.csv", row.names = FALSE)			   

#######################
## Population values ##
#######################

# Find municipality of data point from lat. and long.
matrix2 <-  sf::st_intersects(mapHome, mapMunc, sparse = FALSE) # matrix showing which municipality point belongs to.
colnames(matrix2) <- muncPopData$muncName # assign names of municipalities
rownames(matrix2) <- mapHome$ourID
thisData <- thisData %>% mutate(muncName = colnames(matrix2)[apply(matrix2,1,which.max)]) # Find name of municipality corresponding to Lat. and Long.

# Find population of municipality
thisData <- left_join(thisData, muncPopData, by="muncName")

# Find population of postcode
#N.B. This is 2015 birth cohort data and is only being used as a test of principle for coding!!
popPostcode <- popPostcode %>% mutate(postcode = as.character(postcode)) %>%
							   rename(popPostcode = total)
thisData <- left_join(thisData, popPostcode, by="postcode") 
thisData <-  thisData %>% replace_na(list(popPostcode = 1)) # We do not have 2015 birth cohort population data for all postcodes. Assume the missing ones are very small (1). 

################
## Plot cases ##
################

## Greece cases
map2 <- qtm(mapGreece) +  		
			tm_shape(mapHome) +
				tm_symbols(col="popGroup", 
						   size="nCases",
						   palette=ECDCcol,
						   shape=20,
						   title.col = c("Population group"),
						   title.size = c("Number of cases")) + 
			tm_layout(legend.position = c("right", "centre"), 
						legend.title.size = 1,
						legend.text.size = 0.75,
						attr.outside = TRUE,
						main.title = "Distribution of measles cases by population group", 
						main.title.size = 1) 

print(map2) 
	
## Weekly static
plot_anim <- qtm(mapGreece) + 		
				tm_shape(mapHome) +
					tm_dots(col="popGroup", size="nCases") +
				tm_facets(by = "week", free.coords = FALSE, nrow=6, ncol = 10)  # plot cases by week
				#tm_facets(along = "week")  # plot cases by week. Used for animation

## Weekly animation				
# Require ImageMagick for animation
#tmap_animation(plot_anim, filename = "plot_anim.gif", delay = 25)


## Athens cases
map2A <- qtm(mapAthens) +  		
			tm_shape(mapHome %>% filter(placeOfResidence %in% athensCodes)) +
				tm_symbols(col="popGroup",
						   size="nCases",
						   palette=ECDCcol,
						   shape=20,
						   title.col = c("Population group"),
						   title.size = c("Number of cases")) + 
			tm_layout(legend.show = FALSE,	
					  title = "Athens",
					  title.size = 1,
					  title.position = c("right", "bottom"))

print(map2A) 
	
#######################
## Plot MMR coverage ##
#######################
# MCV1 in Greece
map3 <- qtm(mapGreece) +  		
			tm_shape(mapCoverage, is.master = TRUE) +
				tm_polygons(col = "classCoverage1", title = "MCV1 coverage", palette = "Greens", border.col = "white", alpha=0.8) +
			tm_shape(mapHome %>% filter(age <= 3)) + 
				tm_dots(col="black", shape=20, size="nCases") +
			tm_layout(legend.position = c("right", "centre"), 
						legend.title.size = 1.2,
						legend.text.size = 0.75,
						attr.outside = TRUE,
						main.title = "MCV1 coverage and incident cases in children aged 3 and under", 
						main.title.size = 1.2) 

print(map3) 
	
# MCV1 in Athens
map3A <- qtm(mapAthens) +
			tm_shape(mapCoverageAthens, is.master = TRUE) +
				tm_polygons(col = "classCoverage1", title = "MCV1 coverage", palette = "Greens", border.col = "white") +
			tm_shape(mapHome %>% filter(placeOfResidence %in% athensCodes & age <= 3)) +
				tm_dots(col="black", shape=20, size="nCases") +
			tm_layout(legend.show = FALSE,	
					  title = "Athens",
					  title.size = 1,
					  title.position = c("right", "bottom"))	
	
print(map3A)	

# MCV2 in Greece
map4 <- qtm(mapGreece) +  		
			tm_shape(mapCoverage, is.master = TRUE) +
				tm_polygons(col = "classCoverage2", title = "MCV2 coverage", palette = "Greens", border.col = "white", alpha=0.8) +
			tm_shape(mapHome %>% filter(age <= 3)) + 
				tm_dots(col="black", shape=20, size="nCases") +
			tm_layout(legend.position = c("right", "centre"), 
						legend.title.size = 1.2,
						legend.text.size = 0.75,
						attr.outside = TRUE,
						main.title = "MCV2 coverage and incident cases in children aged 3 and under", 
						main.title.size = 1.2) 

print(map4) 
	
# MCV2 in Athens
map4A <- qtm(mapAthens) +
			tm_shape(mapCoverageAthens, is.master = TRUE) +
				tm_polygons(col = "classCoverage2", title = "MCV2 coverage", palette = "Greens", border.col = "white") +
			tm_shape(mapHome %>% filter(placeOfResidence %in% athensCodes & age <= 3)) +
				tm_dots(col="black", shape=20, size="nCases") +
			tm_layout(legend.show = FALSE,	
					  title = "Athens",
					  title.size = 1,
					  title.position = c("right", "bottom"))	
	
print(map4A)	

# RnVaccine in Greece
map5 <- qtm(mapGreece) +  		
			tm_shape(mapCoverage, is.master = TRUE) +
				tm_polygons(col = "classRnVaccine", title = "Rn", palette = "Reds", border.col = "white", alpha=0.8) +
			tm_layout(legend.position = c("right", "centre"), 
						legend.title.size = 1.2,
						legend.text.size = 0.75,
						attr.outside = TRUE,
						main.title = "Effective reproduction number according to 2015 MCV1 coverage", 
						main.title.size = 1.2) 

print(map5) 
	
# RnVaccine in Athens
map5A <- qtm(mapAthens) +
			tm_shape(mapCoverageAthens, is.master = TRUE) +
				tm_polygons(col = "classRnVaccine", title = "Rn from MCV1 coverage", palette = "Reds", border.col = "white") +
				tm_layout(legend.show = FALSE,	
					  title = "Athens",
					  title.size = 1,
					  title.position = c("right", "bottom"))	
	
print(map5A)	

# Plot ratio of cases to vaccine doses by municipality
muncRatio <- muncRatio %>% mutate(KALCODE = as.factor(muncCode)) %>%
							 select(-muncCode)
							 
mapMuncRatio <- left_join(mapMunc, muncRatio, by="KALCODE") %>%   # combine geographic data with population
					filter(oneDose != "NA")

map5 <- qtm(mapMunc) + 
			tm_shape(mapMuncRatio, is.master = TRUE) +
				tm_polygons(col = "classLogRatio1", title = "Log ratio of cases to number vaccinated", palette = "Greens", border.col = "white") # not strictly no.of doses
print(map5) 
	

###########################
## Calculate shortest    ##
## distance between each ##
## pair of cases         ##
###########################

# List of locations for transmitting (include hospitals where patient)
locVars <- c("home", "school", "hospLoc1", "hospLoc2", "hospLoc3", "expLoc1", "expLoc2", "HCWLoc")
distanceData <- thisData %>% select(ourID, homeTown, homeLat, homeLon, hosp1, hosp1Lat, hosp1Lon, hosp2, hosp2Lat, hosp2Lon, hosp3, hosp3Lat, hosp3Lon, schoolCode, schoolLat, schoolLon, exp1Code, exp1Lat, exp1Lon, exp2Code, exp2Lat, exp2Lon, HCWHospCode, HCWHospLat, HCWHospLon, onsetDate, hospitalDate, hospitalisationLag) %>%
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
						     arrange(ourID) %>%
							 separate(locLat, c("place", "lat", "lon"), sep="[*]") %>%
							 filter(lat != "NA") 

distanceData <- distanceData %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon)) # total list of cases by location
splitDistance <- split(distanceData, as.numeric(distanceData$ourID)) # list of locations for each case (including co-ordinates, place, dates and type of location)
caseIndices <- distanceData %>% select(ourID) %>% distinct(ourID) 							 

# List of locations for being infected
infVars <- c("home", "school", "expLoc1", "expLoc2", "HCWLoc")
infLocData <- thisData %>% select(ourID, homeTown, homeLat, homeLon, schoolCode, schoolLat, schoolLon, exp1Code, exp1Lat, exp1Lon, exp2Code, exp2Lat, exp2Lon, HCWHospCode, HCWHospLat, HCWHospLon, onsetDate, hospitalDate, hospitalisationLag) %>%
							 unite(home, c(homeTown, homeLat, homeLon), sep = "*", remove = TRUE) %>%
							 unite(school, c(schoolCode, schoolLat, schoolLon), sep = "*", remove = TRUE) %>%
							 # add notification place?
							 unite(expLoc1, c(exp1Code, exp1Lat, exp1Lon), sep = "*", remove = TRUE) %>%
							 unite(expLoc2, c(exp2Code, exp2Lat, exp2Lon), sep = "*", remove = TRUE) %>%
							 unite(HCWLoc, c(HCWHospCode, HCWHospLat, HCWHospLon), sep = "*", remove = TRUE) %>%
							 gather(key = location, value = locLat, infVars) %>%   # gather so that there is one entry for each reference to a hospital
						     arrange(ourID) %>%
							 separate(locLat, c("place", "lat", "lon"), sep="[*]") %>%
							 filter(lat != "NA") 

infLocData <- infLocData %>% mutate(lat = as.numeric(lat), lon = as.numeric(lon)) # no. of decimal places
splitInfDistance <- split(infLocData, as.numeric(infLocData$ourID)) # list of locations where each case may have been infected (including co-ordinates, place, dates and type of location)

## Find shortest distance between all locations of infecter and pre-infection locations of infectee
distIndexCase <- function(indexCase){

  # Co-ordinates for thisCase
  thisCase <- splitDistance[[toString(indexCase)]] %>% select(location,lon,lat)
  
  # Reduce to cases for which distance has not yet been calculated
  otherCases <- as_tibble(caseIndices %>% filter(caseIndices > indexCase)) # N.B. when looking at a smaller number of cases for debugging, this yeilds a moving window
  otherCases <- as.character(otherCases$ourID)  
   
  # Compute distance between this location for thisCase and each location of otherCase
  distanceBetweenCases <- function(otherCase){
	  coordsOtherCase <- as.matrix(splitInfDistance[[otherCase]] %>% select(lon,lat)) 
      mat <- distm(thisCase[,c('lon','lat')], coordsOtherCase[,c('lon','lat')], fun=distVincentyEllipsoid) / 1000 # distances in km using Vincenty
	  colnames(mat) <- as.list(splitInfDistance[[otherCase]]$location) # label columns with locations of infectee
	  rownames(mat) <- as.list(thisCase$location)
	  distances <- melt(mat) # distances between each pair of locations
	  #cat(c("Calculating distance between",indexCase,"and",otherCase, "\n"))
	  return(suppressMessages({distances %>% top_n(-1) %>% filter(1:n() == 1)})) # return the row of the minimum distance (only the first in the case of a tie)
   }
   
   # Find minimum distance to each other case
   distBetween <- lapply(otherCases, distanceBetweenCases)
   minDistances <- as_tibble(suppressMessages({melt(distBetween)})) %>% select(Var1, Var2, value) # shortest distance between each pair of cases, including location (e.g. home, hospital, school)
   minDistances <- minDistances %>% mutate(indexCase = indexCase, otherCase = otherCases)
   minDistances <- minDistances %>%
					unite("Combine", c(indexCase, otherCase, Var1, Var2, value))
					
  # minDistances <- bind_cols(case1 = rep(indexCase, nrow(minDistances)),case2 = otherCases, dist = minDistances)
   return(minDistances)
   
}


# Matrix of distances between each pair of cases
distanceMatrix <- lapply(caseIndices$ourID[-length(caseIndices$ourID)],distIndexCase)

distanceMatrix <- bind_rows(distanceMatrix)%>% 
					separate(Combine, c("thisCase", "otherCase", "loc1","loc2","distance"), sep="_") %>% ## where 'loc1' is the location of the transmitter and 'loc2' is the location of the receiver
					gather(variable, value, -(thisCase:otherCase)) %>%
					unite(temp, otherCase, variable) %>%
					spread(temp, value)%>%
					#arrange(thisCase) %>% # workaround to avoid implicit sorting in 'spread'
					mutate(thisCase = as.factor(thisCase))  %>%
					mutate_at(vars(matches("distance")), as.double) # since gathering/uniting/spreading has coerced to character. SLOW!!

## Start here!!


distNames <- paste(sort(thisData$ourID[-1]), "_distance", sep="") # N.B. These are now in order of ourID
dist1 <- distanceMatrix %>% select(thisCase,distNames) # distances only, re-order columns
add_row(dist1, thisCase = as.character(max(as.numeric(thisData$ourID)))) # row for highest ID case (all distances calculated by other cases)

dist2 <- dist1 %>% mutate(ourID = as.numeric(levels(dist1$thisCase))) %>%  # re-order by ourID
				   arrange(ourID)	%>%
				   add_row(ourID = max(as.numeric(thisData$ourID))) %>%
				   mutate('1_distance' = as.numeric(rep(NA, nrow(thisData)))) %>%
				   select('1_distance', everything()) %>%
				   select(- c(ourID, thisCase))

# N.B. Do we want to re-format these as the distance matrix?
loc1Names <- paste(sort(thisData$ourID[-1]), "_loc1", sep="")
loc1 <- distanceMatrix %>% select(thisCase, loc1Names) # location of transmitter
 

loc2Names <- paste(sort(thisData$ourID[-1]), "_loc2", sep="")
loc2 <- distanceMatrix %>% select(thisCase, loc2Names) # location of receiver



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
thisData <- thisData %>% mutate(RnTime = colSums(nLike1)) # calculation of effective reproduction number


# Likelihood based on distance alone (symmetrical, no arrow of time)
beta1 <- 0.5 #scaling factor for distance decay 
rho <- 2 # Sensitivity analysis. Alter to 0.5, 1 and 5
scalePop1 <- 0.5
scalePop2 <- 0.5
popValue  <- thisData$muncPop # Choose which population metric to use here e.g. popPostcode, muncPop
popMatrix <- (scalePop1 * popValue) %o% (scalePop2 * popValue) # outer product of population of home municipality

# Choose distance function here
#like2 <- 1/(dist2)^rho # power law 
#like2 <- popMatrix / 1#(dist2)^rho  # gravity model
like2 <- popMatrix * exp(-beta1 * dist2)  # gravity model

like2[is.infinite(like2)] <- 1  # replace infinite values with 1 (pole)
like2[lower.tri(like2)] = t(like2)[lower.tri(like2)] # make symmetrical 
diag(like2) <- 1

minRow2 <- apply(like2, 1, FUN=min)
nLike2 <- like2/pmax(minRow2,rowSums(like2)) # normalise and adjust for potential zero sum here
nLike2[is.na(nLike2)] <- 0
thisData <- thisData %>% mutate(RnPlace = colSums(nLike2)) # calculation of effective reproduction number based on distance

# Combined likelihood based on distance and time
like3 <- like1*like2
minRow3 <- apply(like3, 1, FUN=min)
nLike3 <- like3/pmax(minRow3,rowSums(like3)) # normalise and adjust for potential zero sum here
nLike3[is.nan(nLike3)] <- 0

thisData <- thisData %>% mutate(RnCombined = colSums(nLike3)) # calculation of effective reproduction number

## Tidy formatting of data
# Likelihood matrix
likeMatrix <- gather(as_tibble(nLike3)) ## !! Double-check this transposition !!

# Location matrix
locMatrix <- as_tibble(loc1) %>%
			 mutate(ourID = as.integer(levels(thisCase)),
					"1_loc1" = rep("NA", nrow(loc1))) %>%
			 select(-thisCase) %>%
			 select(ourID, "1_loc1", everything()) %>%
			 arrange(ourID) %>%
			 rbind(c(2315, rep(NA, ncol(loc1)))) 

locMatrix1 <- as_tibble(t(locMatrix %>% select(-ourID)))
			 
tidyData <- gather(locMatrix1) %>%
				mutate(toID = rep(thisData$ourID, times=nrow(thisData)),
					   fromID = rep(thisData$ourID, each=nrow(thisData)),
					   likelihood = likeMatrix$value) %>%
				rename(location = value) %>%
				select(toID, fromID, location, likelihood) %>%

infected <- thisData %>% select(ourID, age1, popGroup, placeOfNotification) %>%
						 rename(toID = ourID,
							    toAge = age1,
								toPopGroup = popGroup,
								toPlace = placeOfNotification)

infected <- right_join(infected, tidyData) %>%
			select(toID, toAge, toPopGroup, toPlace)

infector <- tidyData %>% select(fromID, location) %>%
						 rename(ourID = fromID) 
						 
infector <- left_join(infector, distanceData) %>% 
				select(ourID, location, geo) 

infector <- left_join(infector, thisData) %>% 
				select(ourID, age1, popGroup, location, geo) %>%
				rename(fromID = ourID,
					   fromAge = age1,
					   fromPopGroup = popGroup,
					   fromPlace = geo)

tidyData1 <- bind_cols(infected, infector, tidyData %>% select(likelihood))




# Secondary infections before / after hospital admission. **Do not use this!**
#preAd <- outer(dates1, thisData$hospitalDate, FUN="-")  # Time lag between date of hospitalisation and symptom onset in secondary case
#preAd <- as.numeric(preAd, units="days")
#dim(preAd) <- c(length(dates1),length(thisData$hospitalDate))
#preAd <- as_tibble(preAd)
#postAd <- preAd
#preAd[preAd <= -11] <- 1 # select for secondary infections before hospital admission (assume 14 days before rash appears. Term this symptom onset)
#preAd[preAd > -11] <- 0 
#preAd <- preAd * nLike3

#postAd[postAd > -11] <- 1 # select for secondary infections following hospital admission 
#postAd[postAd <= -11] <- 0 
#postAd <- postAd * nLike3

#thisData <- thisData %>% mutate(RnPreAd = colSums(preAd)) # effective reproduction number for cases caused prior to hospital admission
#thisData <- thisData %>% mutate(RnPostAd = colSums(postAd)) # effective reproduction number for cases caused post hospital admission

#thisData %>% filter(is.na(hospitalDate) =="FALSE") %>% 
#			select(RnCombined, RnPreAd, RnPostAd)

thisData$popGroup <- as.factor(thisData$popGroup)
# Plot effective reproduction number over time
ggplot(data = thisData, aes(x = onsetDate, y = RnCombined)) +
				geom_point(aes(fill = popGroup), pch=21, alpha = 0.8) +
				#geom_smooth(span=0.15) +
				geom_hline(yintercept = 1, lty = 2, lwd=1.2) +
				labs(
					x = "Date of onset of symptoms",
					y = "Effective reproduction number of case j (Rj)") +
				theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=14)) 
				  
# Plot histogram of Rn 
ggplot(data = thisData, aes(x=RnCombined)) +
  geom_histogram(aes(fill=popGroup), alpha = 0.8, bins=50) +
  labs(	x = "Effective reproduction number", y = "Count") +
  theme(legend.position="right",  text = element_text(size=14))

median(thisData$RnCombined)
var(thisData$RnCombined)	

# Bubble plot to show population group and number of attributable infections
# Figure S1
ggplot(data = thisData, aes(x = onsetDate, y = popGroup)) + 
  geom_point(aes(fill = popGroup, size = RnCombined), pch=21, alpha = 0.4) +
  scale_size_continuous(range = c(0.5, 10)) +  # Adjust the range of points size 
  #scale_y_discrete("Population group", labels = c(1 = "Roma", 2 = "Migrant", 3 = "Healthcare worker", 4 = "General population"))
  labs(x = "Date of onset of symptoms",
	   y = "Population group",
	   size = "Effective reproduction number") +
  theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14)) 
	
# Plot population of municipality by population group
ggplot(data = thisData, aes(x=muncPop, fill=popGroup)) +
geom_density(alpha=0.4) + 
  labs(	x = "Population of municipality", y = "Density") +
  theme(legend.position="right",  text = element_text(size=18))
  

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

## Code population groups (This is now done before map plotting, still included here in case bugs crop up later...)
#thisData <- thisData %>%      
#  mutate(popGroup = case_when(isRoma==1    ~ 1,    ## Roma
#							  isMigrant==1 ~ 2,		## Migrant
#							  isHCW==1     ~ 3,		## Healthcare worker
#							  TRUE         ~ 4))	## General population

# Create matrix of population groups
popMatrix <- data.matrix(nLike3)
popMatrix <- cbind(thisData$popGroup, nLike3) # Define matrix with population group in first row and column and likelihood of i being infected by j in body
popMatrix <- rbind(c(0,thisData$popGroup), popMatrix) 

# Sum the number of inferred cases in each row for each population group
nPopGroups <- length(unique(thisData$popGroup))   # No. of population groups. 
a <- rep(-1, (nPopGroups*nrow(nLike3)))
dim(a) <- c(nrow(nLike3),nPopGroups)

for (i in 1:nPopGroups){
a[,i] <- rowSums(popMatrix[-1,which(popMatrix[,1]==i)]) # Group population groups of 'causes' i.e. this is the sum of all cases in each group caused by the person on this row.
}


WAIFWPop <- rep(-1, nPopGroups^2)  
dim(WAIFWPop) <- c(nPopGroups,nPopGroups)

b <- cbind(thisData$popGroup,a) # re-bind with the population group of the case on each row

for(j in 1:nPopGroups){
d <- aggregate(b[,(j+1)], by=list(Category=b[,1]), FUN=sum) # sum all the cases which are caused by each population group in turn
WAIFWPop[,j] <- d[c(1:nPopGroups),2] # Cases are caused by column.
}

WAIFWPop[is.na(WAIFWPop)] <- 0

groupPopCases <- rowSums(WAIFWPop)  # total no. of cases in each population group
groupPopCauses <- colSums(WAIFWPop)  # total no. of cases caused by each population group

## Plot WAIFW matrix by population group
melted_WAIFWPop <- melt(WAIFWPop)
melted_WAIFWPop <- as.data.frame(melted_WAIFWPop)

# Rename population groups
popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")
melted_WAIFWPop$Var1 <- factor(melted_WAIFWPop$Var1,   #infectee
levels = c(1:length(popLabs)),
labels = popLabs)
melted_WAIFWPop$Var2 <- factor(melted_WAIFWPop$Var2,	# infecter
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
  

## Plot Rn matrix	
WAIFWPop1 <- sweep(WAIFWPop, MARGIN=1, STATS = groupPopCases, FUN =`/`) # number of secondary cases caused by each group
melted_WAIFWPop1 <- melt(WAIFWPop1)
melted_WAIFWPop1 <- as.data.frame(melted_WAIFWPop1)

# Rename population groups
melted_WAIFWPop1$Var1 <- factor(melted_WAIFWPop1$Var1,
levels = c(1:length(popLabs)),
labels = popLabs)
melted_WAIFWPop1$Var2 <- factor(melted_WAIFWPop1$Var2,
levels = c(1:length(popLabs)),
labels = popLabs)

head(melted_WAIFWPop1)

#create new variable from counts
lower <- 0
upper <- 2
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
  scale_fill_manual(values=colours1, guide = guide_legend(keywidth = 3, title = "Rn", label.theme = element_text(face="plain", angle=0))) + 
   labs(title = "Effective reproduction number attributable to each population group", 
    x = "Population group of infected cases",
    y = "Population group of infecting cases") +
  theme(text = element_text(size=16, face="bold"),
        axis.text.x = element_text(),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))
   
# Rn by NUTS-3 and population group. Use place of notification!!
# Create matrix of NUTS-3 regions by population group
popMatrix <- data.matrix(nLike3)
popMatrix <- cbind(thisData$popGroup, thisData$popGroup, nLike3) # Define matrix with population group in first row and column and likelihood of i being infected by j in body
popMatrix <- rbind(c(0,thisData$popGroup), popMatrix) 

# Sum the number of inferred cases in each row for each population group
nPopGroups <- length(unique(thisData$popGroup))   # No. of population groups. 
a <- rep(-1, (nPopGroups*nrow(nLike3)))
dim(a) <- c(nrow(nLike3),nPopGroups)

for (i in 1:nPopGroups){
a[,i] <- rowSums(popMatrix[-1,which(popMatrix[,1]==i)]) # Group population groups of 'causes' i.e. this is the sum of all cases in each group caused by the person on this row.
}


WAIFWPop <- rep(-1, nPopGroups^2)  
dim(WAIFWPop) <- c(nPopGroups,nPopGroups)

b <- cbind(thisData$popGroup,a) # re-bind with the population group of the case on each row

for(j in 1:nPopGroups){
d <- aggregate(b[,(j+1)], by=list(Category=b[,1]), FUN=sum) # sum all the cases which are caused by each population group in turn
WAIFWPop[,j] <- d[c(1:nPopGroups),2] # Cases are caused by column.
}

WAIFWPop[is.na(WAIFWPop)] <- 0

groupPopCases <- rowSums(WAIFWPop)  # total no. of cases in each population group
groupPopCauses <- colSums(WAIFWPop)  # total no. of cases caused by each population group

## Plot WAIFW matrix by population group
melted_WAIFWPop <- melt(WAIFWPop)
melted_WAIFWPop <- as.data.frame(melted_WAIFWPop)

# Rename population groups
popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")
melted_WAIFWPop$Var1 <- factor(melted_WAIFWPop$Var1,   #infectee
levels = c(1:length(popLabs)),
labels = popLabs)
melted_WAIFWPop$Var2 <- factor(melted_WAIFWPop$Var2,	# infecter
levels = c(1:length(popLabs)),
labels = popLabs)

head(melted_WAIFWPop)





nNUTS3 <- length(unique(thisData$placeOfResidence)) # no. of municipalities
WAIFWNUTS3 <- rep(-1, (nPopGroups*nNUTS3))  
dim(WAIFWNUTS3) <- c(nNUTS3,nPopGroups)

b <- cbind(as.factor(thisData$placeOfResidence),a) # re-bind with the NUTS-3 region of the case on each row

for(j in 1:nPopGroups){
d <- aggregate(b[,(j+1)], by=list(Category=b[,1]), FUN=sum) # sum all the cases which are caused by each population group in turn
WAIFWNUTS3[,j] <- d[c(1:nNUTS3),2] # Cases are caused by column.
}

WAIFWNUTS3[is.na(WAIFWNUTS3)] <- 0

melted_NUTS3 <- melt(WAIFWNUTS3)
melted_NUTS3 <- as.data.frame(melted_NUTS3)

# Rename NUTS-3 regions
melted_NUTS3$geo <- factor(melted_NUTS3$Var1,   #infectee
						levels = c(1:length(unique(thisData$placeOfResidence))),
						labels = unique(thisData$placeOfResidence)) 
					 
# Rename population groups
melted_NUTS3$popGroup <- factor(melted_NUTS3$Var2,	# infecter
						levels = c(1:length(popLabs)),
						labels = popLabs)
						
melted_NUTS3 <- melted_NUTS3 %>%
				rename(nCasesAttributed = value) %>%
				select(geo, popGroup, nCasesAttributed)
				
head(melted_NUTS3)

# Number of cases in each population group in each NUTS-3 region
NUTS3Gp <- thisData %>% mutate(geo = placeOfResidence) %>%
						select(-placeOfResidence) %>%
						group_by(geo, popGroup) %>%
					    summarise(nCases = n()) 
						
NUTS3Gp$popGroup <- factor(NUTS3Gp$popGroup,	# infecter
						levels = c(1:length(popLabs)),
						labels = popLabs)

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
thisData %>% filter(isRoma == "Y") %>%
			 group_by(popGroup) %>% # This doesn't make sense when filtering by isRoma
			 #group_by(lagGroup) %>%
		     summarise(meanRn = mean(RnCombined), medianRn = median(RnCombined), sdRn = sd(RnCombined))

# Density plots of distribution of Rn		
# ** Repeat for each risk factor **  
ggplot(thisData %>% filter(isMigrant == 1), aes(x=RnCombined)) + #, fill=(isHospitalised))) +
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
#popAge <- thisData %>% filter(isHospitalised != "NA" & popGroup!=2) %>% group_by(popGroup, ageGroup) %>% # omitting migrants to allow rescaled plot
		  summarise(meanRnTime = mean(RnTime), medianRnTime = median(RnTime), sdRnTime = sd(RnTime),
					meanRnPlace = mean(RnPlace), medianRnPlace = median(RnPlace), sdRnPlace = sd(RnPlace),
					meanRnCombined = mean(RnCombined), medianRnCombined = median(RnCombined), sdRnTime = sd(RnCombined))

popLabs <- c("Roma", "Migrant", "Healthcare worker", "General population")
popAge$popGroup <- factor(popAge$popGroup,
					levels = c(1:length(popLabs)),
					labels = popLabs)
	

lower <- 0
upper <- 7
sep = "-"
above.char = "+"
byThis <- 0.25
 
popAgeLabs <- c(paste(seq(lower, upper - byThis, by = byThis),
                      seq(lower + byThis, upper, by = byThis),
                      sep = sep))
		   
popAge$RnCombGroup <- cut(popAge$meanRnCombined, breaks = seq(lower,upper,by=byThis), labels = popAgeLabs, include.lowest=TRUE, right=FALSE)
				   

colours1 <- colorRampPalette(ECDCcol)(length(unique(popAge$RnCombGroup)))
N <- nlevels(popAge$RnCombGroup)
	
ggplot(data = popAge, aes(x=ageGroup, y=popGroup, fill=RnCombGroup)) + 
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
  theme(text = element_text(size=14, face="bold"),
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
  theme(text = element_text(size=14, face="bold"),
        axis.text.x = element_text(angle=60, hjust=1.1),
		axis.title.x = element_text(vjust = 1.5),
		legend.position = "right") 
 
# Effective reproduction number

# Comparison of likelihood (time, space, combined
i <- 1419
plot(nLike3[i,], pch=20, col=ECDCcol[1])
points(nLike1[i,], pch=20, col=ECDCcol[2])
points(nLike2[i,], pch=20, col=ECDCcol[3])


#####################################
## Code for finding lat. and long. ##
#####################################

# List all locations 
distinctPlace <- data0 %>% distinct(homeTown, .keep_all=TRUE) %>% rename(place = homeTown, lat = homeLat, lon = homeLon) %>% select(place, lat, lon)
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(hosp1Name, .keep_all=TRUE) %>% rename(place = hosp1Name, lat = hosp1Lat, lon = hosp1Lon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(hosp2Name, .keep_all=TRUE) %>% rename(place = hosp2Name, lat = hosp2Lat, lon = hosp2Lon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(hosp3Name, .keep_all=TRUE) %>% rename(place = hosp3Name, lat = hosp3Lat, lon = hosp3Lon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(HCWHospName, .keep_all=TRUE) %>% rename(place = HCWHospName, lat = HCWHospLat, lon = HCWHospLon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(exp1Name, .keep_all=TRUE) %>% rename(place = exp1Name, lat = exp1Lat, lon = exp1Lon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(exp2Name, .keep_all=TRUE) %>% rename(place = exp2Name, lat = exp2Lat, lon = exp2Lon) %>% select(place, lat, lon)) 
distinctPlace <- bind_rows(distinctPlace, data0 %>% distinct(schoolCode, .keep_all=TRUE) %>% rename(place = schoolCode, lat = schoolLat, lon = schoolLon) %>% select(place, lat, lon)) 

# Obtain lat. and long. from Google. It is possible to obtain further information on boundaries of the location etc.
#distinctPlace <- distinctPlace %>% mutate(geocode(distinctPlace$homeTown))
a <- as_tibble(geocode(distinctPlace, output = "all", source="google")) %>%   # can choose between source="google" and source = "dsk"
			   mutate(postal_code = as.character(postal_code))
distinctPlace2 <- bind_cols(distinctPlace,a)


#### OLD ####
# Manual check for missing values (iterate since often due to Google API query limit)
distinctPlace1 <- distinctPlace2 %>% filter(is.na(distinctPlace2$lat))
b <- as_tibble(geocode(distinctPlace1$homeTown, output = "more") %>% 
			   select(lat, lon, postal_code)) 
distinctPlace1 <- bind_cols(distinctPlace1 %>% select(homeTown),b)

distinctPlace2 <- left_join(distinctPlace2, distinctPlace1, by="homeTown") %>% 
				   transmute(homeTown, lat = ifelse(is.na(lat.y), lat.x, lat.y),
									   lon = ifelse(is.na(lon.y), lon.x, lon.y),
									   postal_code = ifelse(is.na(postal_code.y), postal_code.x, postal_code.y)) 

write.csv(distinctPlace2, file = "latLong.csv", row.names = FALSE)			   

# Alternative method for obtaining postcodes from lat. and lon.
postcodes <- revgeo(distinctPlace$lon, distinctPlace$lat, output="frame") %>% select(zip)
postcodes <- bind_cols(distinctPlace, postcodes)
write.csv(postcodes, file = "postcodes.csv", row.names = FALSE)

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
						  
hospData %>% group_by(code) 
            
hospitals <- hospData %>% group_by(code) %>%
						  summarise(nCases = n())
  
## Histogram of cases per hospital
ggplot(data = hospitals, aes(x=nCases)) +
  geom_bar() +
  labs(	x = "", y = "Count") +
  theme(legend.position="none")# 
  
## Subset of hospitals with some record of 'hospital exposure'
theseHospitals <- hospData %>% filter(isHCW==1, !duplicated(code)) %>% select(code)
							   #filter(hospitalExposure == 1 | !is.na(exp1Name) | isHCW==1, !duplicated(code)) %>% select(code)

## Plot time series of hospital outbreak by population group
#thisHosp <- hospData %>% filter(code == "LARI-PNOS", !duplicated(ourID))  %>% select(ourID, onsetDate, isRoma, isHCW, popGroup, exp1Name, HCWHospName, RnTime, RnPlace, RnCombined)# %>% print(n=60)						  
theseHosps <- semi_join(hospData, theseHospitals) %>%
		      select(ourID, code, onsetDate, isRoma, isHCW, popGroup, hospitalExposure, exp1Name, HCWHospName, RnTime, RnPlace, RnCombined)# %>% print(n=60)	
						 
theseHosps$popGroup <- factor(theseHosps$popGroup) 

ggplot(data = theseHosps, aes(x=onsetDate, fill=popGroup) ) +
 geom_bar(width= 0.99) +
  labs(	x = "Date of symptom onset", y = "Count") 

# Bubble plot to show population group and number of attributable infections
ggplot(data = theseHosps, aes(x = onsetDate, y = code)) + 
  geom_point(aes(fill = popGroup, size = RnCombined, colour = (hospitalExposure==1 | !is.na(exp1Name))), pch=21, alpha = 0.8) +
  scale_colour_manual(name = 'Hospital exposure proposed', values = setNames(c('black','white'),c(T, F))) +
  scale_size(range = c(1, 10))  # Adjust the range of points size
    
#####################
## Useful snippets ##
#####################


## Convert municipality names from Greek to Latin script (needs to be sense-checked)
latinNames <- stri_trans_general(mapMunc$LEKTIKO, "Greek-Latin")
write.csv(latinNames, file = "muncNames.csv", row.names = FALSE)	
write.csv(mapMunc$KALCODE, file = "muncCodes.csv", row.names = FALSE)	
    