# Create a map with two different shapefiles laid over each other #
# Sept 2016 #
# Stephanie Kirmer #
# skirmer@uchicago.edu #

library(rgeos)
library(rgdal)
library(ggplot2)
library(ggmap)
library(dplyr)
library(sp)
library(spdep)
library(raster)
library(RColorBrewer)
library(maptools)
library(classInt)
library(zipcode)
library(broom)

# This script is designed to walk you through building a map that overlays two different groupings of data on the same map.
# It assumes you are starting with one normal shapefile, and one dataset without coordinates that is going to be combined with a shapefile
# to add coordinates that do not correspond 1:1 with your first shapefile.

# The example I am using is neighborhoods in the city of Chicago and income according to zipcodes in the City of Chicago. Data is
# from the US Census, the IRS, and the City of Chicago Data Portal.


#1. Load first shapefile (Neighborhoods)
shapefile1 <- readOGR(dsn="/export/home/skirmer/R/Spatial/shapefiles/Boundaries_-_Neighborhoods/Neighborhoods_2012b.shp", layer="Neighborhoods_2012b", verbose=F)

head(shapefile1)
proj4string(shapefile1)

#2. Format for ggplot
shapefile1@data$id <- rownames(shapefile1@data) #create an id based on row number
shapefile1_df<- tidy(shapefile1, region="id") #fortify using the ID so things can be reconnected
shapefile1_df2 <- merge(shapefile1_df, shapefile1, by="id") #get the labels back using the ID as merge connection

head(shapefile1_df2)

#3. Test Run Plot
ggplot(shapefile1_df2, aes(x=long, y=lat, group=PRI_NEIGH))+
  theme(panel.background=element_rect(fill="white", color="black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())+
  geom_polygon(aes(x=long, y=lat, group=group), fill="light gray")+
  geom_path(aes(x=long, y=lat, group=group), color="black", size=0.2)+
  coord_quickmap()+
  labs(title="Chicago Neighborhoods")




#4. Load second file (file2) - is not a shapefile and has no lat/long coding
file2 <- read.csv("/export/home/skirmer/R/Spatial/shapefiles/zip_income_2014.csv")

#add lat/long from a shapefile that can overlay the first shapefile - there needs to be overlap, obviously, but not 1:1 groupings.
il_zips <- readOGR(dsn="/export/home/skirmer/R/Spatial/shapefiles/cb_2015_us_zcta510_500k/cb_2015_us_zcta510_500k.shp", layer="cb_2015_us_zcta510_500k", verbose=F)

#This is just for narrowing down to the same locations for the example.
data(zipcode)
head(zipcode[zipcode$city =="Chicago",])
zipcodes_to_use <- zipcode$zip[zipcode$state=="IL" & zipcode$city == "Chicago"]
il_zips <- il_zips[il_zips$ZCTA5CE10 %in% zipcodes_to_use,]

#5. Ensure the new shapefile shares projection with the first one
coords <- proj4string(shapefile1)
il_zips.trans <- spTransform(il_zips, CRS=coords)

#6. Combine the non-shapefile data with the new shapefile
il_zips.trans <- il_zips.trans[il_zips.trans$ZCTA5CE10 %in% zipcodes_to_use,]
il_zips.trans@data$id <- rownames(il_zips.trans@data) #create an id based on row number
il_zips_df<- tidy(il_zips.trans, region="id") #fortify using the ID so things can be reconnected
il_zips_df2 <- merge(il_zips_df, il_zips.trans, by="id") #get the labels back using the ID as merge connection
file2 <- merge(file2, il_zips_df2[,c("long", "lat", "ZCTA5CE10")], by.x="zipcode", by.y="ZCTA5CE10")

#7. Convert new file into a shapefile of its own
file2 <- SpatialPointsDataFrame(coords = file2[,c("long","lat")], data=file2, proj4string=CRS(proj4string(shapefile1))) #projection needs to match starting file

#8. Format for ggplot
file2@data$id = rownames(file2@data) #create an id based on row number
file2_df <- tidy(file2, region="id") #fortify using the ID so things can be reconnected
file2_df2 <- merge(file2_df, file2, by="id") #get the labels back using the ID as merge connection

#9. Summarize data according to grouping from first shapefile (only if the second file is not already grouped by that)
inczip <- summarize(group_by(file2_df2, zipcode.x),
                    totalincome = sum(A02650.x, na.rm=T),
                    returns = sum(N1.x, na.rm=T),
                    mean_income = totalincome/returns*1000)

inczip$zipcode <- as.factor(inczip$zipcode.x)
inczip$avg_cat[inczip$mean_income >= 200000] <- 6
inczip$avg_cat[inczip$mean_income < 200000 & inczip$mean_income >= 100000] <- 5
inczip$avg_cat[inczip$mean_income < 100000 & inczip$mean_income >= 75000] <- 4
inczip$avg_cat[inczip$mean_income < 75000 & inczip$mean_income >= 50000] <- 3
inczip$avg_cat[inczip$mean_income < 50000 & inczip$mean_income >= 25000] <- 2
inczip$avg_cat[inczip$mean_income <25000] <- 1
inczip <- as.data.frame(inczip)

income_zip_map <- merge(il_zips_df2, inczip, by.x="ZCTA5CE10", by.y="zipcode.x", all.x=T)
income_zip_map <- arrange(income_zip_map, order)

head(income_zip_map)

#10. Create labels for the map (not required of course)
names <- summarize(group_by(shapefile1_df2, PRI_NEIGH),
                    mean_x = mean(long),
                   mean_y = mean(lat))


#11. Run Plot
#color palette is from color brewer. Order of polygon arguments is important to ensure grey for no-data but not hiding the colors desired.

ggplot()+
  theme(panel.background=element_rect(fill="white", color="black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())+
  geom_polygon(data=shapefile1_df2, aes(x=long, y=lat, group=PRI_NEIGH), fill="light gray")+
  geom_polygon(data=income_zip_map, aes(long, lat, group=zipcode, fill=avg_cat))+
  geom_path(data=shapefile1_df2, aes(x=long, y=lat, group=PRI_NEIGH), color="black", size=0.2)+
  geom_text(data=names, aes(x=mean_x, y=mean_y, label=PRI_NEIGH), size=3, check_overlap=TRUE)+
  coord_equal()+
  annotate("text", x = 1122322, y = 1816675, label = "Colors are assigned to zip \ncodes, black line boundaries\n represent neighborhoods. \nText labels identify neighborhoods.")+
  scale_fill_distiller("Mean Income", type="div", palette="RdBu", direction = 1,
                       labels=c("Under $25k", "$25k-$50k", "$50k-$75k", "$75k-$100k", "$100k-$200k", "$200k+"))+
  labs(title="Illinois Neighborhoods Plus Zip Code Income Levels")

