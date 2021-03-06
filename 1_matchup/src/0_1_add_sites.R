
library(LAGOSNE)
library(tidyverse)
library(feather)
library(sf)
library(maps)
library(magrittr)
library(purrr)
library(data.table)
library(ggthemes)
library(dplyr)
library(dataRetrieval)
library(ggplot2)

##################################################################################
#### Combiningin-situ  datasets from USGS sediment and with Aquasat
###############################################################################

wqp.lagos <- read_feather('D:/Dropbox/projects/TSS/data/wqp_lagos_unity.feather')

## add stae county code then convert time local to GMT
usgs_site <- read_csv("data/discrete_sites.csv") %>%
  mutate(SITE_NO = as.character(SITE_NO))

## load unique SSC sites combined usgs and wqp
usgs.tss <- read_csv("data/USGS_tis.csv") %>%
  mutate(SITE_NO = as.character(SITE_NO),
         dt = as.POSIXct(DATETIME, format="%m/%d/%Y %H:%M", tz="GMT") ,
         date =as.Date(dt) ) %>%
  left_join(usgs_site, by="SITE_NO") %>%
  filter(SSC > 0, SSC < 10^5, !is.na(SSC)) %>%
  mutate(SiteID = paste("USGS-",SITE_NO, sep="")) %>%
  rename(tss = SSC, date_unity = dt,  MonitoringLocationName=STATION_NM,HUCEightDigitCode = HUC_8, lat=LATITUDE, long = LONGITUDE, StateCode = STATE,
         CountyCode = COUNTY_NAME, resultCount=sample_count) %>%
  mutate(MonitoringLocationTypeName = "Stream",
         ResolvedMonitoringLocationTypeName = "Stream",
         CountryCode="US",
         source="USGSsed",
         Constituent = "tss",
         HorizontalCoordinateReferenceSystemDatumName="NAD83",
         HUCEightDigitCode =  as.character(HUCEightDigitCode),
         date_only=FALSE,
         epsg = 4269) %>%
  select(SiteID, date_unity, date_only, tss, source)
  
wqp.lagos.usgs <- bind_rows(wqp.lagos, usgs.tss)

#write_feather(wqp.lagos.usgs, "output/wqp_lagos_usgs_unity.feather")

#############################################################
## combine new usgs data, and lagos with wqp_long_with_methods

wqp_long <- read_csv("D:/Dropbox/projects/TSS/data/wqp_long_with_methods.csv") %>%
  mutate(date_unity = ymd_hms(ifelse(is.na(date_time),
                                     paste(date,'00:00:00'),
                                     as.character(date_time)),
                                     tz='UTC'),
         source = "WQP") %>%
  mutate(date = as.character(date),
         time = as.character(time)) %>%
  select(-date_time)


# load new usgs sites and convert to long format 
usgs_long <- read_csv("data/USGS_tis.csv") %>%
  filter(SSC > 0, SSC < 10^5, !is.na(SSC)) %>%
  rename(harmonized_value = SSC) %>%
  mutate(SiteID = paste("USGS-",SITE_NO, sep=""),
         date_unity = as.POSIXct(DATETIME, format="%m/%d/%Y %H:%M", tz="UTC")) %>%
  mutate(
         date = as.character(format(date_unity, '%Y-%m-%d')),
         time = as.character(format(date_unity, '%H:%M:%S')),
         source="USGSsed",
         date_only=FALSE,
         characteristicName = "Suspended Sediment Concentration (SSC)",
         harmonized_parameter = "tss" ) 


cols <- intersect(names(wqp_long), names(usgs_long))

usgs_long <- usgs_long %>%
  select(cols)

## join lagos sites to wqp_lagos_unity and convert to long
lagos_long <- wqp.lagos %>%
  filter(source =="LAGOS") %>%
  mutate(date = as.character(format(date_unity, '%Y-%m-%d'))) %>%
  mutate(date_unity = paste(date,'00:00:00')) %>%
  mutate(date_unity = ymd_hms(ifelse(date_only == TRUE,
                                     paste(date,'00:00:00'),
                                     as.character(date_unity)),
                              tz='UTC')) %>%

  rename(chl.a = chl_a, p.sand = p_sand) %>%
  gather(chl.a:tss, key='harmonized_parameter', value='harmonized_value') %>%
  filter(!is.na(harmonized_value)) %>%
  mutate(characteristicName = ifelse(harmonized_parameter=="chl.a", "Chlorophyll a", NA),
         characteristicName = ifelse(harmonized_parameter=="doc", "Organic carbon", characteristicName),
         characteristicName = ifelse(harmonized_parameter=="secchi", "Water transparency, Secchi disc", characteristicName),
         characteristicName = ifelse(harmonized_parameter=="tis", "Fixed suspended solids", characteristicName),
         characteristicName = ifelse(harmonized_parameter=="tss", "Total suspended solids" , characteristicName),
         characteristicName = ifelse(harmonized_parameter=="p.sand", "percent sand" , characteristicName),
         time = as.character(format(date_unity, '%H:%M:%S')),
         source = "LAGOS"
  )
  
cols2 <- intersect(names(wqp_long), names(lagos_long))

lagos_long <- lagos_long %>%
  select(cols2)

wqp_lagos_usgs_long <- bind_rows(wqp_long, usgs_long, lagos_long)

#write_feather(wqp_lagos_usgs_long, "output/wqp_lagos_usgs_long_methods.feather")

#### make inventory data

#Read in inventory and exclude nutrient data (not a goal of this project just yet)
inv_simple <- #read_feather('1_wqdata/out/wqp_inventory.feather') %>%
  read_feather("data/wqp_inventory.feather") %>%
  rename(SiteID=MonitoringLocationIdentifier,lat=LatitudeMeasure,
         long=LongitudeMeasure,orig_datum = HorizontalCoordinateReferenceSystemDatumName,
         obs_count=resultCount) %>% 
  filter(!is.na(lat)) %>% #Exclude sites with missing geospatial data
  filter(!is.na(long)) %>%
  mutate(source='WQP')

#write_feather(inv_simple,'1_wqdata/out/inv_simple.feather')

#Make a table of all the projections
table(inv_simple$orig_datum)

#Setup a datum table to transform different datums into WGS84 for GEE and to 
#match LAGOS. Note that we make the somewhat reisky decision that others and unknowns
# are mostly NAD83 (since > 60% of the data is originally in that projection)
datum_epsg <- tibble(orig_datum=c('NAD27','NAD83','OTHER','Unknown','UNKWN','WGS84','WGS72'),
                     epsg=c(4267,4269,4269,4269,4269,4326,4322))

# Get distinct lat longs and sites
inv_uniques <- inv_simple %>%
  distinct(SiteID,lat,long,orig_datum, .keep_all = T)

## project inventory using the above lookup table
## Have to loop over each projection
projected <- list() 
for(i in 1:nrow(datum_epsg)){
  #Select each datum iteratively
  d1 <- datum_epsg[i,]
  # Join inventory to that datum only and then feed it that CRS and then 
  # transform to wgs84
  inv_reproj <- inv_uniques %>%
    filter(orig_datum == d1$orig_datum) %>%
    inner_join(d1, by='orig_datum') %>%
    st_as_sf(.,coords=c('long','lat'),crs=d1$epsg) %>%
    st_transform(4326)
  
  projected[[i]] <- inv_reproj
}
# Note that we lose some sites in GUAM which is intentional
# Add back in lat longs which are all WGS84 now
inv_wgs84 <- do.call('rbind',projected) %>%
  mutate(lat = st_coordinates(.)[,2],
         long=st_coordinates(.)[,1]) %>%
  as.data.frame(.) %>%
  select(-geometry,-orig_datum) %>%
  mutate(source='WQP')

#Load in lagos data
cols <- intersect(names(inv_wgs84), names(lagos.sites))

lagos.sites <- lagosne_load("1.087.1")$locus %>%
  rename(SiteID = lagoslakeid,lat=nhd_lat,long=nhd_long,
         MonitoringLocationName = gnis_name) %>%
  mutate(SiteID = paste('lagos',SiteID,sep='-'),
         source='LAGOS',
         CountryCode = "US",epsg=4326)  %>%
  separate(county_zoneid, into = c('junk', 'CountyCode'), sep="_", remove=T ) %>%
  separate(state_zoneid, into = c('junk2', 'StateCode'), sep="_", remove=T ) %>%
  mutate(StateCode = str_pad(StateCode, 2, pad="0"),
         CountyCode = str_pad(CountyCode, 3, pad="0"),
         ResolvedMonitoringLocationTypeName = "Lake")

cols <- intersect(names(inv_wgs84), names(lagos.sites))

lagos.sites <- lagos.sites %>%
  select(cols)

#load in new usgs SSC data
library(cdlTools)
usgs_site <- read_csv("data/discrete_sites.csv") %>%
  mutate(SITE_NO = as.character(SITE_NO))

## load unique SSC sites combined usgs and wqp
usgs.inv <- read_csv("data/USGS_tis.csv") %>%
  mutate(SITE_NO = as.character(SITE_NO),
         dt = as.POSIXct(DATETIME, format="%m/%d/%Y %H:%M", tz="GMT") ,
         date =as.Date(dt) ) %>%
  left_join(usgs_site, by="SITE_NO") %>%
  filter(SSC > 0, SSC < 10^5, !is.na(SSC)) %>%
  mutate(SiteID = paste("USGS-",SITE_NO, sep="")) %>%
  rename(tss = SSC, date_unity = dt,  MonitoringLocationName=STATION_NM,HUCEightDigitCode = HUC_8, lat=LATITUDE, long = LONGITUDE, StateCode = STATE,
         CountyCode = COUNTY_NAME, resultCount=sample_count) %>%
  select(SiteID, date, tss, lat, long, date_unity, resultCount, 
         MonitoringLocationName, HUCEightDigitCode, StateCode, CountyCode ) %>%
  mutate(MonitoringLocationTypeName = "Stream",
         ResolvedMonitoringLocationTypeName = "Stream",
         CountryCode="US",
         source="USGSsed",
         Constituent = "tss",
         HorizontalCoordinateReferenceSystemDatumName="NAD83",
         HUCEightDigitCode =  as.character(HUCEightDigitCode),
         StateCode = fips(StateCode, to="FIPS"),
         StateCode = str_pad(StateCode, 2, pad="0"),
         date_only=FALSE,
         epsg = 4269) %>%
  select(-tss)

cols2 <- intersect(names(inv_wgs84), names(usgs.inv))

usgs.inv <- usgs.inv %>%
  select(cols2)

#Bind lagos and wqp data together
lagos.wqp.inv <- bind_rows(inv_wgs84, usgs.inv, lagos.sites) %>%
  distinct(SiteID, lat, long, source, .keep_all = T) %>%
  arrange(lat,long)


