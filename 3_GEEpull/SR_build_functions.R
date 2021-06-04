
##########################################################################
### Functions used in making and analyzing river surface reflectance database
###########################################################################

# Function to clean and filter surface reflectance pull
pull_clean <- function(df, dswe = 2, minPix = 5, maxcScore = 900, maxClouds = 50,
                       maxRef = 10000, minRef = 1) {
  
  data <- df %>%
    dplyr::filter_if(is.numeric,all_vars(!is.na(.))) %>%
    dplyr::filter( 
      pixelCount > minPix,
      Cloud_Cover < maxClouds,
      cScore < maxcScore,
      dswe <= dswe,
      hillshadow > 0, 
      blue <= maxRef,
      blue >= minRef,
      green <= maxRef,
      green >= minRef,
      red <= maxRef,
      red >= minRef,
      nir <= maxRef,
      nir >= minRef,
      swir1 <= maxRef,
      swir1 >= minRef,
      swir2 <= maxRef,
      swir2 >= minRef) %>%
    separate(date, into =c("date", "time"), sep="T", remove=T) %>%
    dplyr::mutate(
      sat = factor(sat, levels = c(8,7,5), labels = c( 'LC08','LE07', 'LT05')),
      date_time= ymd_hms(paste0(date, time, sep=" ")),
      date = ymd(date),
      month = as.numeric(month(date)),
      year = year(date),
      hour = hour(date_time)) %>%
    mutate(season = case_when(
      month %in%  9:11 ~ "Fall",
      month %in%  c(12,1,2)  ~ "Winter",
      month %in%  3:5  ~ "Spring",
      TRUE ~ "Summer")) %>%
    mutate(decade= cut(year, breaks = c(1983,1990,1995,2000,2005,2010,2015,2020),
                       labels = c(1990,1995,2000,2005,2010,2015,2020) )) %>%
    rename(LS_ID = 'system:index') %>% 
    dplyr::group_by(ID) %>%
    dplyr::mutate(count =n(),
                  max_year=max(year, na.rm=T),
                  min_year = min(year, na.rm=T),
                  n_years = (max_year - min_year)) %>%
    ungroup() 
  
  return(data)
}

##########################################################################
# Function that calculates band metrics and water color

pull_transform <- function(df, maxRGB=10000, RGB=F) {
  
  if(RGB == T) { 
    
    maxRGB <- maxRGB
    
    data <- df %>%
      filter_at(vars(red, green, blue), all_vars(.< maxRGB))
    
  }else{ 
    data <- df
    
    maxRGB <- df  %>%
      dplyr::select(red, green, blue) %>%
      dplyr::summarise(maxRGB = max(., na.rm=T)) 
    
    maxRGB <- maxRGB$maxRGB
  }
  
  data<- data %>%
    dplyr::mutate(NR = nir/red,
                  BR = blue/red,
                  GR = green/red,
                  SR = swir1/red,
                  BG = blue/green,
                  RG = red/green, 
                  NG = nir/green,
                  SG = swir1/green,
                  BN = blue/nir,
                  GN = green/nir,
                  RN = red/nir,
                  SN = swir1/nir,
                  BS = blue/swir1,
                  GS = green/swir1,
                  RS = red/swir1,
                  NS = nir/swir1,
                  R.GN = red/ (green + nir),
                  R.GB = red/ (green + blue),
                  R.GS = red/ (green + swir1),
                  R.BN = red/ (blue + nir),
                  R.BS = red/ (blue + swir1),
                  R.NS = red/ (nir + swir1),
                  G.BR = green/ (blue + red),
                  G.BN = green / (blue + nir),
                  G.BS = green / (blue + swir1),
                  G.RN = green / ( red + nir),
                  G.RB = green / (red + blue),
                  G.NS = green / (nir + swir1),
                  B.RG = blue / (red + green),
                  B.RN = blue / (red + nir),
                  B.RS = blue / (red + swir1),
                  B.GN = blue / (green + nir),
                  B.GS = blue / (green + swir1),
                  B.NS = blue / (nir + swir1),
                  N.RG = nir / (red + green),
                  N.RB = nir / (red + blue),
                  N.RS = nir / (red + swir1),
                  N.GB = nir / (green + blue),
                  N.GS = nir / (green + swir1),
                  N.BS = nir / (blue  + swir1),
                  
                  GR2 = (green + red) / 2,
                  GN2 = (green + nir) / 2,
                  #blooms
                  BR_G = (blue - red) / green,
                  NS_NR = (nir - swir1) / (red - swir1),
                  fai = nir - (red + (swir1-red)*((830-660)/(1650-660))),
                  # fai = (nir - red) + (red -swir) * (830-660)/(1648-660)
                  N_S= nir - swir1,
                  N_R = nir- red,
                  #
                  ndvi = ((nir-red)/(nir+red)),
                  ndwi = ((green- swir1)/(green + swir1)),
                  ndssi = ((blue - nir)/ (blue + nir)),
                  gn.gn= ((green- nir)/ (green + nir)),
                  hue = rgb2hsv(r=red, g=green, b=blue, maxColorValue = maxRGB)[1,],
                  saturation = rgb2hsv(r=red, g=green,  b=blue, maxColorValue = maxRGB)[2,],
                  bright = rgb2hsv(r=red, g=green,  b=blue, maxColorValue = maxRGB)[3,],
                  bright_tot = (red + green + nir +blue),
                  dw = chroma(R=red, G=green, B=blue),
                  hexcolor = rgb(r=red, g=green, b=blue, maxColorValue = maxRGB)) 

  return(data)
}

#################################################################

# Function used in function above to calculate dominant wavelength in chromaticity 
# colorspace.

chroma <- function(R, G, B) {
  require(colorscience)
  require(tidyverse)
# Converst R,G, and B spectral reflectance to dominant wavelength based
# on CIE chromaticity color space

# see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
# Classification of Inland Water With the Forel-Ule
# Scale: A Case Study of Lake Taihu

# chromaticity.diagram.color.fill()
Xi <- 2.7689*R + 1.7517*G + 1.1302*B
Yi <- 1.0000*R + 4.5907*G + 0.0601*B
Zi <- 0.0565*G + 5.5943*B

x <-  Xi / (Xi + Yi +  Zi)
y <-  Yi / (Xi + Yi +  Zi)
z <-  Zi / (Xi + Yi +  Zi)

# calculate hue angle
alpha <- atan2( (x - (1/3)), (y - (1/3))) * 180/pi

# make look up table for hue angle to wavelength conversion
cie <- cccie31 %>%
  dplyr::mutate(a = atan2( (x - (1/3)), (y - (1/3))) * 180/pi) %>%
  dplyr::filter(wlnm <= 700) %>%
  dplyr::filter(wlnm >=380) 

# find nearest dominant wavelength to hue angle
wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']

return(wl)
}

##############################################################
# Function to add column of chromaticity coordinates and hue angle to dataframe
hueangle <- function(df, R, G, B) {
  
  require(colorscience)

  R = as.data.frame(df)[,R]
  G = as.data.frame(df)[,G]
  B = as.data.frame(df)[,B]
  
  # Add hue angle
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  chroma_x <-  Xi / (Xi + Yi +  Zi)
  chroma_y <-  Yi / (Xi + Yi +  Zi)
  chroma_z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  hue_angle <- atan2( (chroma_x - (1/3)), (chroma_y - (1/3))) * 180/pi
  
  return(cbind(df, chroma_x, chroma_y, chroma_z, hue_angle))
}



#####################################################################
# Function to calculate purity in chromaticity colorspace. Purity is similar
# to saturation is HSV colorspace

purity <- function(R, G, B) {
  require(colorscience)
  # Converst R,G, and B spectral reflectance to dominant wavelength based
  # on CIE chromaticity color space
  
  # see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
  # Classification of Inland Water With the Forel-Ule
  # Scale: A Case Study of Lake Taihu
  
  
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2( (x - (1/3)), (y - (1/3))) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    dplyr::mutate(a = atan2( (x - (1/3)), (y - (1/3))) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380) #%>%

  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']

  wl_xy <- wl %>%
    as_data_frame() %>%
    left_join(cie %>%
    dplyr::filter(wlnm %in% wl), by= c("value"= "wlnm"))
  
  xd <- wl_xy %>% pull(x)
    
  yd <- wl_xy %>% pull(y)
  
  P <- sqrt((x - (1/3))^2 + (y - (1/3))^2) / sqrt((xd - (1/3))^2 + (yd - (1/3))^2) 
  
  return(P)
}


#########################################################
# 
ls_correction <- function(data) {

# bands are quite different but dw is not 
sr_57 <- data %>%
  filter(sat %in% c("LE07", "LT05")) %>%
  filter(between(date, "1999-01-01", "2012-05-01" )) %>%
  # filter to site with enough data
  filter(n_years > 10) %>%
  select(ID, date, sat, count, n_years, blue, red, green, nir, swir1, swir2) %>%
  gather(blue:swir2, green, key='band', value='value') 

# plot distribution compare means
# ggplot(data=sr_57, aes(x=value, color=sat)) +
#   geom_density() +
#   facet_wrap(~band, scales="free")


# do ranking plotting percentiles, joining, and correcting
sr_57_rank  <- sr_57 %>%
  droplevels() %>%
  filter(sat =="LT05") %>%
  group_by(band) %>%
  nest() %>%
  mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
          ret = purrr::invoke_map(tibble, ret)) %>%
  unnest(ret) %>%
  dplyr::select(-data) %>%
  pivot_longer(
    cols= contains("%")
  ) %>%
  mutate(quant = parse_number(name)/100) %>%
  rename(value_5 = value) %>%
  inner_join(sr_57 %>%
               droplevels() %>%
               filter(sat =="LE07") %>%
               group_by(band) %>%
               nest() %>%
               mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
                       ret = purrr::invoke_map(tibble, ret)) %>%
               unnest(ret) %>%
               dplyr::select(-data) %>%
               pivot_longer(
                 cols= contains("%")
               ) %>%
               mutate(quant = parse_number(name)/100) %>%
               rename(value_7 = value) %>%
               dplyr::select(-name),
             by=c("band", "quant")) 

poly_5_trunc <- function(df){
  lm(value_7 ~ poly(value_5, 2, raw=T), data = df %>%
       filter(!quant %in% c(0, 1))  )
}
poly_5_all <- function(df){
  lm(value_7 ~ poly(value_5, 2, raw=T), data = df)
}

## polynomial correction fit
poly_57 <- sr_57_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  nest(-band) %>%
  mutate( model = purrr::map(data, poly_5_trunc)) %>%
  mutate( model_all = purrr::map(data, poly_5_all)) %>%
  mutate( pred = purrr::map2(model, data, predict)) %>%
  mutate( pred_all = purrr::map2(model_all, data, predict)) %>%
  unnest(c(pred, pred_all, data))  %>%
  dplyr::select(-model, -model_all)

png("D:/Dropbox/projects/rivercolor/figs/poly_5_correction.png",
    units="in", width = 6, height=4, res = 300)

ggplot(poly_57) +
   geom_point( aes(x=value_5, y= value_7))+
   geom_point( aes(x=pred_all, y= value_7), color="red", alpha=0.5)+
   geom_point( aes(x=pred, y= value_7), color="blue", alpha=0.5)+
   geom_abline(aes(slope=1, intercept=0)) +
   facet_wrap(~band, scales="free") +
  theme_bw() 

dev.off()
#ggsave("figs/poly_5_correction.png", units="in", width = 6, height=4, dpi = 300)

#
write_csv(poly_57, "D:/Dropbox/projects/rivercolor/out/ls_57_quantile_correction.csv")

##
coef_5 <- sr_57_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  filter(!quant %in% c(0, 1)) %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_5, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LT05") %>%
  mutate(fit = "98_quant")

coef_5_all <- sr_57_rank %>%
  ungroup() %>%
  #  filter(band != "dw") %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_5, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LT05") %>%
  mutate(fit = "all_quant")

########################

sr_78 <- data %>%
  filter(sat %in% c("LE07", "LC08")) %>%
  filter(date > "2013-04-11" ) %>%
  # filter to site with enough data
  filter(n_years > 10) %>%
  select(ID, date, sat, count, n_years, blue, red, green, nir, swir1, swir2) %>%
  gather(blue:swir2, key='band', value='value') 

# plot distribution compare means
# ggplot(data=sr_78, aes(x=value, color=sat)) +
#   geom_density() +
#   facet_wrap(~band, scales="free")

# do ranking plotting percentiles, joining, and correcting
sr_78_rank  <- sr_78 %>%
  droplevels() %>%
  filter(sat =="LC08") %>%
  group_by(band) %>%
  nest() %>%
  mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
          ret = purrr::invoke_map(tibble, ret)) %>%
  unnest(ret) %>%
  dplyr::select(-data) %>%
  pivot_longer(
    cols= contains("%")
  ) %>%
  mutate(quant = parse_number(name)/100) %>%
  rename(value_8 = value) %>%
  inner_join(sr_78 %>%
               droplevels() %>%
               filter(sat =="LE07") %>%
               group_by(band) %>%
               nest() %>%
               mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
                       ret = purrr::invoke_map(tibble, ret)) %>%
               unnest(ret) %>%
               dplyr::select(-data) %>%
               pivot_longer(
                 cols= contains("%")
               ) %>%
               mutate(quant = parse_number(name)/100) %>%
               rename(value_7 = value) %>%
               dplyr::select(-name),
             by=c("band", "quant"))  

# ggplot(sr_78_rank, aes(value_8, value_7)) +
#   geom_point() +
#   geom_abline(aes(slope=1, intercept=0), color="red") +
#   geom_smooth(method="lm", color="blue") +
#   facet_wrap(~band, scales = "free")

poly_8_trunc <- function(df){
  lm(value_7 ~ poly(value_8, 2), data = df %>%
       filter(!quant %in% c(0, 1))  )
}
poly_8_all <- function(df){
  lm(value_7 ~ poly(value_8, 2), data = df)
}

poly_78 <- sr_78_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  nest(-band) %>%
  mutate( model = purrr::map(data, poly_8_trunc)) %>%
  mutate( model_all = purrr::map(data, poly_8_all)) %>%
  mutate( pred = purrr::map2(model, data, predict)) %>%
  mutate( pred_all = purrr::map2(model_all, data, predict)) %>%
  unnest(c(pred, pred_all, data)) %>%
  dplyr::select(-model, -model_all)

png("D:/Dropbox/projects/rivercolor/figs/poly_8_correction.png",
    units="in", width = 6, height=4, res = 300)

ggplot(poly_78) +
  geom_point( aes(x=value_8, y= value_7))+
  geom_point( aes(x=pred_all, y= value_7), color="red", alpha=0.5)+
  geom_point( aes(x=pred, y= value_7), color="blue", alpha=0.5)+
  geom_abline(aes(slope=1, intercept=0)) +
  facet_wrap(~band, scales="free") +
  theme_bw()
dev.off()

#ggsave("figs/poly_8_correction.png", units="in", width = 6, height=4, dpi = 300)

#
write_csv(poly_78, "D:/Dropbox/projects/rivercolor/out/ls_78_quantile_correction.csv")

#
coef_8 <- sr_78_rank %>%
  ungroup() %>%
  filter(!quant %in% c(0, 1)) %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_8, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LC08") %>%
  mutate(fit = "98_quant")

#
coef_8_all <- sr_78_rank %>%
  ungroup() %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_8, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LC08") %>%
  mutate(fit = "all_quant")

#
coef_7 <- tibble(band = c("blue", "red", "green", "nir", "swir1", "swir2"), intercept = 0, coef1=1, coef2=0, sat= "LE07")

#
corr_coef <- bind_rows(coef_5, coef_7, coef_8, coef_5_all, coef_8_all) %>%
  ungroup()

return(corr_coef)

}
