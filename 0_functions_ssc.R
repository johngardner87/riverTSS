

#######################################################
### functions
#######################################################

RBIcalc <- function(Q){

  # Size
  myLen <- length(Q)
  # Add previous record in second column
  Qprev <- c(NA,Q[-myLen])
  # Create dataframe.
  myData <- as.data.frame(cbind(Q,Qprev))
  # delta (absolute)
  myData[,"AbsDelta"] <- abs(myData[,"Q"] - myData[,"Qprev"])
  # SumQ
  SumQ <- sum(myData[,"Q"],na.rm=TRUE)
  # Sum Delta
  SumDelta <- sum(myData[,"AbsDelta"], na.rm=TRUE)
  #
  RBIsum <- SumDelta / SumQ
  
  # Return RBI value for data submitted.
  return(RBIsum)
  #
} 
########################################

RBIcalc2 <- function(x, col){
  
  Q <- x %>%
    pull(col)
  
  # Size
  myLen <- length(Q)
  # Add previous record in second column
  Qprev <- c(NA,Q[-myLen])
  # Create dataframe.
  myData <- as.data.frame(cbind(Q,Qprev))
  # delta (absolute)
  myData[,"AbsDelta"] <- abs(myData[,"Q"] - myData[,"Qprev"])
  # SumQ
  SumQ <- sum(myData[,"Q"],na.rm=TRUE)
  # Sum Delta
  SumDelta <- sum(myData[,"AbsDelta"], na.rm=TRUE)
  
  RBIsum <- SumDelta / SumQ
  
  # Return RBI value for data submitted.
  return(RBIsum)
  
} 
##########################################

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
###########################################

mkmk <- function(x, col) {
  
  require(modifiedmk)
  x = x %>%
    pull(col)
  
  y <- mkttest(x)
  
  out = tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
               s = y["S"], n=length(x))
  
  return(out)
}
#########################################

mmk <- function(x, col) {
  
  require(modifiedmk)
  x = x %>%
    pull(col)
  
  y <- mmky(x)
  
  out = tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["new P-value"], p.old = y["old P.value"],
                z.new = y["Corrected Zc"], z.old = y["Original Z"], new.var = y["new.variance"],
               old.var = y["old.variance"],  n=length(x), N = y["N/N*"])
  
  return(out)
}
##################################################

mk_name <- function(x, col) {
  
  require(modifiedmk)
  x = x %>%
    pull(col)
  
  y <- mkttest(x)
  
  out <- tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
               s = y["S"], n=length(x)) %>%
    rename_at(vars(tau:n), function(x) paste0(x, "_",col))
  
  return(out)
}
####################################################

mk_name2 <- function(x, col) {
  
  require(modifiedmk)
  x = x %>%
    pull(col)
  
  y <- mkttest(x)
  
  out <- tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
                s = y["S"], n=length(x), fit=col)
  
  return(out)
}
###############################################################

powerf <- function(data) {
  
  mod <- lm(log10(data$tss_mean_prof) ~ log10(data$mean_flow) )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="power")
  
  return(out)
}
##################################################

linearf<- function(data) {
  
  mod <- lm(data$tss_mean_prof ~ data$mean_flow )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="linear")
  
  return(out)
}
###############################################

expf<- function(data) {
  
  mod <- lm(log(data$tss_mean_prof) ~ (data$mean_flow) )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared ,n=nrow(data), fit="exponential")
  
  return(out)
}
################################################

logf <- function(data) {
   
  mod <- lm(log(data$tss_mean_prof) ~ log(data$mean_flow) )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="logarithmic") 
  
  return(out)
}
######################################################

lin_sp <- function(data) {
  
  mod <- lm((data$mean) ~ data$Pthlngt)
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=broom::glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="linear")
  
  return(out)
}
###########################################################

exp_sp <- function(data) {
  
  mod <- lm(log(data$mean) ~ data$Pthlngt)
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=broom::glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="exponential")
  
  return(out)
}
###############################################################

pow_sp <- function(data) {
  
  mod <- lm((data$mean) ~ log10(data$Pthlngt))
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=broom::glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="power")
  
  return(out)
}
#################################################################

log_sp <- function(data) {
  
  mod <- lm(log(data$mean) ~ log(data$Pthlngt))
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=broom::glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=nrow(data), fit="logarithmic")
  
  return(out)
}

###############################################################
MK_sp <- function(data) {
  
  require(modifiedmk)
  
  x <- data %>%
    arrange(desc(Pthlngt)) %>%
    pull(mean)
  
  y <- mkttest(x)
  
  out <- tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
                s = y["S"], n=length(x))

  return(out)
}
#############################

MKwbm_dam <- function(data) {
  
  require(modifiedmk)
  
  x <- data %>%
    arrange(desc(Pathlength)) %>%
    pull(Qc_gL)
  
  y <- mkttest(x)
  
  out <- tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
                s = y["S"], n=length(x))
  
  return(out)
}
##################################

MKwbm <- function(data) {
  
  require(modifiedmk)
  
  x <- data %>%
    arrange(desc(Pathlength)) %>%
    pull(Qc_SemiPri)
  
  y <- mkttest(x)
  
  out <- tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
                s = y["S"], n=length(x))
  
  return(out)
}
###################################

lm_name <- function(data, x, y) {
  
  data <- data %>%
    as.data.frame() %>%
    filter(!is.na(x) & !is.na(y))
  
  mod <- lm(data[,y] ~ data[,x] )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=length(x), fit="linear")
  
  return(out)
}
##############################################

power_name <- function(data, x, y) {
  
  data <- data %>%
    as.data.frame() %>%
    filter(!is.na(x) & !is.na(y))
  
  mod <- lm(log10(data[,y]) ~ log10(data[,x]) )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=length(x), fit="power")
  
  return(out)
}
##################################################

exp_name <- function(data, x, y) {
  
  data <- data %>%
    as.data.frame() %>%
    filter(!is.na(x) & !is.na(y))
  
  mod <- lm(log(data[,y]) ~ data[,x] )
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=length(x), fit="exponential")
  
  return(out)
}
##########################################

log_name <- function(data, x, y) {
  
  data <- data %>%
    as.data.frame() %>%
    filter(!is.na(x) & !is.na(y))
  
  mod <- lm(data[,y] ~ log(data[,x]))
  
  out <- tibble(slope = coef(mod)[2], yint=coef(mod)[1], pvalue=glance(mod)$p.value,
                r2=summary(mod)$r.squared, n=length(x), fit="logarithmic")
  
  return(out)
}
#####################################