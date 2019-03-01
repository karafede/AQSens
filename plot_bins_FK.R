
rm(list = ls())

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(broom)
library(minpack.lm)
library(mgcv)
library(threadr)
library(dygraphs)
library(ggpmisc)

# saving the original graphical parameters
op <- par(no.readonly = TRUE)
# Restoring graphical parameters on exit of function, even if an error occurs
# on.exit(par(op)) # it reset the par(mfrow) allways plotting on the upper left plot
par(mfrow=c(2,2))

# load bin counts from OPC-N3
# setwd
setwd("L:/ERLAP/Diffusion/Box Sync/AirSensEUR/Fieldtests/Shiny/JRC_11/General_data")
# setwd(choose.dir())
# load fitting functions repository
source("L:/ERLAP/Diffusion/Box Sync/AirSensEUR/Fieldtests/Shiny/151016 Sensor_Toolbox.R")
source("L:/ERLAP/Diffusion/Box Sync/AirSensEUR/Fieldtests/Shiny/Functions4ASE.R")
# source(choose.files())

#load assignation of bins for the OPCN3
bins_diameters <- read.csv("Bins_OPCN3_config.csv")

#########################################
# load DMPS and APS data ################
#########################################
# The counts of the OPC-N3 and PM10, PM2.5 ..., temperature and humidity are in General.Rdata.
# When loaded, it creates a dataFrame: General.df with all AirSensEUR sensor data
load("General.Rdata")

# The reference counts of the DMPS and APS are in RefData.Rdata
# It creates a dataframe RefData when loaded with all reference data
load("RefData.Rdata")
RefData <- RefData %>%
    filter(date %in% General.df$date)

#########################################
## Quick dynamic time series ############
#########################################
# Build timeseries for plots
time_series_sensor <- data_frame_to_timeseries(General.df %>%
                                                   select(date,
                                                          Particulate_Matter_1,
                                                          Particulate_Matter_25,
                                                          Particulate_Matter_10))

ts_sensor <- cbind(time_series_sensor$Particulate_Matter_25, 
                      time_series_sensor$Particulate_Matter_10,
                      time_series_sensor$Particulate_Matter_1)

plot_sensors <- dygraph(ts_sensor) %>% 
    dySeries("..1",label = "PM2.5", color = "red") %>%
    dySeries("..2",label = "PM10", color = "blue") %>%
    dySeries("..3",label = "PM1", color = "black") %>%
    dyAxis("y", label = "PM<sub></sub> (&#956;g m<sup>-3</sup>)") %>% 
    dyRangeSelector()
plot_sensors


#############################################################################
### Fitting the overlap between DMPS and APS darta to find density value ####
#############################################################################

### input parameters to be inserted by the USER
# Start and end date for selection of reference data
Start_Time <- as.POSIXct("2018-09-30 01:00", tz = "UTC")
Stop_Time   <- as.POSIXct("2018-10-07 01:00", tz = "UTC") #"2018-10-07"  #"2018-09-30 01:00"

#### TO BE integrated in the Function4ASE.R ################################

Final_function_OPC <- function(Begin, End, 
                               RefData = RefData, General.df = General.df, Sensor_dates = NULL, 
                               verbose = TRUE, Min_APS = 0.62, Max_APS =  0.800) {
    
    # Refdata dataframe with timestamp and reference counts
    # DateBegin,   The start and ending selected dates in POSIXCt format. Default values are NULL, 
    # DateEnd      It is possible to have only one date limit, DateBegin or DateEnd. DateBegin or DateEnd are not discarded.
    # Sensor_dates  vector of POSIXCt: the dates for which the PM sensor provides binned counts. 
    #               In case the sensor PM data series is not complete, dates with missing PM sensor data will be discared from RefData
    
    # Min_APS,      Numeric, Minimum and maximum diameters of reference counts to be selected 
    # Max_APS       for reference counts. Default values are NULL. It is possible to have only one diameter limit,  Min_APS or Max_APS.
    #               Min_APS and Max_APS are discarded.
    
    # gather data in list and dataframe
    # Compute mean of reference count
    # browser()
    Dist_Ref.d_TS <- Distribution_Ref_TS(RefData = RefData, Begin, End, 
                                         Min_APS = Min_APS, Max_APS = Max_APS) #
    
    # Determining of density with 5 last DMPS diameter, and 3rd to 6th APS counts
    # Diameters overlapping of DMPS and APS
    Diam_APS_to_keep  <- Dist_Ref.d_TS$Diam_APS
    Diam_DMPS_to_keep <- Dist_Ref.d_TS$Diam_DMPS[(length(Dist_Ref.d_TS$Diam_DMPS)-3):length(Dist_Ref.d_TS$Diam_DMPS)]
    
    DataXY.DMPS <- Dist_Ref.d_TS$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep)
    DataXY.DMPS <- data.frame(x = log10(DataXY.DMPS$diameters),
                              y = log10(DataXY.DMPS$counts))
    Model <- lm(y ~ x, data = DataXY.DMPS, model = TRUE, x = TRUE, y = TRUE)
    Model.i <- list(Tidy = tidy(Model), 
                    Augment = augment(Model), 
                    Glance = glance(Model), 
                    Call = Model$call, Coef = coef(Model),
                    Model = Model)
    
    DataXY.APS <- Dist_Ref.d_TS$counts_Ref %>% filter(diameters %in% Diam_APS_to_keep)
    DataXY.APS <- data.frame(x = log10(DataXY.APS$diameters),
                             y = log10(DataXY.APS$counts))
    
    # find optima value for the "density" of the particles
    density_fun <- optimise(f_Error, c(0.7, 3), DataXY.APS = DataXY.APS, Model.i = Model.i) # limit ranges for the density
    if (verbose) {
        
        Plot.DMPS       <- Plot_Dist_Ref(Dist_Ref.d_TS$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep), Model.i = Model.i)
        Plot.DMPS.APS   <- Plot.DMPS + geom_point(data = DataXY.APS, 
                                                      aes(x, y), stat = "identity", col = "red")
        Plot.DMPS.APS.d <- Plot.DMPS.APS + geom_point(data = DataXY.APS, 
                                                  aes(log10(10^x/density_fun$minimum), y), stat = "identity", shape = 17, col = "red", size = 3)
        
        Plot.DMPS.APS.d
    }
    
    

    ############################################################################################################################################ 
    ####### make difference MinAPS/MaxAPS for density determination and for fitting
    ############################################################################################################################################
    ### we have the correct density in Model.i and the correct diameters of APS in Dist_ref
    # return a list with "RefData_filtered" "counts_Ref"   "Diam_DMPS"         "Diam_APS" with the full diameters
    Dist_Ref.full <- Distribution_Ref_TS(RefData = RefData, DateBegin = Start_Time, DateEnd = Stop_Time, 
                                         Min_APS = 0.777, Max_APS =  14.86000) #
    
    ## correction of APS diamters with density in Dist_Ref.d determined in Model.i  <- Density_Ref_log
    Dist_Ref.full$counts_Ref$diameters <- Diam_APS_Corr(Dist_Ref.full$Diam_APS, Dist_Ref.full$counts_Ref$diameters, density_fun)
    if (verbose) Plot_Dist_Ref_log(Dist_Ref.full$counts_Ref)
    
    # Fitting of distribution of reference counts versus diameter in log/log with fitted density (using corrected diameters)
    Model.i.Gam  <- Density_Ref_log(Dist_Ref.full$counts_Ref, Density = NULL, Mod_type = 'GAM_GAUSS' )
    
    # Selected timestamp of OPC counts, make mean calculation
    # Change names of OPC Sensor with diameters and return Refdata with names of Bin, diamters and counts
    Dist_OPC_Sensor <- Distribution_OPC_Sensor(General.df = General.df, DateBegin = Start_Time, DateEnd = Stop_Time, 
                                               bins_diameters = bins_diameters) 
    
    # Correction of OPS diameter with density
    # use of Gam distribution to predict the reference counts at the OPC diameters
    # return predicted values of the OPC
    Predict_Dist_OPC_Sensor <- Density_OPC_predict(Dist_OPC_Sensor$counts_OPC, Mod_type = 'GAM_GAUSS', 
                                                   density = density_fun$minimum, Model.i  = Model.i.Gam)
    
    # Put together timeInterval, density, predict of reference counts at OPC corrected diameters, OPC counts,  OPCHUM, Temp, VOL, sampling rate
    
    # find Bin assignation for each diameter according to the Max APS diameter
    Bins <- as.data.frame(t(bins_diameters))
    Bins$Bins <- as.character(rownames(Bins))
    rownames(Bins) <- NULL
    names(Bins) <- c("diameters", "Bins")
    # Bins <- Bins %>%
    #     filter(diameters < max(Dist_Ref.full$Diam_APS, na.rm = T))
    

    # Determine the Bins used by the OPC and mantain their order (missing Bins will be replaced by NA)
    # Dist_OPC_Sensor$counts_OPC$ID <- 1: nrow(Dist_OPC_Sensor$counts_OPC)
    Bins$ID <- 1: nrow(Bins)
    
    OPC_all <- cbind(Dist_OPC_Sensor$counts_OPC,
                 Predict_Dist_OPC_Sensor$predict_counts)
    names(OPC_all) <- c("diameters", "counts", "pred_counts")

    OPC_all <- Bins %>%
        left_join(OPC_all, by = "diameters") %>%
        select(- diameters)
    
    df_OPC <- as.data.frame(t( round(OPC_all$counts, digits = 4)))
    names(df_OPC) <- paste(c(OPC_all$Bins))
    
    df_OPC_pred <- as.data.frame(t (round(OPC_all$pred_counts, digits = 4)))
    names(df_OPC_pred) <- paste0(c(OPC_all$Bins), "_pred")
    

    # data frame (TS) ###
    df_OPC_Ref <- cbind(Start_Time, Stop_Time, 
                        density = density_fun$minimum, 
                        df_OPC_pred, 
                        df_OPC,
                        Dist_OPC_Sensor$Met_OPC)
    
    return(df_OPC_Ref)

}

interval <- 60  # minutes
TS <- seq(from = Start_Time, by = interval*60, to = Stop_Time)

#Determining number of columns in the Final dataframe
df_OPC1 <- Final_function_OPC(RefData, General.df, Begin = Start_Time, End = Stop_Time, 
                              Sensor_dates = NULL, Min_APS = 0.62, Max_APS =  0.800, verbose = TRUE)
Columns <- length(names(df_OPC1))
Rows    <- length(TS) -1
df_OPC  <- as.data.frame(matrix( rep(0, Columns * Rows), nrow = Rows, ncol = Columns))
names(df_OPC) <- names(df_OPC1)
rm(df_OPC1)


# loop over DateBegin and DateEnd
for (i in 1:(length(TS)-1)) {
    
    Start_Time <- TS[i]
    Stop_Time   <- TS[i+1]
    
    df_OPC[i,] <- Final_function_OPC(RefData = RefData, General.df = General.df, Begin = Start_Time, End = Stop_Time, 
                                     Sensor_dates = NULL, Min_APS = 0.62, Max_APS =  0.800, verbose = TRUE)

}

df_OPC$Start_Time <- as.POSIXct(df_OPC$Start_Time, tz = "UTC", origin = strptime("1970-01-01", "%Y-%m-%d"))
df_OPC$Stop_Time   <- as.POSIXct(df_OPC$Stop_Time  , tz = "UTC", origin = strptime("1970-01-01", "%Y-%m-%d"))


########################################################################
########################################################################

# Start_Time <- as.POSIXct("2018-09-30 01:00", tz = "UTC")
# Stop_Time   <- as.POSIXct("2018-10-07 01:00", tz = "UTC") #"2018-10-07"  #"2018-09-30 01:00"
# 
# library(purrr)
# DateBegin_list <- TS 
# DateEnd_list <- seq(from = Start_Time +interval*60, by = interval*60, to = Stop_Time + interval*60)
# # AAA <- map2_dfr(.x = DateBegin_list, .y = DateEnd_list, Final_function_OPC)
# 
# list_date <- list(Begin = DateBegin_list, End = DateEnd_list)
# map_list  <- pmap(list_date, Final_function_OPC)
# 
# Lapply <- lapply(list_date, Final_function_OPC)
# 
# 
# map2(as.list(DateBegin_list), as.list(DateEnd_list), Final_function_OPC)

###############################################################################
###############################################################################


# make a loop for all the plots
 
# get Temp, RH, Vol and Samping
df_OPC_met <- df_OPC[ ,c("Stop_Time", "OPCTemp", "OPCHum", "OPCVol","OPCTsam")]
df_OPC_bins <-  df_OPC[, !names(df_OPC) %in% c("Start_Time", "Stop_Time", "density", "OPCTemp", "OPCHum", "OPCVol","OPCTsam")]
df_OPC_density <- df_OPC[, c("Start_Time", "density", "OPCTemp", "OPCHum")]


# attach one bin per per loop
for (i in 1:ncol(df_OPC)+1) {
pairs( cbind(df_OPC_density, df_OPC_bins[i]),
      lower.panel = panel.smooth, 
      upper.panel = panel.cor,nrow(),
      diag.panel  = panel.hist, 
     # labels = Labels, 
      main = paste0("Corr Density"), 
      cex.labels = 2) # cex.cor = 1.3
    print(i)
}




df_OPC[, c("Start_Time", "density", "OPCTemp", "OPCHum")]
########################################################################################################################
###### OLD STUFF #######################################################################################################
########################################################################################################################

########################################################################################################################
########################################################################################################################
##### Building the Time Series #########################################################################################
########################################################################################################################
########################################################################################################################

Dist_Ref.d_TS <- Distribution_Ref_TS(RefData, Date_Begin, Date_End, Min_APS = 0.62, Max_APS =  0.800) #

# Diameters overlapping of DMPS and APS
Diam_APS_to_keep  <- Dist_Ref.d_TS$Diam_APS
Diam_DMPS_to_keep <- Dist_Ref.d_TS$Diam_DMPS[(length(Dist_Ref.d_TS$Diam_DMPS)-3):length(Dist_Ref.d_TS$Diam_DMPS)]


DataXY.DMPS <- Dist_Ref.d_TS$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep)

DataXY.DMPS <- data.frame(x = log10(DataXY.DMPS$diameters),
                          y = log10(DataXY.DMPS$counts))
Model <- lm(y ~ x, data = DataXY.DMPS, model = TRUE, x = TRUE, y = TRUE)
Model.i <- list(Tidy = tidy(Model), 
                Augment = augment(Model), 
                Glance = glance(Model), 
                Call = Model$call, Coef = coef(Model),
                Model = Model)

DataXY.APS <- Dist_Ref.d_TS$counts_Ref %>% filter(diameters %in% Diam_APS_to_keep)
DataXY.APS <- data.frame(x = log10(DataXY.APS$diameters),
                         y = log10(DataXY.APS$counts))

# find optima value for the "density" of the particles
density_fun <- optimise( f_Error, c(0.9, 2) ) # limit ranges for the density
density_TS  <- density_fun$minimum

###########################################################################
# Fitting distribution to reference in order to estimate the model #######
###########################################################################

Dist_Ref.d_TS <- Distribution_Ref_TS(RefData, Date_Begin, Date_End, Min_APS = 0.673, Max_APS =  1.20) 

### we have the correct density in Model.i and the correct diameters of APS in Dist_ref
# return a list with "RefData_filtered" "counts_Ref"   "Diam_DMPS"         "Diam_APS" with the full diameters
Dist_Ref.full <- Distribution_Ref_TS(RefData, Date_Begin, Date_End, Min_APS = 0.777, Max_APS =  14.86000) #

## correction of APS diamters with density in Dist_Ref.d determined in Model.i  <- Density_Ref_log
Dist_Ref.full$counts_Ref$diameters <- Diam_APS_Corr(Dist_Ref.full$Diam_APS, Dist_Ref.full$counts_Ref$diameters, density_fun)

# Plot of distribution of reference counts versus diameter in log/log with fitted density (using corrected diameters)
Model.i.Gam  <- Density_Ref_log(Dist_Ref.full$counts_Ref, Density = NULL, Mod_type = 'GAM_GAUSS' )


#######################
### load OPCN3 ########
#######################

# Change names of OPC Sensor with diameters and return Refdata with names of Bin, diamters and counts
Dist_OPC_Sensor <- Distribution_OPC_Sensor(General.df = General.df, Date_Begin = Date_Begin, Date_End = Date_End, 
                                           bins_diameters = bins_diameters, Dist_Ref.full = Dist_Ref.full) 

# return predicted values of the OPC
Predict_Dist_OPC_Sensor <- Density_OPC_predict(Dist_OPC_Sensor$counts_OPC, Mod_type = 'GAM_GAUSS', density = density_fun$minimum)



# find Bin assignation for each diameter according to the Max APS diameter
Bins <- as.data.frame(t(bins_diameters))
Bins$Bins <- as.character(rownames(Bins))
rownames(Bins) <- NULL
names(Bins) <- c("diameters", "Bins")
Bins <- Bins %>%
    filter(diameters < max(Dist_Ref.full$Diam_APS, na.rm = T))



df_OPC_pred <- as.data.frame(t (round(Predict_Dist_OPC_Sensor$predict_counts, digits = 4)))
names(df_OPC_pred) <- paste0(c(Bins$Bins), "_pred")


df_OPC <- as.data.frame(t( round(Dist_OPC_Sensor$counts_OPC$counts, digits = 4)))
names(df_OPC) <- paste(c(Bins$Bins))


# data frame (TS) ###
df_OPC_Ref <- cbind(Date_Begin, Date_End, density_TS, df_OPC_pred, df_OPC)
df_OPC_Ref





########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################



# Remove Bin.APS and Bin.DMPS from column names RefData with diameters and 
# return a list with "RefData_filtered" "counts_Ref"   "Diam_DMPS"         "Diam_APS" ([1] 0.626 0.673 0.723 0.777 0.835 0.898 0.965 1.037 1.114 1.197)
Dist_Ref.d <- Distribution_Ref(RefData, Date_Begin, Date_End, Min_APS = 0.62, Max_APS =  0.800) #

# Diameters overlapping of DMPS and APS
Diam_APS_to_keep  <- Dist_Ref.d$Diam_APS
Diam_DMPS_to_keep <- Dist_Ref.d$Diam_DMPS[(length(Dist_Ref.d$Diam_DMPS)-3):length(Dist_Ref.d$Diam_DMPS)]

# plot log of DMPS
plot <-  ggplot() + 
    theme_bw() +
    geom_point(data = Dist_Ref.d$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep), 
               aes(log10(diameters), log10(counts)), stat = "identity", fill = "gray") +
    #   geom_line(data = augmented, aes(exp(x), exp(.fitted),  col = "Modelled"), size = 2) +
    theme(axis.title.x = element_text(colour  = "black", size = 15),
          axis.text.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 15, colour = "black")) +
    theme(axis.title.y = element_text(colour = "black", size = 15),
          axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 15, colour = "black")) +
    xlab(expression(paste("log10 of diameter (µm)"))) + 
    ylab(expression(paste("log10 of counts / log of diameters"))) 

plot

DataXY.DMPS <- Dist_Ref.d$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep)
DataXY.DMPS <- data.frame(x = log10(DataXY.DMPS$diameters),
                          y = log10(DataXY.DMPS$counts))
Model <- lm(y ~ x, data = DataXY.DMPS, model = TRUE, x = TRUE, y = TRUE)
Model.i <- list(Tidy = tidy(Model), 
                Augment = augment(Model), 
                Glance = glance(Model), 
                Call = Model$call, Coef = coef(Model),
                Model = Model)
Plot.DMPS <- Plot_Dist_Ref_log(Dist_Ref.d$counts_Ref %>% filter(diameters %in% Diam_DMPS_to_keep), Model.i = Model.i)
Plot.DMPS
                               
DataXY.APS <- Dist_Ref.d$counts_Ref %>% filter(diameters %in% Diam_APS_to_keep)
DataXY.APS <- data.frame(x = log10(DataXY.APS$diameters),
                         y = log10(DataXY.APS$counts))
Plot.DMPS.APS <- Plot.DMPS + geom_point(data = DataXY.APS, 
                                    aes(x, y), stat = "identity", col = "red") 
                                      
Plot.DMPS.APS

# find optima value for the "density" of the particles
density <- optimise( f_Error, c(0.9, 2) )   # limit ranges for the density

# correct diameters of the APS using the fitted density
DataXY.APS.d <- DataXY.APS %>% mutate(x = log10(10^x/sqrt(density$minimum)))
Plot.DMPS.APS.d <- Plot.DMPS.APS + geom_point(data = DataXY.APS.d, 
                                                aes(x, y), stat = "identity", shape = 17, col = "red", size = 3)
Plot.DMPS.APS.d

DataXY.APS$Predicted_APS <- predict(Model.i$Model, DataXY.APS.d)
Error <- sum(DataXY.APS$y - DataXY.APS$Predicted_APS)^2


#######################################
# Fitting distribution to reference ###
#######################################

### input parameters to be inserted by the USER
# Start and end date for selection of reference data
Date_Begin <- as.POSIXct("2018-09-30 01:00", tz = "UTC")
Date_End   <- as.POSIXct("2018-10-07 01:00", tz = "UTC") #"2018-10-07"  #"2018-09-30 01:00"

### Date for which we have OPC sensor data
#Mod_type <- "Normal"
#Mod_type <- "Normal_density"
#Mod_type <- "GAM_GAUSS"

# Remove Bin.APS and Bin.DMPS from column names RefData with diameters and 
# return a list with "RefData_filtered" "counts_Ref"   "Diam_DMPS"         "Diam_APS" ([1] 0.626 0.673 0.723 0.777 0.835 0.898 0.965 1.037 1.114 1.197)
Dist_Ref.d <- Distribution_Ref(RefData, Date_Begin, Date_End, Min_APS = 0.673, Max_APS =  1.20) #

# Plot of distribution of reference counts versus diameter in log/log
# lognormal distribution
Plot_Dist_Ref_log(Dist_Ref.d$counts_Ref)
# normal distribution
Plot_Dist_Ref(Dist_Ref.d$counts_Ref)

## Determine density that fits DPMS and ApS
# Model.i.d  <- Density_Ref_log(Dist_Ref.d$counts_Ref, Density = 1.5, Mod_type = 'Normal_density', Dist_Ref.d$Diam_APS )

## correction of APS diamters with density in Dist_Ref.d determined in Model.i  <- Density_Ref_log
# Dist_Ref.d$counts_Ref$diameters <- Diam_APS_Corr(Dist_Ref.d$Diam_APS, Dist_Ref.d$counts_Ref$diameters, Model.i.d )
Dist_Ref.d$counts_Ref$diameters <- Diam_APS_Corr(Dist_Ref.d$Diam_APS, Dist_Ref.d$counts_Ref$diameters, density)


## Plot of distribution of reference counts versus diameter in log/log with fitted density
# Plot_Dist_Ref_log(Dist_Ref.d$counts_Ref, Model.i = Model.i.d)
Plot_Dist_Ref_log(Dist_Ref.d$counts_Ref)



### we have the correct density in Model.i and the correct diameters of APS in Dist_ref
# return a list with "RefData_filtered" "counts_Ref"   "Diam_DMPS"         "Diam_APS" with the full diameters
Dist_Ref.full <- Distribution_Ref(RefData, Date_Begin, Date_End, Min_APS = 0.777, Max_APS =  14.86000) #
Plot_Dist_Ref_log(Dist_Ref.full$counts_Ref)

## Fitting of model with known density
# Model.i.full  <- Density_Ref_log(Dist_Ref.full$counts_Ref, Density = NULL, Mod_type = 'Normal')
## use corrected diameters
# Dist_Ref.full$counts_Ref$diameters <- Diam_APS_Corr(Dist_Ref.full$Diam_APS, Dist_Ref.full$counts_Ref$diameters, Model.i.d )
# Plot_Dist_Ref_log(Dist_Ref.full$counts_Ref, Model.i = Model.i.full)


# Plot of distribution of reference counts versus diameter in log/log with fitted density (using corrected diameters)
Model.i.Gam  <- Density_Ref_log(Dist_Ref.full$counts_Ref, Density = NULL, Mod_type = 'GAM_GAUSS' )
Plot_Dist_Ref_log(Dist_Ref.full$counts_Ref, Model.i = Model.i.Gam)

################################################################################################
#################### OPCN3 #####################################################################
################################################################################################

# Change names of OPC Sensor with diameters and return Refdata
Dist_OPC_Sensor <- Distribution_OPC_Sensor(General.df, Date_Begin, Date_End, bins_diameters) 

# diameter correction with the density
Dist_OPC_Sensor$counts_OPC$diameters <- Dist_OPC_Sensor$counts_OPC$diameters/sqrt(density$minimum)

# Predict_Dist_OPC_Sensor <- Density_OPC_predict(Dist_OPC_Sensor$counts_OPCN3, Mod_type = 'GAM_GAUSS')
Predict_Dist_OPC_Sensor <- Density_OPC_predict(Dist_OPC_Sensor$counts_OPC, Mod_type = 'GAM_GAUSS')

# Plot of distribution of OPC Sensor counts versus diameter in log/log
Plot_Distribution_OPC_Sensor(Dist_OPC_Sensor$counts_OPC, Predict_Dist_OPC_Sensor$predict_OPC)


################################################################################################
################################################################################################
### analyze other covariates RH and Temperature ################################################
################################################################################################

# load Temperature and RH from OPC sensor

# remove raws with NA values (from the OPC)


General.df_met <- General.df %>%
    select(date, 
           OPCTemp,
           OPCHum,
           OPCVol,
           OPCTsam)

# make all Bins as.numeric
General.df[,grep("Bin",names(General.df))] <- sapply(names(General.df[,grep("Bin",names(General.df))]), function(x) { as.numeric(General.df[,x])})

# transpose data by bin and filter diameters as the MAX of the APS
bins_diameters <- gather(bins_diameters, "bins", "diameters")
bins_diameters <- bins_diameters %>%
    filter(diameters < max(Dist_Ref.full$Diam_APS, na.rm = T))

# select only data from OPC bins
bins_OPC <- General.df[ , names(General.df) %in% bins_diameters$bins]

bins_OPC <- bins_OPC %>%
    gather("bins", "counts")

# merge with diameters
bins_OPC <- bins_OPC %>%
    left_join(bins_diameters, by = c("bins"))
bins_OPC <- bins_OPC %>%
    dplyr::select(-bins) 

#==============================================#
# calculate predicted value from the Gam model
#==============================================#
Predict_bins_OPC <- Density_OPC_predict(bins_OPC, Mod_type = 'GAM_GAUSS')
Predict_bins_OPC <- as.data.frame(Predict_bins_OPC)
#### to be added to the raw data....


bins_OPC$diameter <- as.factor(bins_OPC$diameter)

bins_OPC <- bins_OPC %>%
    dplyr::group_by(diameter) %>%
    dplyr::mutate(ID = row_number()) %>%
    tidyr::spread(diameter, counts) %>%
    select(-ID)
    


# addepredicted values

General.df_met <- cbind(General.df_met,
                        bins_OPC)

# remove all NA and NaN value
General.df_met <- na.omit(General.df_met)






##################################################################
##################################################################

Date_Begin <- "2018-10-07"
Date_End <- "2018-09-30"
Bin <- "Bin0"

General.df_met <- General.df %>%
    select(date, 
        Temperature,
        Relative_humidity)
General.df_met <- cbind(General.df_met,
                        bins_OPCN3)

General.df_met <- General.df_met %>%
    filter(date <= Date_Begin & date >= Date_End )

plot <- ggplot(General.df_met, aes(y = get(Bin), x = Relative_humidity)) +
    theme_bw() +
    geom_point(size = 2, show.legend = FALSE) +
    ylab(expression(paste("counts"))) +
    geom_smooth(data = General.df_met, method="lm", formula = y ~ x) + 
    theme(strip.text.y = element_text(angle = 0)) +
    theme(axis.title.x = element_text(colour="black", size=15),
          axis.text.x  = element_text(angle=0, vjust=0.5, hjust = 0.5, size=15, colour = "black")) +
    theme(axis.title.y = element_text(colour="black", size=15),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=15, colour = "black")) +
    ggtitle(expression(paste(R^2, " Bin0", " vs Relative Humidity"))) +
    theme(plot.title = element_text(lineheight=.8, face="bold", size = 15, hjust = 0.5)) +
    
    stat_poly_eq(data = General.df_met, formula = y ~ x, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.y = 150,
                 # label.y = max(General.df_met$Bin, na.rm = T),
                 parse = TRUE) 

plot


max(General.df_met$Bin0, na.rm = T)




# gather data according to bins
General.df_met <- gather(General.df_met, key= "bins" , value="counts", -Temperature, - Relative_humidity, na.rm = T)


# merge with diameters
General.df_met <- General.df_met %>%
    left_join(bins_diameters, by = c("bins"))
# remove bins
General.df_met <- General.df_met %>%
    select(-bins) %>%
    mutate(diameters = as.factor(diameters))



# colour_vector <- c("red", "blue", "black", "green", "cornflowerblue", "chocolate4", "darkblue",
#                    "darkgoldenrod3", "darkorange", "darkolivegreen4",
#                    "darkmagenta", "darkgreen", "darkcyan", "gray", "blue", "black", "green", 
#                    "cornflowerblue", "chocolate4", "darkblue", "darkgoldenrod3", "darkorange", "darkolivegreen4", 
#                    "goldenrod4", "darkred")
# colour_vector <- colour_vector[1:length( unique(General.df_met$diameters))]
# 
# 
# 
# # super plot of Counts by DIAMETERS, TEMPERATURE
# plot <- ggplot(General.df_met, aes(y = counts , x = Temperature)) +
#     theme_bw() +
#     geom_point( aes(colour  = diameters), size = 2) +
#     scale_color_manual(values = colour_vector)  # guide="none"
#  plot
# 
# # save plot
#  png("L:/ERLAP/Diffusion/Box Sync/AirSensEUR/Fieldtests/Shiny/JRC_11/General_data/Temp_BINS.jpg",
#      width = 1800, height = 1050, units = "px", pointsize = 30,
#      bg = "white", res = 150)
#  print(plot)
#  dev.off()
#  
# 
#  # super plot of Counts by DIAMETERS, RH
#  plot <- ggplot(General.df_met, aes(y = counts , x = Relative_humidity)) +
#      theme_bw() +
#      geom_point( aes(colour  = diameters), size = 2) +
#      scale_color_manual(values = colour_vector)  # guide="none"
#  plot
#  
#  png("L:/ERLAP/Diffusion/Box Sync/AirSensEUR/Fieldtests/Shiny/JRC_11/General_data/RH_BINS.jpg",
#      width = 1800, height = 1050, units = "px", pointsize = 30,
#      bg = "white", res = 150)
#  print(plot)
#  dev.off()

 
 # correlation by Bin with all covariates

 library("PerformanceAnalytics")
 AAA <- cor(General.df_met, use = "pairwise.complete.obs")
 chart.Correlation(AAA, histogram=F, pch=19, line = F)

################################################################################################
################################################################################################
################################################################################################
############################## OLD STUFF #######################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################


#################################
# load APS data #################
#################################

APS <- read.csv("APS data.csv", header = F)
names(APS) <- as.matrix(APS[1, ])

APS <- APS %>%
    select(-date)

# transpose data by bin
APS <- gather(APS, "diameters", "counts") 

# sum all counts in each bin
counts_APS <- APS %>%
    group_by(diameters) %>%
    summarise(counts = sum(counts, na.rm = T)) %>%
    mutate(diameters = as.numeric(diameters))


# combine OPCN3 and APS together

plot <-  ggplot() + 
    theme_bw() +
    geom_point(data = counts_APS, aes(diameters, counts), stat="identity", fill="gray") +
    geom_point(data = counts_OPCN3, aes(diameters, counts), stat="identity", fill="steelblue") +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    scale_fill_manual() +
    theme(axis.title.x = element_text(colour="black", size=15),
          axis.text.x  = element_text(angle=0, vjust=0.5, hjust = 0.5, size=15, colour = "black")) +
    theme(axis.title.y = element_text(colour="black", size=15),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=15, colour = "black")) +
    xlab(expression(paste("diameter (µm)"))) + 
    ylab(expression(paste("tot # of counts"))) 
plot


#################################
# load DMPS data ################
#################################

DMPS <- read.csv("DMPS data.csv", header = F)
names(DMPS) <- as.matrix(DMPS[1, ])

DMPS <- DMPS %>%
    select(-`date`)

# transpose data by bin
DMPS <- gather(DMPS, "diameters", "counts") 

# sum all counts in each bin
counts_DMPS <- DMPS %>%
    group_by(diameters) %>%
    summarise(counts = sum(counts, na.rm = T)) %>%
    mutate(diameters = as.numeric(diameters))

# ...also tranform nm into um
# counts_DMPS <- counts_DMPS %>%
#     mutate(diameters = diameters/1000)
    # filter(diameters >= min(bins_diameters$diameters))


# combine OPCN3, APS and DMPS together

plot <-  ggplot() + 
    theme_bw() +
    geom_point(data = counts_APS, aes(diameters, counts), stat="identity", fill="gray") +
    geom_point(data = counts_DMPS, aes(diameters, counts), stat="identity", fill="green") +
    geom_point(data = counts_OPCN3, aes(diameters, counts), stat="identity", fill="steelblue") +
 geom_smooth(data = counts_OPCN3, aes(diameters, counts), alpha = .3) +
    geom_smooth(data = counts_DMPS, aes(diameters, counts), alpha = .3) +
    geom_smooth(data = counts_APS, aes(diameters, counts),alpha = .3) +

     scale_x_continuous(trans='log10') +
     scale_y_continuous(trans='log10') +
    theme(axis.title.x = element_text(colour="black", size=15),
          axis.text.x  = element_text(angle=0, vjust=0.5, hjust = 0.5, size=15, colour = "black")) +
    theme(axis.title.y = element_text(colour="black", size=15),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=15, colour = "black")) +
    xlab(expression(paste("diameter (µm)"))) + 
    ylab(expression(paste("tot # of counts"))) 
    
plot

# with log values.....

plot <-  ggplot() + 
    theme_bw() +
    geom_point(data = counts_APS, aes(log(diameters), log(counts)), stat="identity", fill="gray") +
    geom_point(data = counts_DMPS, aes(log(diameters), log(counts)), stat="identity", fill="green") +
    geom_point(data = counts_OPCN3, aes(log(diameters), log(counts)), stat="identity", fill="steelblue") +
    geom_smooth(data = counts_OPCN3, aes(log(diameters), log(counts)), alpha = .3) +
    geom_smooth(data = counts_DMPS, aes(log(diameters), log(counts)), alpha = .3) +
    geom_smooth(data = counts_APS, aes(log(diameters), log(counts)),alpha = .3) +

    theme(axis.title.x = element_text(colour="black", size=15),
          axis.text.x  = element_text(angle=0, vjust=0.5, hjust = 0.5, size=15, colour = "black")) +
    theme(axis.title.y = element_text(colour="black", size=15),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=15, colour = "black")) +
    xlab(expression(paste("diameter (µm)"))) + 
    ylab(expression(paste("tot # of counts"))) 

plot





