rm(list = ls())

library(parallel)
library(zoo)
library(rChoiceDialogs)

General.Rdata.File <- "D:\\Profil_Michel\\Desktop\\Bureau\\Diffusion\\Box Sync\\AirSensEUR\\Fieldtests\\Shiny\\NL_01min\\General_data\\General.Rdata"
if(file.exists(General.Rdata.File)) load(General.Rdata.File) else load(rChoiceDialogs::jchoose.files(default = getwd(),
                                                                                                     caption = "Select General Rdata files", 
                                                                                                     multi = FALSE,
                                                                                                     filters = "*.Rdata", 
                                                                                                     modal = canUseJavaModal()))
no_cores <- detectCores()-1
c1 <- makeCluster(no_cores)
c1 <- makeCluster(1)
#clusterExport(c1, "General.df")

clusterEvalQ(c1, {
     
     library(zoo)
     utmaxmin <- function(y, threshold = 1) {
         # y            : vector of numeric (NA possible) on whcih to calculate Median - threshold * MAD amd Median + threshold * MAD
         #                where MAD is the Median Absolute Deviation i. e. median of the absolute deviations from the median,
         # Threshold    : numeric, default 1scale factor.
         # Return a nimric vector c(Median - threshold * cMedian(abs(y - Median(y))), Median + threshold * cMedian(abs(y - Median(y))) )
          Median <- stats::median(y, na.rm = TRUE)
          MAD    <- stats::mad(y, constant = threshold, na.rm = TRUE)
          return(c(Median - MAD, Median + MAD))
     }
     
     My.rm.Outliers <- function( y,  General.df, window, threshold, 
                                 date, ymin = NULL, ymax = NULL, ThresholdMin = NULL) {
         ### Extracting index of the values outside [ymin, ymax] and the outliers tolerance with rollapply (zmax and zmin together)
         
         # y                     = name of column of General.df with time series on which ouliers are detected
         # General.df            = Data.frame with times series in columns including vector on which to detect outliers
         # window                = width of the window used to compute median and average in number of data
         # threshold             = coefficient that muliplied by the difference between median and average that is exceeded result in outlier
         # date                  = Character, name of columns if General.df giving x axis, date values in Posixct class
         # ymin                  = minimum values for y, for example to remove negative values
         # ymax                  = maximum values for y, for example to limit highest values
         # ThresholdMin          = minimum values for zmin, the minimum values that evidence outliers when exceeded
         # utmaxmin            
         
         # Return                = dataframe with date then logical values Low_values (<min), High_values(>ymax), 
         #                         OutliersMin(< Median - threshold * MAD, OutliersMax (> Median + threshold *MAD, 
         #                         and values zmin (Median - threshold * MAD) and zmax (Median + threshold * MAD)
         
         
         # Removing values lower than ymin
         if (!is.null(ymin)) Low_values  <- (General.df[,y] < ymin) else Low_values  <- rep(FALSE, nrow(General.df))
         if (!is.null(ymax)) High_values <- (General.df[,y] > ymax) else High_values <- rep(FALSE, nrow(General.df))
         
         # max and max limits
         zmax.zmin <- zoo::rollapply(zoo(General.df[,y]), width = window, FUN = utmaxmin, threshold = threshold, align = "center", partial = TRUE)
         
         # Changing negative values of the minimum of interval of tolerance by ThresholdMin
         if (!is.null(ThresholdMin)) {zmax.zmin[which(zmax.zmin[,2] < ThresholdMin),2] <- ThresholdMin}
         
         # OutliersMax and OutliersMin
         OutliersMin <- General.df[,y] < zmax.zmin[,1]
         OutliersMax <- General.df[,y] > zmax.zmin[,2]
         
         # data frame to return
         df <- data.frame(date        = General.df[,date], 
                          Low_values  = Low_values, 
                          High_values = High_values, 
                          OutliersMin = OutliersMin, 
                          OutliersMax = OutliersMax, 
                          zmin        = zmax.zmin[,1], 
                          zmax        = zmax.zmin[,2]) 
         return(df)
     }
     
     roll <- function(y, General.df) {
          zoo::rollapply(zoo(General.df[,y]), width = 19, FUN = utmaxmin, threshold = 14, align = "center", partial = TRUE)
     }
     
})

list.gas.sensors <- c("Carbon_monoxide","Nitric_oxide","Nitrogen_dioxide","Ozone")
tm1 <- system.time(
     L_rr <- parallel::parLapply(c1, list.gas.sensors, roll, General.df = General.df)
)
tm2 <- system.time(
     L_rr <- parallel::parLapply(c1,list.gas.sensors, My.rm.Outliers, General.df = General.df, window = 19 , threshold = 14, date = "date")
)
stopCluster(c1)
tm3 <- system.time(My.rm.Outliers( y = list.gas.sensors[1],  General.df = General.df, window = 19 , threshold = 14, date = "date"))
