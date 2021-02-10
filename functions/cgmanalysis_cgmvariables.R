
source(paste0(path_cgm_repo,"/package/cgm_format.R"))

## Derived from cgmanalysis:: cgmvariables function
cgmanalysis_cgmvariables <- function (inputdirectory, outputdirectory = tempdir(), 
                                      outputname = "REDCap Upload", log_file,
                                      aboveexcursionlength = 35, belowexcursionlength = 10, magedef = "1sd", 
                                      congan = 1, 
                                      daystart = 6, dayend = 12, 
                                      evestart = 12, eveend = 18, 
                                      nightstart = 18, nightend = 23, 
                                      format = "rows", printname = T) 
{
  cat("\n","Reading Files:",file=log_file,append=TRUE)

  #Full names set to FALSE
  files <- base::list.files(path = inputdirectory, full.names = FALSE)
  cat("\n","Total Files:",length(files),file=log_file,append=TRUE)
  
  cgmupload <- base::as.data.frame(base::matrix(nrow = 0, ncol = base::length(files)))
  base::colnames(cgmupload) <- base::rep("Record", base::length(files))
  dateparseorder <- c("mdy HM", "mdy HMS", "mdY HM", "mdY HMS", 
                      "dmy HM", "dmy HMS", "dmY HM", "dmY HMS", "Ymd HM", "Ymd HMS", 
                      "ymd HM", "ymd HMS", "Ydm HM", "Ydm HMS", "ydm HM", "ydm HMS")
  allhours <- 0:23
  
  for (f in 1:base::length(files)) {
    # table <- utils::read.csv(files[f], stringsAsFactors = FALSE, 
    #                          na.strings = c("NA", ""))
    
    # Added sep = "\t" for data extraction
    # Added path to extract file name easily
    table <- utils::read.csv(paste0(inputdirectory,files[f]), stringsAsFactors = FALSE, 
                             na.strings = c("NA", ""),sep = "\t")
    
    
    # table$subjectid <- table$subjectid[1]
    # subjectid = table$subjectid[1]
    
    # Only patient ID : https://stackoverflow.com/questions/39086400/extracting-a-string-between-other-two-strings-in-r
    # subjectid <- str_match(files[f],"^(.*?)_")[,2]
    
    # If we need to incorporate second visit on same day
    subjectid = str_extract(files[f],"^.*_") %>% str_replace(.,"_$","")
    
    cgmupload["file_name",f] <- files[f]
    # cgmupload["subject_id", f] <- table$subjectid[1]
    cgmupload["subject_id", f] <- subjectid
    
    
      
    
    tryCatch({
    
    
      table$subjectid <- subjectid
      table <- cgm_format(table)
      
      
      table <- unique(table)
      if (printname == T) {
        cat("\n",basename(files[f]),file=log_file,append=TRUE)
      }
      
    table$timestamp <- base::as.POSIXct(lubridate::parse_date_time(table$timestamp, 
                                                                   dateparseorder, tz = "Asia/Kolkata"))
    table$sensorglucose <- base::as.numeric(table$sensorglucose)
    
    
    
    
    interval <- pracma::Mode(base::diff(base::as.numeric(table$timestamp)))
    interval <- base::abs(interval)
    cgmupload["date_cgm_placement", f] <- base::as.character(min(table$timestamp, 
                                                                 na.rm = T))
    totaltime <- base::as.numeric(base::difftime(base::max(table$timestamp, 
                                                           na.rm = T), base::min(table$timestamp, na.rm = T), 
                                                 units = "secs"))
    cgmupload["percent_cgm_wear", f] <- base::round(((base::length(which(!is.na(table$sensorglucose)))/(totaltime/interval)) * 
                                                       100), 2)
    
    x_rle <- rle(is.na(table$sensorglucose))
    x_rle$values <- ifelse(x_rle$values,"T","F")
    cgmupload["missingness_seq",f] <- paste0(x_rle$lengths,x_rle$values,collapse="_")
    
    table <- table[, -c(1)]
    table <- table[stats::complete.cases(table), ]
    cgmupload["num_days_good_data", f] <- base::round(base::length(table$sensorglucose)/(86400/interval))
    cgmupload["total_sensor_readings", f] <- base::as.numeric(base::length(base::which(!is.na(table$sensorglucose))))
    cgmupload["average_sensor", f] <- base::mean(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                                 na.rm = T)
    cgmupload["estimated_a1c", f] <- base::round((46.7 + 
                                                    (base::mean(table$sensorglucose[base::which(!is.na(table$sensorglucose))])))/28.7, 
                                                 digits = 1)
    cgmupload["gmi", f] <- base::round(3.31 + (0.02392 * 
                                                 base::mean(table$sensorglucose[base::which(!is.na(table$sensorglucose))])), 
                                       digits = 1)
    cgmupload["q1_sensor", f] <- base::as.numeric(base::summary(table$sensorglucose[base::which(!is.na(table$sensorglucose))])[2])
    cgmupload["median_sensor", f] <- base::as.numeric(base::summary(table$sensorglucose[base::which(!is.na(table$sensorglucose))])[3])
    cgmupload["q3_sensor", f] <- base::as.numeric(base::summary(table$sensorglucose[base::which(!is.na(table$sensorglucose))])[5])
    cgmupload["standard_deviation", f] <- stats::sd(table$sensorglucose[base::which(!is.na(table$sensorglucose))])
    cgmupload["cv", f] <- (stats::sd(table$sensorglucose[base::which(!is.na(table$sensorglucose))]))/base::mean(table$sensorglucose[base::which(!is.na(table$sensorglucose))])
    cgmupload["min_sensor", f] <- base::min(table$sensorglucose[base::which(!is.na(table$sensorglucose))])
    cgmupload["max_sensor", f] <- base::max(table$sensorglucose[base::which(!is.na(table$sensorglucose))])
    BGover120 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGover120[BGover120 <= 120] <- 0
    BGover120[BGover120 > 120] <- 1
    BG120.rle <- base::rle(BGover120)
    excursions120 <- base::as.numeric(BG120.rle$lengths[base::which(BG120.rle$values == 
                                                                      1)])
    cgmupload["excursions_over_120", f] <- base::length(base::which(excursions120 > 
                                                                      ((aboveexcursionlength * 60)/interval)))
    cgmupload["min_spent_over_120", f] <- base::sum(BGover120) * 
      (interval/60)
    cgmupload["percent_time_over_120", f] <- ((base::sum(BGover120) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGover140 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGover140[BGover140 <= 140] <- 0
    BGover140[BGover140 > 140] <- 1
    BG140.rle <- base::rle(BGover140)
    excursions140 <- base::as.numeric(BG140.rle$lengths[base::which(BG140.rle$values == 
                                                                      1)])
    cgmupload["excursions_over_140", f] <- base::length(base::which(excursions140 > 
                                                                      ((aboveexcursionlength * 60)/interval)))
    cgmupload["min_spent_over_140", f] <- base::sum(BGover140) * 
      (interval/60)
    cgmupload["percent_time_over_140", f] <- ((base::sum(BGover140) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGover180 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGover180[BGover180 <= 180] <- 0
    BGover180[BGover180 > 180] <- 1
    BG180.rle <- base::rle(BGover180)
    excursions180 <- base::as.numeric(BG180.rle$lengths[base::which(BG180.rle$values == 
                                                                      1)])
    cgmupload["excursions_over_180", f] <- base::length(base::which(excursions180 > 
                                                                      ((aboveexcursionlength * 60)/interval)))
    cgmupload["min_spent_over_180", f] <- base::sum(BGover180) * 
      (interval/60)
    cgmupload["percent_time_over_180", f] <- ((base::sum(BGover180) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGover200 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGover200[BGover200 <= 200] <- 0
    BGover200[BGover200 > 200] <- 1
    BG200.rle <- base::rle(BGover200)
    excursions200 <- base::as.numeric(BG200.rle$lengths[base::which(BG200.rle$values == 
                                                                      1)])
    cgmupload["excursions_over_200", f] <- base::length(base::which(excursions200 > 
                                                                      ((aboveexcursionlength * 60)/interval)))
    cgmupload["min_spent_over_200", f] <- base::sum(BGover200) * 
      (interval/60)
    cgmupload["percent_time_over_200", f] <- ((base::sum(BGover200) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    cgmupload["avg_excur_over_140_per_day", f] <- as.numeric(cgmupload["excursions_over_140", 
                                                                       f])/as.numeric(cgmupload["num_days_good_data", f])
    cgmupload["avg_excur_over_200_per_day", f] <- as.numeric(cgmupload["excursions_over_200", 
                                                                       f])/as.numeric(cgmupload["num_days_good_data", f])
    BGover250 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGover250[BGover250 <= 250] <- 0
    BGover250[BGover250 > 250] <- 1
    BG250.rle <- base::rle(BGover250)
    excursions250 <- base::as.numeric(BG250.rle$lengths[base::which(BG250.rle$values == 
                                                                      1)])
    cgmupload["excursions_over_250", f] <- base::length(base::which(excursions250 > 
                                                                      ((aboveexcursionlength * 60)/interval)))
    cgmupload["min_spent_over_250", f] <- base::sum(BGover250) * 
      (interval/60)
    cgmupload["percent_time_over_250", f] <- ((base::sum(BGover250) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGunder54 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGunder54[BGunder54 < 54] <- 1
    BGunder54[BGunder54 >= 54] <- 0
    BG54.rle <- base::rle(BGunder54)
    excursions54 <- base::as.numeric(BG54.rle$lengths[base::which(BG54.rle$values == 
                                                                    1)])
    cgmupload["excursions_under_54", f] <- base::length(base::which(excursions54 > 
                                                                      ((belowexcursionlength * 60)/interval)))
    cgmupload["min_spent_under_54", f] <- base::sum(BGunder54) * 
      (interval/60)
    cgmupload["percent_time_under_54", f] <- ((base::sum(BGunder54) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGunder60 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGunder60[BGunder60 < 60] <- 1
    BGunder60[BGunder60 >= 60] <- 0
    BG60.rle <- base::rle(BGunder60)
    excursions60 <- base::as.numeric(BG60.rle$lengths[base::which(BG60.rle$values == 
                                                                    1)])
    cgmupload["excursions_under_60", f] <- base::length(base::which(excursions60 > 
                                                                      ((belowexcursionlength * 60)/interval)))
    cgmupload["min_spent_under_60", f] <- base::sum(BGunder60) * 
      (interval/60)
    cgmupload["percent_time_under_60", f] <- ((base::sum(BGunder60) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGunder70 <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGunder70[BGunder70 < 70] <- 1
    BGunder70[BGunder70 >= 70] <- 0
    BG70.rle <- base::rle(BGunder70)
    excursions70 <- base::as.numeric(BG70.rle$lengths[base::which(BG70.rle$values == 
                                                                    1)])
    cgmupload["excursions_under_70", f] <- base::length(base::which(excursions70 > 
                                                                      ((belowexcursionlength * 60)/interval)))
    cgmupload["min_spent_under_70", f] <- base::sum(BGunder70) * 
      (interval/60)
    cgmupload["percent_time_under_70", f] <- ((base::sum(BGunder70) * 
                                                 (interval/60))/(base::length(table$sensorglucose) * 
                                                                   (interval/60))) * 100
    BGinrange <- base::as.numeric(table$sensorglucose[base::which(!is.na(table$sensorglucose))], 
                                  length = 1)
    BGinrange <- ifelse(BGinrange %in% 70:180, 1, 0)
    cgmupload["min_spent_70_180", f] <- base::sum(BGinrange) * 
      (interval/60)
    cgmupload["percent_time_70_180", f] <- ((base::sum(BGinrange) * 
                                               (interval/60))/(base::length(table$sensorglucose) * 
                                                                 (interval/60))) * 100
    # DAY - CHECK -----
    daytime_indexes <- base::which(base::as.numeric(base::format(table$timestamp, 
                                                                 "%H")) %in% c(daystart:(dayend-1)))
    if(length(daytime_indexes>0)){
      daytime_sensor <- table$sensorglucose[daytime_indexes]
      xaxis <- base::seq(from = 0, length.out = base::length(daytime_sensor), 
                         by = (interval/60))
      xaxis[base::which(is.na(daytime_sensor))] <- NA
      xaxis <- xaxis[!is.na(xaxis)]
      daytime_sensor <- daytime_sensor[!is.na(daytime_sensor)]
      aucs <- pracma::cumtrapz(xaxis, daytime_sensor)
      cgmupload["daytime_auc", f] <- aucs[base::length(daytime_sensor)]
      BGinrange <- ifelse(daytime_sensor %in% 70:180, 1, 0)
      cgmupload["min_spent_70_180_day", f] <- base::sum(BGinrange, 
                                                        na.rm = T) * (interval/60)
      cgmupload["percent_time_70_180_day", f] <- (base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                           (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor < 54, 1, 0)
      cgmupload["min_spent_under_54_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_54_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor < 60, 1, 0)
      cgmupload["min_spent_under_60_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_60_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor < 70, 1, 0)
      cgmupload["min_spent_under_70_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_70_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor > 180, 1, 0)
      cgmupload["min_spent_over_180_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_180_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor > 200, 1, 0)
      cgmupload["min_spent_over_200_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_200_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(daytime_sensor > 250, 1, 0)
      cgmupload["min_spent_over_250_day", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_250_day", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(daytime_sensor) * 
                                                                                             (interval/60)) * 100
      cgmupload["daytime_avg_sensor_glucose", f] <- base::mean(stats::na.omit(daytime_sensor))
      cgmupload["daytime_min_sensor_glucose", f] <- base::min(daytime_sensor)
      cgmupload["daytime_max_sensor_glucose", f] <- base::max(daytime_sensor)
      cgmupload["daytime_sd", f] <- stats::sd(daytime_sensor)
      
    }
    # EVE - CHECK -----
    evetime_indexes <- base::which(base::as.numeric(base::format(table$timestamp, 
                                                                 "%H")) %in% c(evestart:(eveend-1)))
    if(length(evetime_indexes>0)){
      evetime_sensor <- table$sensorglucose[evetime_indexes]
      xaxis <- base::seq(from = 0, length.out = base::length(evetime_sensor), 
                         by = (interval/60))
      xaxis[base::which(is.na(evetime_sensor))] <- NA
      xaxis <- xaxis[!is.na(xaxis)]
      evetime_sensor <- evetime_sensor[!is.na(evetime_sensor)]
      aucs <- pracma::cumtrapz(xaxis, evetime_sensor)
      cgmupload["evetime_auc", f] <- aucs[base::length(evetime_sensor)]
      BGinrange <- ifelse(evetime_sensor %in% 70:180, 1, 0)
      cgmupload["min_spent_70_180_eve", f] <- base::sum(BGinrange, 
                                                        na.rm = T) * (interval/60)
      cgmupload["percent_time_70_180_eve", f] <- (base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                           (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor < 54, 1, 0)
      cgmupload["min_spent_under_54_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_54_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor < 60, 1, 0)
      cgmupload["min_spent_under_60_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_60_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor < 70, 1, 0)
      cgmupload["min_spent_under_70_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_under_70_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor > 180, 1, 0)
      cgmupload["min_spent_over_180_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_180_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor > 200, 1, 0)
      cgmupload["min_spent_over_200_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_200_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(evetime_sensor > 250, 1, 0)
      cgmupload["min_spent_over_250_eve", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_over_250_eve", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(evetime_sensor) * 
                                                                                             (interval/60)) * 100
      cgmupload["evetime_avg_sensor_glucose", f] <- base::mean(stats::na.omit(evetime_sensor))
      cgmupload["evetime_min_sensor_glucose", f] <- base::min(evetime_sensor)
      cgmupload["evetime_max_sensor_glucose", f] <- base::max(evetime_sensor)
      cgmupload["evetime_sd", f] <- stats::sd(evetime_sensor)
      
    }
    
    
    
    
    # NIGHT- CHECK ------
    nighttime_indexes <- base::which(base::as.numeric(base::format(table$timestamp, 
                                                                 "%H")) %in% c(nightstart:(nightend-1)))
    if(length(nighttime_indexes>0)){
      nighttime_sensor <- table$sensorglucose[nighttime_indexes]
      xaxis <- base::seq(from = 0, length.out = base::length(nighttime_indexes), 
                         by = (interval/60))
      cgmupload["day_night_sensor_ratio", f] <- base::round(base::length(daytime_sensor)/base::length(nighttime_sensor), 
                                                            1)
      xaxis[base::which(is.na(nighttime_sensor))] <- NA
      xaxis <- xaxis[!is.na(xaxis)]
      nighttime_sensor <- nighttime_sensor[!is.na(nighttime_sensor)]
      aucs <- pracma::cumtrapz(xaxis, nighttime_sensor)
      cgmupload["nighttime_auc", f] <- aucs[base::length(nighttime_sensor)]
      BGinrange <- ifelse(nighttime_sensor %in% 70:180, 1, 
                          0)
      cgmupload["min_spent_70_180_night", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_70_180_night", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor < 54, 1, 0)
      cgmupload["min_spent_under_54_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_54_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor < 60, 1, 0)
      cgmupload["min_spent_under_60_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_60_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor < 70, 1, 0)
      cgmupload["min_spent_under_70_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_70_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor > 180, 1, 0)
      cgmupload["min_spent_over_180_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_180_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor > 200, 1, 0)
      cgmupload["min_spent_over_200_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_200_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(nighttime_sensor > 250, 1, 0)
      cgmupload["min_spent_over_250_night", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_250_night", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * 
                                                                                               (interval/60)) * 100
      cgmupload["nighttime_avg_sens_glucose", f] <- base::mean(stats::na.omit(nighttime_sensor))
      cgmupload["nighttime_min_sens_glucose", f] <- base::min(nighttime_sensor)
      cgmupload["nighttime_max_sens_glucose", f] <- base::max(nighttime_sensor)
      cgmupload["nighttime_sd", f] <- stats::sd(nighttime_sensor)
      
    }
    
    # SLEEP- CHECK ------
    sleeptime_indexes <- base::which(base::as.numeric(base::format(table$timestamp, 
                                                                   "%H")) %in% allhours[base::which(!(0:23 %in% c(daystart:(nightend-1))))])
    if(length(sleeptime_indexes>0)){
      sleeptime_sensor <- table$sensorglucose[sleeptime_indexes]
      xaxis <- base::seq(from = 0, length.out = base::length(sleeptime_indexes), 
                         by = (interval/60))
      cgmupload["day_sleep_sensor_ratio", f] <- base::round(base::length(daytime_sensor)/base::length(sleeptime_sensor), 
                                                            1)
      xaxis[base::which(is.na(sleeptime_sensor))] <- NA
      xaxis <- xaxis[!is.na(xaxis)]
      sleeptime_sensor <- sleeptime_sensor[!is.na(sleeptime_sensor)]
      aucs <- pracma::cumtrapz(xaxis, sleeptime_sensor)
      cgmupload["sleeptime_auc", f] <- aucs[base::length(sleeptime_sensor)]
      BGinrange <- ifelse(sleeptime_sensor %in% 70:180, 1, 
                          0)
      cgmupload["min_spent_70_180_sleep", f] <- base::sum(BGinrange, 
                                                          na.rm = T) * (interval/60)
      cgmupload["percent_time_70_180_sleep", f] <- (base::sum(BGinrange, 
                                                              na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                             (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor < 54, 1, 0)
      cgmupload["min_spent_under_54_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_54_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor < 60, 1, 0)
      cgmupload["min_spent_under_60_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_60_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor < 70, 1, 0)
      cgmupload["min_spent_under_70_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_under_70_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor > 180, 1, 0)
      cgmupload["min_spent_over_180_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_180_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor > 200, 1, 0)
      cgmupload["min_spent_over_200_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_200_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      BGinrange <- ifelse(sleeptime_sensor > 250, 1, 0)
      cgmupload["min_spent_over_250_sleep", f] <- base::sum(BGinrange, 
                                                            na.rm = T) * (interval/60)
      cgmupload["percent_time_over_250_sleep", f] <- (base::sum(BGinrange, 
                                                                na.rm = T) * (interval/60))/(base::length(sleeptime_sensor) * 
                                                                                               (interval/60)) * 100
      cgmupload["sleeptime_avg_sens_glucose", f] <- base::mean(stats::na.omit(sleeptime_sensor))
      cgmupload["sleeptime_min_sens_glucose", f] <- base::min(sleeptime_sensor)
      cgmupload["sleeptime_max_sens_glucose", f] <- base::max(sleeptime_sensor)
      cgmupload["sleeptime_sd", f] <- stats::sd(sleeptime_sensor)
      
    }
    
    
    
    
    
     sensorBG <- base::as.numeric(table$sensorglucose, length = 1)
    xaxis <- base::seq(from = 0, length.out = base::length(sensorBG), 
                       by = (interval/60))
    xaxis[base::which(is.na(sensorBG))] <- NA
    xaxis <- xaxis[!is.na(xaxis)]
    sensorBG <- sensorBG[!is.na(sensorBG)]
    aucs <- pracma::cumtrapz(xaxis, sensorBG)
    cgmupload["total_auc", f] <- aucs[base::length(sensorBG)]
    cgmupload["average_auc_per_day", f] <- base::as.numeric(cgmupload["total_auc", 
                                                                      f])/base::as.numeric(cgmupload["num_days_good_data", 
                                                                                                     f])
    
### AUC Computation -----
### AUC 180 ----
    sensorover180 <- table$sensorglucose
    sensorover180 <- sensorover180[sensorover180 > 180]
    sensorover180 <- sensorover180[!is.na(sensorover180)]
    xaxis <- base::seq(from = 0, length.out = base::length(sensorover180), 
                       by = (interval/60))
    if (base::length(sensorover180) > 1) {
      aucs <- pracma::cumtrapz(xaxis, sensorover180)
      aucs <- (aucs[base::length(sensorover180)]) - (xaxis[base::length(xaxis)] * 
                                                       180)
      cgmupload["auc_over_180", f] <- aucs
    }
    else {
      cgmupload["auc_over_180", f] <- 0
    }
    cgmupload["average_auc_180", f] <- base::as.numeric(cgmupload["auc_over_180", 
                                                                  f])/base::as.numeric(cgmupload["num_days_good_data", 
                                                                                                 f])
### AUC 250 ----
    sensorover250 <- table$sensorglucose
    sensorover250 <- sensorover250[sensorover250 > 250]
    sensorover250 <- sensorover250[!is.na(sensorover250)]
    xaxis <- base::seq(from = 0, length.out = base::length(sensorover250), 
                       by = (interval/60))
    if (base::length(sensorover250) > 1) {
      aucs <- pracma::cumtrapz(xaxis, sensorover250)
      aucs <- (aucs[base::length(sensorover250)]) - (xaxis[base::length(xaxis)] * 
                                                       250)
      cgmupload["auc_over_250", f] <- aucs
    }
    else {
      cgmupload["auc_over_250", f] <- 0
    }
    cgmupload["average_auc_250", f] <- base::as.numeric(cgmupload["auc_over_250", 
                                                                  f])/base::as.numeric(cgmupload["num_days_good_data", 
                                                                                                 f])
    
    
### AUC 54 ----
    sensorunder54 <- table$sensorglucose
    sensorunder54 <- sensorunder54[sensorunder54 < 54]
    sensorunder54 <- sensorunder54[!is.na(sensorunder54)]
    xaxis <- base::seq(from = 0, length.out = base::length(sensorunder54), 
                       by = (interval/60))
    if (base::length(sensorunder54) > 1) {
      aucs <- pracma::cumtrapz(xaxis, sensorunder54)
      # Reversed direction of computation
      aucs <- (xaxis[base::length(xaxis)] * 
                 54) - (aucs[base::length(sensorunder54)]) 
      cgmupload["auc_under_54", f] <- aucs
    }
    else {
      cgmupload["auc_under_54", f] <- 0
    }
    cgmupload["average_auc_54", f] <- base::as.numeric(cgmupload["auc_under_54", 
                                                                  f])/base::as.numeric(cgmupload["num_days_good_data", 
                                                                                                 f])
    
    
    
### AUC 70 ----
    sensorunder70 <- table$sensorglucose
    sensorunder70 <- sensorunder70[sensorunder70 < 70]
    sensorunder70 <- sensorunder70[!is.na(sensorunder70)]
    xaxis <- base::seq(from = 0, length.out = base::length(sensorunder70), 
                       by = (interval/60))
    if (base::length(sensorunder70) > 1) {
      aucs <- pracma::cumtrapz(xaxis, sensorunder70)
      # Reversed direction of computation
      aucs <- (xaxis[base::length(xaxis)] * 
                 70) - (aucs[base::length(sensorunder70)]) 
      cgmupload["auc_under_70", f] <- aucs
    }
    else {
      cgmupload["auc_under_70", f] <- 0
    }
    cgmupload["average_auc_70", f] <- base::as.numeric(cgmupload["auc_under_70", 
                                                                 f])/base::as.numeric(cgmupload["num_days_good_data", 
                                                                                                f])
    
    
    
    
    
# MAGE ----
    table$smoothed <- base::as.numeric(zoo::rollapply(zoo::zoo(table$sensorglucose), 
                                                      9, function(x) c(1, 2, 4, 8, 16, 8, 4, 2, 1) %*% 
                                                        (x/46), fill = NA))
    table$smoothed[1:4] <- base::mean(stats::na.omit(table$sensorglucose[1:4]))
    table$smoothed[(base::length(table$smoothed) - 3):base::length(table$smoothed)] <- base::mean(table$sensorglucose[(base::length(table$sensorglucose) - 
                                                                                                                         3):base::length(table$sensorglucose)])
    sd <- stats::sd(table$sensorglucose)
    tpoints <- pastecs::turnpoints(table$smoothed)
    peaks <- base::which(tpoints[["peaks"]] == TRUE)
    pits <- base::which(tpoints[["pits"]] == TRUE)
    
    
    if (tpoints[["firstispeak"]] == TRUE && base::length(peaks) != 
        base::length(pits)) {peaks <- peaks[2:base::length(peaks)]
    } else if (tpoints[["firstispeak"]] == FALSE && base::length(peaks) != 
             base::length(pits)) {pits <- pits[1:(base::length(pits) - 1)]}
    
    peak_times <- table$timestamp[peaks]
    pits_times <- table$timestamp[pits]
    
    cgmupload["n_peak", f] <- length(peaks)
    cgmupload["avg_smoothedpeak_size", f] <- mean(stats::na.omit(base::as.numeric(table$smoothed[peaks])))
    cgmupload["sd_smoothedpeak_size", f] <- sd(stats::na.omit(base::as.numeric(table$smoothed[peaks])))
    cgmupload["avg_peak_time_diff", f] <- mean(base::diff(stats::na.omit(base::as.numeric(peak_times))))/(3600)
    cgmupload["avg_peak_size", f] <- mean(stats::na.omit(base::as.numeric(table$sensorglucose[peaks])))
    cgmupload["sd_peak_size", f] <- sd(stats::na.omit(base::as.numeric(table$sensorglucose[peaks])))
    
    cgmupload["n_pits", f] <- length(pits)
    cgmupload["avg_smoothedpits_size", f] <- mean(stats::na.omit(base::as.numeric(table$smoothed[pits])))
    cgmupload["sd_smoothedpits_size", f] <- sd(stats::na.omit(base::as.numeric(table$smoothed[pits])))
    cgmupload["avg_pits_time_diff", f] <- mean(base::diff(stats::na.omit(base::as.numeric(pits_times))))/(3600)
    cgmupload["avg_pits_size", f] <- mean(stats::na.omit(base::as.numeric(table$sensorglucose[pits])))
    cgmupload["sd_pits_size", f] <- sd(stats::na.omit(base::as.numeric(table$sensorglucose[pits])))
    
    
    
    
    differences <- table$sensorglucose[peaks] - table$sensorglucose[pits]
    differences_abs <- abs(differences)
    
    time_differences <- stats::na.omit(base::as.numeric(peak_times) - base::as.numeric(pits_times))
    cgmupload["avg_time_diff",f] <- mean(time_differences[time_differences>0])/3600
    cgmupload["avg_time_diff_abs",f] <- mean(abs(time_differences))/3600
    
    if (magedef == "1sd") {
      cgmupload["r_mage", f] <- base::mean(stats::na.omit(differences[base::which(differences > 
                                                                                    sd)]))
      
      cgmupload["r_mage_abs", f] <- base::mean(stats::na.omit(differences_abs[base::which(differences_abs > 
                                                                                    sd)]))
      mage_time_differences <- stats::na.omit(time_differences[base::which(differences > 
                                                                   sd)])
      cgmupload["r_mage_time_diff",f] <- mean(abs(mage_time_differences))/(3600)
      
      mage_abs_time_differences <- stats::na.omit(time_differences[base::which(differences_abs > 
                                                                             sd)])
      cgmupload["r_mage_abs_time_diff",f] <- mean(abs(mage_abs_time_differences))/(3600)
    }
    else if (magedef == "1.5sd") {
      cgmupload["r_mage", f] <- base::mean(stats::na.omit(differences[base::which(differences > 
                                                                                    (sd * 1.5))]))
      
      cgmupload["r_mage_abs", f] <- base::mean(stats::na.omit(differences_abs[base::which(differences_abs > 
                                                                                    (sd * 1.5))]))
      
      mage_time_differences <- stats::na.omit(time_differences[base::which(differences > 
                                                                             (sd*1.5))])
      cgmupload["r_mage_time_diff",f] <- mean(mage_time_differences)/(3600)
      
      mage_abs_time_differences <- stats::na.omit(time_differences[base::which(differences_abs > 
                                                                                 (sd*1.5))])
      cgmupload["r_mage_abs_time_diff",f] <- mean(abs(mage_abs_time_differences))/(3600)
      
      
    }
    else if (magedef == "2sd") {
      cgmupload["r_mage", f] <- base::mean(stats::na.omit(differences[base::which(differences > 
                                                                                    (sd * 2))]))
      
      cgmupload["r_mage_abs", f] <- base::mean(stats::na.omit(differences_abs[base::which(differences_abs > 
                                                                                    (sd * 2))]))
      
      mage_time_differences <- stats::na.omit(time_differences[base::which(differences > 
                                                                             (sd*2))])
      cgmupload["r_mage_time_diff",f] <- mean(mage_time_differences)/(3600)
      
      mage_abs_time_differences <- stats::na.omit(time_differences[base::which(differences_abs > 
                                                                                 (sd*2))])
      cgmupload["r_mage_abs_time_diff",f] <- mean(abs(mage_abs_time_differences))/(3600)
      
    }
    else {
      cgmupload["r_mage", f] <- base::mean(stats::na.omit(differences[base::which(differences > 
                                                                                    magedef)]))
    }
    cgmupload["j_index", f] <- 0.001 * (base::mean(table$sensorglucose, 
                                                   na.rm = T) + stats::sd(table$sensorglucose, na.rm = T))^2
    n <- (congan * 3600)
    conga.times <- table$timestamp + n
    conga.times <- conga.times[!is.na(conga.times)]
    conga.times <- conga.times[base::order(conga.times)]
    conga.times <- conga.times[base::which(conga.times %in% 
                                             table$timestamp)]
    begin.times <- conga.times - n
    congas <- table$sensorglucose[base::which(table$timestamp %in% 
                                                conga.times)] - table$sensorglucose[base::which(table$timestamp %in% 
                                                                                                  begin.times)]
    cgmupload[base::paste0("conga_", congan), f] <- stats::sd(congas, 
                                                              na.rm = T)
    table$time <- lubridate::round_date(table$timestamp, 
                                        "5 minutes")
    table$time <- base::strftime(table$time, format = "%H:%M", 
                                 tz = "UTC")
    moddtable <- base::data.frame(base::matrix(ncol = 2, 
                                               nrow = base::length(unique(table$time))))
    base::colnames(moddtable) <- c("time", "mean_differences")
    moddtable$time <- base::unique(table$time)
    for (r in 1:nrow(moddtable)) {
      moddtable$mean_differences[r] <- base::mean(base::abs(base::diff(table$sensorglucose[base::which(table$time == 
                                                                                                         moddtable$time[r])])))
    }
    cgmupload["modd", f] <- base::mean(stats::na.omit(moddtable$mean_differences))
    a <- 1.084
    b <- 5.381
    y <- 1.509
    table$gluctransform <- y * ((base::log(table$sensorglucose)^a) - 
                                  b)
    table$rBG <- 10 * (table$gluctransform^2)
    rl <- table$rBG[base::which(table$gluctransform < 0)]
    rh <- table$rBG[base::which(table$gluctransform > 0)]
    cgmupload["lbgi", f] <- base::mean(stats::na.omit(rl))
    cgmupload["hbgi", f] <- base::mean(stats::na.omit(rh))
  },
  error=function(e){})
  }
  cgmupload <- base::cbind(`Variable / Field Name` = rownames(cgmupload), 
                           cgmupload)
  if (format == "rows") {
    cgmupload <- base::as.data.frame(base::t(cgmupload))
    cgmupload <- cgmupload[-1, ]
  }
  
  # Added code for creating outputdirectory if it doesn't exist
  if(!dir.exists(outputdirectory)){
    dir.create(outputdirectory)
  }
  
  filename <- base::paste0(outputdirectory, "/", outputname,sep = "")
  
  utils::write.csv(cgmupload, file = paste0(filename,".csv"), row.names = FALSE)
  
  
  return(cgmupload)
}




