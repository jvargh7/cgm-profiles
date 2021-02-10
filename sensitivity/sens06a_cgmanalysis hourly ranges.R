source(paste0(path_cgm_repo,"/package/cgm_format.R"))


percent_time <- function(cgm_table,hours,interval,glycemic_ranges){
  
  
  
  hours_index <- base::which(base::as.numeric(base::format(cgm_table$timestamp, 
                                                               "%H")) %in% c(hours[1]:hours[2]))
  
  
  hours_sensor <- cgm_table$sensorglucose[hours_index]
  xaxis <- base::seq(from = 0, length.out = base::length(hours_sensor), 
                     by = (interval/60))
  
  
  xaxis[base::which(is.na(hours_sensor))] <- NA
  xaxis <- xaxis[!is.na(xaxis)]
  
  hours_sensor <- hours_sensor[!is.na(hours_sensor)]
  
  
  aucs <- pracma::cumtrapz(xaxis, hours_sensor)
  hours_auc <- aucs[base::length(hours_sensor)]
  
  
  hours_average <- base::mean(stats::na.omit(hours_sensor))
  hours_min <- base::min(hours_sensor)
  hours_max <- base::max(hours_sensor)
  hours_sd <- stats::sd(hours_sensor)
  
  summary_ranges <- map_dfr(glycemic_ranges,function(x){
    
    sr = tryCatch({
      BGinrange <- ifelse(hours_sensor %in% x, 1, 0);
      mins_spent_in_range <- base::sum(BGinrange, 
                                       na.rm = T) * (interval/60);
      
      percent_time_in_range <- (base::sum(BGinrange, 
                                          na.rm = T) * (interval/60))/(base::length(hours_sensor) * 
                                                                         (interval/60)) * 100;
      
      return(data.frame(range = paste0(min(x)," - ",max(x)),
                        mins = mins_spent_in_range,
                        percent = percent_time_in_range
      ))
    },
    error = function(e){
      return(data.frame(range = paste0(min(x)," - ",max(x)),
                        mins = NA,
                        percent = NA))
    });
    return(sr)
    
  }) %>% 
    dplyr::select(range,percent) %>% 
    pivot_wider(names_from = "range",values_from="percent",names_prefix = "GRpercent") 
  
  
  summary_ranges <- summary_ranges %>% 
    mutate(hours_range = paste0(c(hours[1],hours[2]),collapse="-"),
           hours_range_min = hours[1],
           hours_range_max = hours[2],
           hours_auc = hours_auc,
           hours_average = hours_average,
           hours_min = hours_min,
           hours_max = hours_max,
           hours_sd = hours_sd)
  
  return(summary_ranges)

}



## Derived from cgmanalysis:: cgmvariables function ----------------------
cgmanalysis_hourly <- function (inputdirectory, outputdirectory = tempdir(), 
                                      outputname = "REDCap Upload", log_file,
                                       printname = T) 
{
  cat("\n","Reading Files:",file=log_file,append=TRUE)
  
  #Full names set to FALSE
  files <- base::list.files(path = inputdirectory, full.names = FALSE)
  cat("\n","Total Files:",length(files),file=log_file,append=TRUE)
  
  dateparseorder <- c("mdy HM", "mdy HMS", "mdY HM", "mdY HMS", 
                      "dmy HM", "dmy HMS", "dmY HM", "dmY HMS", "Ymd HM", "Ymd HMS", 
                      "ymd HM", "ymd HMS", "Ydm HM", "Ydm HMS", "ydm HM", "ydm HMS")
  allhours <- 0:23
  
  hours_ranges <- lapply(allhours,function(x){
    if(x==23){
      out = c(x,0)
    } else{
      out = c(x,x+1)
    }
    return(out)
  })
  
  
  # y = files[10]
  
  cgm_hourly_ranges <- map_dfr(files,
                               
                               function(y){
                                 
                                 table <- utils::read.csv(paste0(inputdirectory,y), stringsAsFactors = FALSE, 
                                                          na.strings = c("NA", ""),sep = "\t");
                                 subjectid = str_extract(y,"^.*_") %>% str_replace(.,"_$","");
                                 
                                 summary_hour_ranges <- tryCatch({
                                   
                                   # Return data frame with each hour range in a row
                                   # Columns include percent in glycemic ranges, average, and other hourly statistics
                                   table$subjectid <- subjectid
                                   table <- cgm_format(table)
                                   
                                   
                                   table <- unique(table)
                                   if (printname == T) {
                                     cat("\n",basename(y),file=log_file,append=TRUE)
                                   }
                                   
                                   table$timestamp <- base::as.POSIXct(lubridate::parse_date_time(table$timestamp, 
                                                                                                  dateparseorder, tz = "Asia/Kolkata"))
                                   table$sensorglucose <- base::as.numeric(table$sensorglucose)
                                   
                                   interval <- pracma::Mode(base::diff(base::as.numeric(table$timestamp)))
                                   interval <- base::abs(interval)
                                   
                                   gr = list(c(70:180),c(0:54),c(54:69),c(181:250),c(250:1000))
                                   
                                   # One row for each of the 24 hours
                                   
                                   # hour_range_output per row -------------------
                                   hour_range_output <- map_dfr(hours_ranges,
                                                            function(x){
                                                              # One row per hour range
                                                              # Percent in glycemic ranges in columns
                                                              r_o <- percent_time(cgm_table=table,
                                                                                        hours=x,
                                                                                        interval=interval,
                                                                                        glycemic_ranges=gr);
                                                              return(r_o)
                                                              
                                                            }) %>% 
                                     mutate(file_name = y,
                                            subjectid = subjectid)
                                   
                                  
                                   
                                   return(hour_range_output)
                                   
                                   
                                 },
                                 error=function(e){
                                   error_h_r_o = map_dfr(hours_ranges,
                                                         function(x){
                                                           data.frame(hours_range = paste0(c(x[1],x[2]),collapse="-"),
                                                           hours_range_min = x[1],
                                                           hours_range_max = x[2],
                                                           file_name = y,
                                                           subjectid = subjectid)
                                                           
                                                         })
                                   
                                   return(error_h_r_o)
                                   
                                 });
                                 
                                 
                               })
  

    

  
  utils::write.csv(cgm_hourly_ranges, file = paste0(outputname,".csv"), row.names = FALSE)
  
  
  return(cgm_hourly_ranges)
}

# RUN FOR NEW INDICES ----------------------
inputdirectory = paste0(path_cgm_datasets,"/anonymized/")
outputdirectory = paste0(path_cgm_working,"/cgm_output/")
outputname = paste0("mdrf_cgm_summary_",Sys.Date())

log_file = paste0("mdrf/mdrf_cgmvariables.log")
cat("\n","Date: ",as.character(Sys.Date()),file=log_file,append=TRUE)

source(paste0(path_cgm_repo,"/package/cgmanalysis_cgmvariables.R"))
cgm_summary <- cgmanalysis_cgmvariables(inputdirectory,
                                        outputdirectory,
                                        outputname,log_file=log_file,
                                        magedef = "1sd",
                                        aboveexcursionlength = 15,
                                        belowexcursionlength = 15,run_ranges=FALSE)

source(paste0(path_cgm_repo,"/documentation/label_variable.R"))
cgm_summary <- label_variables(cgm_summary)
cgm_summary <- cgm_summary %>% 
  mutate_at(vars(percent_cgm_wear,num_days_good_data:hbgi),.f= function(x)as.numeric(levels(x)[x])) %>% 
  mutate(date_cgm_placement = as.character(date_cgm_placement)) %>% 
  mutate_if(is.factor,as.character)

saveRDS(cgm_summary,file=paste0(outputdirectory,
                                outputname,".RDS"))



# RUN FOR HOURLY ----------------------
outputname = paste0("mdrf_cgm_hourly_",Sys.Date())
log_file = paste0("mdrf/mdrf_cgm_hourly.log")
cat("\n","Date: ",as.character(Sys.Date()),file=log_file,append=TRUE)
cgm_hourly_summary <- cgmanalysis_hourly(inputdirectory,
                                        outputdirectory,
                                        outputname,log_file=log_file
                                        )

saveRDS(cgm_hourly_summary,paste0(path_cgm_working,"/cgm_output/",outputname,".RDS"))
