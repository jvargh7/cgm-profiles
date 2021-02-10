# Function to approximately match CGM and Demographic datasets

closest_agp <- function(patient_id,
                        date_of_visit,
                        demo_df){
  
  patient_id = as.numeric(patient_id)
  date_of_visit = lubridate::ymd(date_of_visit)
  
  closest_agp_visit <- tryCatch({demo_df %>% 
      arrange(AGPdt) %>% 
      dplyr::filter(mno == patient_id,
                    AGPdt %in% c(date_of_visit+c(-40:40))) %>% 
      mutate(matching_visits = n(),
             day_difference = days(date_of_visit - AGPdt)@day,
             summary_matching = paste0("N: ",matching_visits,
                                       "; Min: ",min(day_difference),
                                       "; Abs min: ",min(abs(day_difference)),
                                       "; Abs max: ",max(abs(day_difference)),
                                       "; Max: ",max(day_difference))) %>% 
      arrange(abs(day_difference)) %>% 
      rename(closest_agpdt = AGPdt,
             closest_agp_number = d_agp_number) %>% 
    .[1,] },
    error = function(e){
      data.frame(mno = NA,
                 closest_agpdt = NA,
                 closest_agp_number = NA,
                 matching_visits = NA,
                 day_difference = NA,
                 summary_matching = NA)
    })
  
  
  closest_agp <- bind_cols(
    patient_id = patient_id,
    date_of_visit = date_of_visit,
    closest_agp_visit)
  
  return(closest_agp)
  
  
}
                        

approximate_matching <- function(
                                 unmatched_cgm=data.frame(),
                                 demo=data.frame(),
                                 path_demo = paste0(path_cgm_datasets,"/mdrf_demographics.RDS")) {
  
  if(nrow(demo)==0) {
    demo <-  readRDS(path_demo) %>% 
      mutate(mno_agpdt = paste0(mno,"_",AGPdt)) %>% 
      dplyr::select(mno,AGPdt,d_agp_number)
  }
  
  unmatched_cgm <- unmatched_cgm %>% 
    ungroup()
  
  approx_matched_cgm <- map2_df(.x=unmatched_cgm$patient_id,
                                .y=unmatched_cgm$date_of_visit,
            .f= function(x=.x,y=.y){closest_agp(x,y,demo)})
  
  # The data frame returned has the following columns:
  # - patient_id
  # - date_of_visit
  # - mno
  # - closest_agpdt
  # - closest_agp_number
  # - matching_visits: Number of matching visits
  # - day_difference: How far is the closest visit
  # - summary_matching: Summary of matching
  return(approx_matched_cgm)
  
  
  
}

## History Matching ------
# history_df <- mdrf_history
# param_df <- biochemistry
# history_id_var = "mno"
# history_date_var = "Record_Date"
# param_id_var = "patient_id"
# param_date_var = "Record_date"


# x = item_df$patient_id
# item_patient_id = x[1]
# y = item_df$record_date
# item_record_date = y[1]


closest_history <- function(item_patient_id,
                        item_record_date,
                        general_df,
                        interval_period=c(-13:13)){
  
  item_patient_id = as.numeric(item_patient_id)
  item_record_date = lubridate::ymd(item_record_date)
  
  closest_history_visit <- tryCatch({general_df %>% 
      arrange(record_date) %>% 
      dplyr::filter(patient_id == item_patient_id,
                    record_date %in% c(item_record_date+interval_period)) %>% 
      mutate(matching_visits = n(),
             day_difference = days(item_record_date - record_date)@day,
             summary_matching = paste0("N: ",matching_visits,
                                       "; Min: ",min(day_difference),
                                       "; Abs min: ",min(abs(day_difference)),
                                       "; Abs max: ",max(abs(day_difference)),
                                       "; Max: ",max(day_difference))) %>% 
      arrange(abs(day_difference)) %>% 
      rename(closest_visit = record_date) %>% 
      .[1,] },
      error = function(e){
        data.frame(mno = NA,
                   closest_visit = NA,
                   matching_visits = NA,
                   day_difference = NA,
                   summary_matching = NA)
      })
  
  
  closest_visit <- bind_cols(
    patient_id = item_patient_id,
    record_date = item_record_date,
    closest_history_visit) %>% 
    rename(patient_id_general = patient_id1)
  
  return(closest_visit)
  
  
}



history_matching2 <- function(
                              general_df = data.frame(),
                              item_df = data.frame(),
                              interval_period = c(-13:13)) {
  item_df <- item_df %>% 
    ungroup()
  
  approx_matched_item <- map2_df(.x=item_df$patient_id,
                                .y=item_df$record_date,
                                .f= function(x=.x,y=.y){closest_history(x,y,general_df,interval_period)})
  
  # The data frame returned has the following columns:
  # - patient_id
  # - record_date
  # - patient_id1
  # - closest_visit
  # - closest_visit_number
  # - matching_visits: Number of matching visits
  # - day_difference: How far is the closest visit
  # - summary_matching: Summary of matching
    
    
  approx_matched_item <- item_df %>% 
    
    mutate(patient_record = paste0(patient_id,"_",as.character(record_date))) %>% 
    
    left_join(approx_matched_item %>% 
                mutate(patient_record = paste0(patient_id,"_",as.character(record_date))) %>% 
                dplyr::select(-patient_id,-record_date),
              by="patient_record")
  
  return(approx_matched_item)
  
  
}



# history_matching = function(history_df = data.frame(),
#                             param_df = data.frame(),
#                             history_date_var=character(),
#                             history_id_var=character(),
#                             param_date_var=character(),
#                             param_id_var=character()) {
#   
#   history_df <- history_df  %>% 
#     dplyr::mutate(id = !!rlang::sym(history_id_var),
#                   date = !!rlang::sym(history_date_var)) %>% 
#     dplyr::mutate(id_date = paste0(id,
#                                    "_",
#                                    date))
#   
#   param_df <- param_df %>% 
#     dplyr::mutate(id = !!rlang::sym(param_id_var),
#                   date = !!rlang::sym(param_date_var)) %>% 
#     dplyr::mutate(id_date = paste0(id,
#                                    "_",
#                                    date))
#   
#   exact_match <- inner_join(history_df,
#                             param_df %>% dplyr::select(-param_id_var,-param_date_var,-id,-date),
#                             by="id_date")
#   
#   unmatched_df <- history_df %>% 
#     anti_join(exact_match,by="id_date")
#   
#   
#   approx_matched <- map2_df(.x=unmatched_df$id,
#                                 .y=unmatched_df$date,
#                                 .f= function(x=.x,y=.y){closest_history(x,y,param_df)})
#   
#   
#   
#   
# }