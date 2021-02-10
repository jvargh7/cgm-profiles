
cgm_history_matching <- readRDS(paste0(path_cgm_working,"/cgm_history_matching.RDS"))

## Remove below for getting paper3_df from matched_df.RDS
# paper3_df <- readRDS(paste0(path_cgm_working,"/matched_df.RDS")) %>% 
paper3_df <- cgm_history_matching[["matched_df"]] %>% 
  ungroup() %>% 
  dplyr::filter(d_current_age >=18 & d_current_age <= 80) %>% 
  ## Remove below for getting paper3_df from matched_df.RDS
  # dplyr::select(mno, file_name,
  dplyr::select(patient_id, file_name,
                average_sensor:max_sensor,
                contains("percent"), contains("average_auc"),
                # -contains("_day"),
                # -contains("_night"),
                # -contains("_eve"),
                # -contains("_sleep"),
                n_peak:sd_pits_size,
                r_mage,
                
                -contains("_under_60"),
                -contains("_over_120"),
                -contains("_over_140"),
                -contains("_over_200")) %>%
  dplyr::mutate(
    percent_time_54_70 = percent_time_under_70 - percent_time_under_54,
    percent_time_180_250 = percent_time_over_180 - percent_time_over_250,
    
    percent_time_54_70_day = percent_time_under_70_day - percent_time_under_54_day,
    percent_time_180_250_day = percent_time_over_180_day - percent_time_over_250_day,
    
    percent_time_54_70_eve = percent_time_under_70_eve - percent_time_under_54_eve,
    percent_time_180_250_eve = percent_time_over_180_eve - percent_time_over_250_eve,
    
    percent_time_54_70_night = percent_time_under_70_night - percent_time_under_54_night,
    percent_time_180_250_night = percent_time_over_180_night - percent_time_over_250_night,
    
    percent_time_54_70_sleep = percent_time_under_70_sleep - percent_time_under_54_sleep,
    percent_time_180_250_sleep = percent_time_over_180_sleep - percent_time_over_250_sleep
    
    
  ) %>% 
  
  
  
  dplyr::filter(complete.cases(.)) %>% 
  
  left_join(cgm_history_matching[["matched_df"]] %>% 
              ungroup() %>% 
              dplyr::filter(d_current_age >=18 & d_current_age <= 80) %>% 
              ## Remove below for getting paper3_df from matched_df.RDS
              # group_by(mno) %>%
              group_by(patient_id) %>%
              arrange(desc(total_sensor_readings)) %>%
              mutate(max_readings = 1:n()) %>%
              ungroup() %>% 
              mutate(diabetes_status = case_when(d_current_age <= 19 ~ d_hba1c_classification_children,
                                                 d_current_age > 19 ~ d_hba1c_classification,
                                                 TRUE ~ NA_character_),
                     diabetes_status = case_when(diabetes_status %in% c("2, Prediabetic range",
                                                                        "3, Normal/Controlled range") ~
                                                   "3, Normal/Controlled range",
                                                 diabetes_status == "1, Diabetic range" ~ "1, Diabetic range",
                                                 TRUE ~ NA_character_),
                     diabcategory = case_when(diabcategory == "TYPE 1 DIABETES" ~ "T1D",
                                              diabcategory == "TYPE 2 DIABETES" ~ "T2D",
                                              diabcategory == "TYPE 2 DM - IRDM" ~ "T2D-IRDM",
                                              TRUE ~ NA_character_),
                     d_paper3_age = cut(d_current_age,breaks=c(18,35,50,80,Inf),  include.lowest=TRUE)) %>% 
              # dplyr::filter(max_readings == 1) %>% 
              dplyr::filter(!is.na(hba)) %>%
              dplyr::select(file_name,record_date,date_of_visit,
                            demo_patient_record,cgm_patient_record,
                            d_current_age,d_paper3_age,
                            diabcategory,diabetes_status,
                            gender, bmi,
                            d_bmi_category,
                            ageonset,d_duration,
                            insulin,
                            hba,bpsys,bpdia,
                            fbs,ppb,
                            gad,
                            cpf,
                            cps,
                            cho:vld,
                            # ure:mic,
                            
                            d_retinopathy:d_pvd_yn,
                            
                            # Add CGM variables without complete cases
                            conga_1,modd,lbgi,hbgi,j_index
                            
                            ),
            by= "file_name")


paper3_df <- paper3_df %>% 
  dplyr::filter(diabcategory=="T2D")  


saveRDS(paper3_df,paste0(path_cgm_working,"/paper3/paper3_df.RDS"))

## Exploratory:
paper3_df %>% 
  group_by(patient_id) %>% 
  tally() %>% 
  View()
