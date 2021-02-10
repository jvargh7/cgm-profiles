

comparison_df <- paper3_df %>% 
  ungroup() %>% 
  mutate(cv_le36 = case_when(cv <= 36 ~ 1,
                             TRUE ~ 0)) %>%
  # # dplyr::filter(data_partition == "Training") %>%
  # dplyr::filter(lubridate::year(record_date) < 2019) %>%
  # Battelino 2019 Clinical Targets for Continuous Glucose Monitoring Data Interpretation: Recommendations From the International Consensus on Time in Range
  select(patient_id,file_name,
        average_sensor,gmi,cv,
        percent_time_over_250,
        percent_time_180_250,
        percent_time_70_180,
        percent_time_54_70,
        percent_time_under_54
        ) 




comparison_output <- comparison_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>% 
  nmf(.,2:8, nrun=100, seed=2019,.options='kv')
saveRDS(comparison_output,paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))

source(paste0(path_cgm_repo,"/package/nmf_amari.R"))

comparison_amari_output <- nmf_amari(comparison_output)

comparison_smallest_r <- comparison_amari_output$rabd_df %>% 
  group_by(r) %>% 
  summarize(d = mean(d)) %>% 
  ungroup() %>% 
  dplyr::filter(d == min(d)) %>% 
  dplyr::select(r) %>% 
  pull()


saveRDS(comparison_amari_output,paste0(path_cgm_working,"/paper3/comparison_amari_output.RDS"))

