library(NMF)



tab1_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
  dplyr::filter(!is.na(hba)) %>% 
  group_by(patient_id) %>% 
  arrange(record_date) %>% 
  
  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>% 
  ungroup() %>% 
  left_join(readRDS(paste0(path_cgm_working,"/cgm_output/mdrf_cgm_summary_2021-02-22.rds")) %>% 
              dplyr::select(file_name,starts_with("grade"),G),
            by = "file_name") %>% 
  mutate_at(vars(starts_with("grade"),G),~as.numeric(.)) %>% 
  select(patient_id,file_name,
         G,average_sensor,gmi,cv,
         percent_time_over_250,
         percent_time_180_250,
         percent_time_70_180,
         percent_time_54_70,
         percent_time_under_54) %>% 
  mutate(log_sensor = log(average_sensor))


tab1_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  cor(.) %>% 
  View()
