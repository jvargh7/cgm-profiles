library(NMF)



tab1_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
  dplyr::filter(!is.na(hba)) %>% 
  group_by(patient_id) %>% 
  arrange(record_date) %>% 
  
  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>% 
  ungroup() 


logtab1_df <- tab1_df %>% 
  ungroup() %>% 
  select(patient_id,file_name,
         average_sensor,gmi,cv,
         percent_time_over_250,
         percent_time_180_250,
         percent_time_70_180,
         percent_time_54_70,
         percent_time_under_54) %>% 
  mutate(average_sensor = log(average_sensor))




logtab1_output <- logtab1_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>% 
  nmf(.,2:ncol(.), nrun=100, seed=2019,.options='kv')
saveRDS(logtab1_output,paste0(path_cgm_working,"/paper 3 review/logtab1_output.RDS"))


# Find hard clusters for train_df
basis_logtab1_r2 = basis(logtab1_output$fit$`2`)
basis_logtab1_r3 = basis(logtab1_output$fit$`3`)
basis_logtab1_r4 = basis(logtab1_output$fit$`4`)

logtab1_df$nmf_r2_log = apply(basis_logtab1_r2,1,function(x) which.max(x))%>% as.numeric()
logtab1_df$nmf_r3_log = apply(basis_logtab1_r3,1,function(x) which.max(x))%>% as.numeric()
logtab1_df$nmf_r4_log = apply(basis_logtab1_r4,1,function(x) which.max(x))%>% as.numeric()


nmf_output <- readRDS(paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))
basis_r2 = basis(nmf_output$fit$`2`)
basis_r3 = basis(nmf_output$fit$`3`)
basis_r4 = basis(nmf_output$fit$`4`)


tab2_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
  mutate(
    d_retinopathy_status = case_when(npdr_Lt == 1 | npdr_rt == 1 ~ "1, NPDR",
                                     pdr_lt == 1 | pdr_rt == 1 ~ "2, PDR",
                                     npdr_Lt == 0 & npdr_rt == 0 & pdr_lt == 0 & pdr_rt == 0 ~ "0, NO",
                                     TRUE ~ NA_character_),
    d_nephropathy_status = case_when(d_mic < 30 ~ "0, NO",
                                     d_mic >=30 & d_mic < 300 ~ "1, MICRO",
                                     d_mic >= 300 ~ "2, MACRO",
                                     TRUE ~ NA_character_)) 

tab2_df$nmf_r2_orig = apply(basis_r2,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r3_orig = apply(basis_r3,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r4_orig = apply(basis_r4,1,function(x) which.max(x))%>% as.numeric()

tab2_df <- tab2_df %>% 
  mutate(cluster = factor(nmf_r3_orig,levels=c(3,1,2), 
                          labels=c("TIR profile","Hypo profile","Hyper profile"))
  ) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(hba)) %>% 
  group_by(patient_id) %>% 
  arrange(record_date) %>% 
  
  mutate(data_partition = case_when(year(record_date) == 2019 ~ "Validation",
                                    year(record_date) < 2019 ~ "Training",
                                    TRUE ~ NA_character_)) %>% 
  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>%
  dplyr::select(-date_of_visit,-record_date) %>% 
  ungroup()


tab2_df$nmf_r2_log = apply(basis_logtab1_r2,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r3_log = apply(basis_logtab1_r3,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r4_log = apply(basis_logtab1_r4,1,function(x) which.max(x))%>% as.numeric()

tab2_df <- tab2_df %>% 
  mutate(cluster_log = factor(nmf_r3_log,levels=c(3,2,1), 
                             labels=c("TIR profile","Hypo profile","Hyper profile")))

caret::confusionMatrix(with(tab2_df,table(cluster,cluster_log)))
