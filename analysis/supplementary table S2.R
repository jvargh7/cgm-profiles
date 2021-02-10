library(NMF)


nmf_df <- paper3_df %>% 
  ungroup() %>% 
  # dplyr::filter(data_partition == "Training") %>%
  # dplyr::filter(lubridate::year(record_date) < 2019) %>%
  select(patient_id,file_name,
         average_sensor:max_sensor,
         r_mage,
         contains("percent"), contains("average_auc"),
         # -contains("_day"),
         # -contains("_night"),
         # -contains("_eve"),
         # -contains("_sleep"),
         n_peak:sd_pits_size) %>% 
  dplyr::select(-contains("_under_70"),
                -contains("_over_180"))




nmf_output <- nmf_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>% 
  nmf(.,2:12, nrun=100, seed=2019,.options='kv')
saveRDS(nmf_output,paste0(path_cgm_working,"/paper3/nmf_output.RDS"))


source(paste0(path_cgm_repo,"/package/nmf_amari.R"))

amari_output <- nmf_amari(nmf_output)

smallest_r <- amari_output$rabd_df %>% 
  group_by(r) %>% 
  summarize(d = mean(d)) %>% 
  ungroup() %>% 
  dplyr::filter(d == min(d)) %>% 
  dplyr::select(r) %>% 
  pull()


saveRDS(amari_output,paste0(path_cgm_working,"/paper3/amari_output.RDS"))

# # Check:------------
# diss = matrix(0,nrow=10,ncol=10)
# for (a in 2:10){
#   for (b in 1:(a-1)){
#     
#     d1 = coef(nmf_output$fit$`3`[[a]]) %>% t()
#     d2 = coef(nmf_output$fit$`3`[[b]]) %>% t()
#     
#     C <- matrix(0,nrow=3,ncol=3)
#     
#     for(j in c(1:3)){  
#       for(k in c(1:3)) {
#         C[j,k] = cor(d1[,j],d2[,k])
#       }
#     }
#     
#     sum1 = 0
#     for(j in c(1:3)){
#       sum1 = 1/3*(1-max(C[j,],na.rm=TRUE)) + sum1
#     }
#     
#     sum2 = 0
#     for(k in c(1:3)){
#       sum2 = 1/3*(1-max(C[,k],na.rm=TRUE)) + sum2
#     }  
#     
#     diss[a,b] = (sum1 + sum2)/2
#   }
# }  




# Generating Random data NMF ------
# shuffle original data
V.random <- nmf_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>%
  randomize()
# estimate quality measures from the shuffled data (use default NMF algorithm)

random_output <- nmf( V.random , 2:12, nrun=100, seed=2019,.options='kv')
saveRDS(random_output,paste0(path_cgm_working,"/paper3/random_output.RDS"))

comparison_V.random <- comparison_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>%
  randomize()
# estimate quality measures from the shuffled data (use default NMF algorithm)
comparison_random_output <- nmf(comparison_V.random, rank = 2:8, nrun=100, seed=2019,.options='kv')
saveRDS(comparison_random_output,paste0(path_cgm_working,"/paper3/comparison_random_output.RDS"))
