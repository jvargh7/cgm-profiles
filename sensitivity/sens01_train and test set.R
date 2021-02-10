library(NMF)



tab1_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
  dplyr::filter(!is.na(hba)) %>% 
  group_by(patient_id) %>% 
  arrange(record_date) %>% 

  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>% 
  ungroup()


train_df = tab1_df %>% sample_frac(0.75)
test_df = anti_join(tab1_df,train_df,by="file_name")

train_df <- train_df %>% 
  ungroup() %>% 
  select(patient_id,file_name,
         average_sensor,gmi,cv,
         percent_time_over_250,
         percent_time_180_250,
         percent_time_70_180,
         percent_time_54_70,
         percent_time_under_54)




train_output <- train_df %>% 
  dplyr::select(-patient_id,-file_name) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>% 
  nmf(.,2:ncol(.), nrun=100, seed=2019,.options='kv')
saveRDS(train_output,paste0(path_cgm_working,"/paper 3 review/train_output.RDS"))


source(paste0(path_cgm_repo,"/package/nmf_amari.R"))

train_amari_output <- nmf_amari(train_output)

train_smallest_r <- train_amari_output$rabd_df %>% 
  group_by(r) %>% 
  summarize(d = mean(d)) %>% 
  ungroup() %>% 
  dplyr::filter(d == min(d)) %>% 
  dplyr::select(r) %>% 
  pull()


saveRDS(train_amari_output,paste0(path_cgm_working,"/paper 3 review/train_amari_output.RDS"))


# Find hard clusters for train_df
basis_train_r2 = basis(train_output$fit$`2`)
basis_train_r3 = basis(train_output$fit$`3`)
basis_train_r4 = basis(train_output$fit$`4`)

train_df$nmf_r2_hard = apply(basis_train_r2,1,function(x) which.max(x))%>% as.numeric()
train_df$nmf_r3_hard = apply(basis_train_r3,1,function(x) which.max(x))%>% as.numeric()
train_df$nmf_r4_hard = apply(basis_train_r4,1,function(x) which.max(x))%>% as.numeric()

# Apply coefficient matrix to test_df
test_mat = test_df %>% 
  ungroup() %>% 
  select(
         average_sensor,gmi,cv,
         percent_time_over_250,
         percent_time_180_250,
         percent_time_70_180,
         percent_time_54_70,
         percent_time_under_54) %>% 
  map_df(.x =.,.f = function(x) { 
    y = (x - min(x))/(max(x) - min(x)); 
    return(y)
  }) %>% 
  as.matrix(.)

coef_r2 = coef(train_output$fit$`2`)
coef_r3 = coef(train_output$fit$`3`)
coef_r4 = coef(train_output$fit$`4`)

basis_test_r2 = test_mat %*% MASS::ginv(coef_r2)
basis_test_r3 = test_mat %*% MASS::ginv(coef_r3)
basis_test_r4 = test_mat %*% MASS::ginv(coef_r4)


# Find hard clusters for test_df
test_df$nmf_r2_hard = apply(basis_test_r2,1,function(x) which.max(x))%>% as.numeric()
test_df$nmf_r3_hard = apply(basis_test_r3,1,function(x) which.max(x))%>% as.numeric()
test_df$nmf_r4_hard = apply(basis_test_r4,1,function(x) which.max(x))%>% as.numeric()



# Append train_df and test_df
traintest_df <- bind_rows(
  train_df %>% 
    mutate(set = "train"),
  test_df %>% 
    mutate(set = "test")
)

# Find extent of agreement with nmf_output
nmf_output <- readRDS(paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))

basis_r2 = basis(nmf_output$fit$`2`)
basis_r3 = basis(nmf_output$fit$`3`)
basis_r4 = basis(nmf_output$fit$`4`)


# tab2_df <- label_variables(tab2_df)


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
  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>%
  dplyr::select(-date_of_visit,-record_date) %>% 
  ungroup() %>% 
  mutate_at(vars(d_retinopathy, d_retinopathy_yn,d_nephropathy, d_nephropathy_yn,d_pvd_yn,
                 pdr, npdr, npdr_or_pdr,macroprotein), 
            function(x) {
              y = as.character(x);
              y = case_when(is.na(y) ~ "Missing",
                            TRUE ~ y)
            }) %>% 
  left_join(traintest_df %>% 
              dplyr::select(patient_id,file_name,nmf_r2_hard,nmf_r3_hard,nmf_r4_hard,set),
            by = c("patient_id","file_name")) %>% 
  mutate(cluster_tt = factor(nmf_r3_hard,levels=c(2,3,1), 
                             labels=c("TIR profile","Hypo profile","Hyper profile")))

with(tab2_df,table(nmf_r2_orig,nmf_r2_hard,set))
with(tab2_df,table(nmf_r3_orig,nmf_r3_hard,set))
with(tab2_df,table(nmf_r4_orig,nmf_r4_hard,set))

with(tab2_df[tab2_df$set=="test",],table(cluster,cluster_tt)) %>% 
  caret::confusionMatrix(.)

with(tab2_df[tab2_df$set=="train",],table(cluster,cluster_tt)) %>% 
  caret::confusionMatrix(.)
