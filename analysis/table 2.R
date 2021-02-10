
nmf_output <- readRDS(paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))

tab2_df <- paper3_df %>% 
  mutate(
    d_retinopathy_status = case_when(npdr_Lt == 1 | npdr_rt == 1 ~ "1, NPDR",
                                     pdr_lt == 1 | pdr_rt == 1 ~ "2, PDR",
                                     npdr_Lt == 0 & npdr_rt == 0 & pdr_lt == 0 & pdr_rt == 0 ~ "0, NO",
                                     TRUE ~ NA_character_),
    d_nephropathy_status = case_when(d_mic < 30 ~ "0, NO",
                                     d_mic >=30 & d_mic < 300 ~ "1, MICRO",
                                     d_mic >= 300 ~ "2, MACRO",
                                     TRUE ~ NA_character_))

basis_r2 = basis(nmf_output$fit$`2`)
basis_r3 = basis(nmf_output$fit$`3`)
basis_r4 = basis(nmf_output$fit$`4`)


# tab2_df <- label_variables(tab2_df)

tab2_df$nmf_r2_hard = apply(basis_r2,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r3_hard = apply(basis_r3,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r4_hard = apply(basis_r4,1,function(x) which.max(x))%>% as.numeric()



basis_r2_v2 = basis(nmf_output2$fit$`2`)
basis_r3_v2 = basis(nmf_output2$fit$`3`)
basis_r4_v2 = basis(nmf_output2$fit$`4`)

tab2_df$nmf_r2_hard_v2 = apply(basis_r2_v2,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r3_hard_v2 = apply(basis_r3_v2,1,function(x) which.max(x))%>% as.numeric()
tab2_df$nmf_r4_hard_v2 = apply(basis_r4_v2,1,function(x) which.max(x))%>% as.numeric()

tab2_df <- tab2_df  %>% 
  mutate(cluster = factor(nmf_r3_hard,levels=c(3,1,2), 
                          labels=c("TIR profile","Hypo profile","Hyper profile")),
         cluster_v2 = factor(nmf_r3_hard_v2,levels=c(1,2,3), 
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
  ungroup() %>% 
  mutate_at(vars(d_retinopathy, d_retinopathy_yn,d_nephropathy, d_nephropathy_yn,d_pvd_yn,
                 pdr, npdr, npdr_or_pdr,macroprotein), 
            function(x) {
              y = as.character(x);
              y = case_when(is.na(y) ~ "Missing",
                            TRUE ~ y)
            }) %>% 
  mutate(mage_sd = r_mage/standard_deviation) 

tab2_df %>% 
  group_by(cluster) %>% 
  summarize_at(vars(gad),
               .f=function(x) {
                 paste0(
                   round(median(x,na.rm=TRUE),1),
                   "(",
                   round(quantile(x,probs=0.25,na.rm=TRUE),1),
                   ",",
                   round(quantile(x,probs=0.75,na.rm=TRUE),1),
                   ")"
                 )
               }) %>% View()


with(tab2_df,table(cluster,insulin))

table(tab2_df$cluster,tab2_df$cv>0.36)/nrow(tab1_df)
table(tab2_df$cluster,tab2_df$cv>0.36)


tab2_df %>% 
  group_by(cluster) %>% 
  summarize_at(.v =
                 vars(
                   percent_time_70_180,
                   percent_time_under_54,
                   percent_time_54_70,
                   percent_time_180_250,
                   percent_time_over_250,

                   
                 ),
               .f=function(x) {
                 paste0(
                   round(median(x),1),
                   " (",
                   round(quantile(x,probs=0.25),1),
                   ", ",
                   round(quantile(x,probs=0.75),1),
                   ")"
                 )
               }) %>% t() %>% write.csv(.,paste0(path_cgm_working,"/table2.csv"))







# NUMBER OF OBSERVATIONS ------------------
with(prevalence_df,table(d_retinopathy_status,d_duration>5,useNA="always"))
with(prevalence_df,table(d_nephropathy_status,d_duration>5,useNA="always"))

with(prevalence_df,table(d_retinopathy_yn,d_duration>5,useNA="always"))
with(prevalence_df,table(d_nephropathy_yn,d_duration>5,useNA="always"))


















# Results : Description of clusters (ynNA) --------------------
glm_retinopathy_ynNA = tab2_df %>% 
  mutate(d_retinopathy_ynNA = case_when(is.na(d_retinopathy_yn) ~ "2, Not tested",
                                 TRUE ~ d_retinopathy_yn),
         d_nephropathy_ynNA = case_when(is.na(d_nephropathy_yn) ~ "2, Not tested",
                                        TRUE ~ d_nephropathy_yn)) %>% 
  nnet::multinom(d_retinopathy_ynNA ~ cluster + ageonset + d_duration + gender,# + hba,
      family = binomial(),data=.)




  
glm_nephropathy_ynNA = tab2_df %>% 
  mutate(d_retinopathy_ynNA = case_when(is.na(d_retinopathy_yn) ~ "2, Not tested",
                                        TRUE ~ d_retinopathy_yn),
         d_nephropathy_ynNA = case_when(is.na(d_nephropathy_yn) ~ "2, Not tested",
                                        TRUE ~ d_nephropathy_yn)) %>% 
  nnet::multinom(d_nephropathy_ynNA ~ cluster + ageonset + d_duration + gender,# + hba,
                 family = binomial(),data=.)

summary(glm_retinopathy_ynNA)
paste0("OR : ",
       round(exp(coef(glm_retinopathy_ynNA)[1,2]),2),
       ", 95% CI:",
       round(exp(summary(glm_retinopathy_ynNA)$coefficients[1,2] - 
                   1.96*summary(glm_retinopathy_ynNA)$standard.errors[1,2] 
       ),2),
       ",",
       
       round(exp(summary(glm_retinopathy_ynNA)$coefficients[1,2] + 
                   1.96*summary(glm_retinopathy_ynNA)$standard.errors[1,2] 
                   ),2)
       
       )

paste0("OR : ",
       round(exp(coef(glm_nephropathy_ynNA)[1,2]),2),
       ", 95% CI:",
       round(exp(summary(glm_nephropathy_ynNA)$coefficients[1,2] - 
                   1.96*summary(glm_nephropathy_ynNA)$standard.errors[1,2] 
       ),2),
       ",",
       
       round(exp(summary(glm_nephropathy_ynNA)$coefficients[1,2] + 
                   1.96*summary(glm_nephropathy_ynNA)$standard.errors[1,2] 
       ),2)
       
)

paste0("OR : ",
       round(exp(coef(glm_retinopathy_ynNA)[1,3]),2),
       ", 95% CI:",
       round(exp(summary(glm_retinopathy_ynNA)$coefficients[1,3] - 
                   1.96*summary(glm_retinopathy_ynNA)$standard.errors[1,3] 
       ),2),
       ",",
       
       round(exp(summary(glm_retinopathy_ynNA)$coefficients[1,3] + 
                   1.96*summary(glm_retinopathy_ynNA)$standard.errors[1,3] 
       ),2)
       
)

paste0("OR : ",
       round(exp(coef(glm_nephropathy_ynNA)[1,3]),2),
       ", 95% CI:",
       round(exp(summary(glm_nephropathy_ynNA)$coefficients[1,3] - 
                   1.96*summary(glm_nephropathy_ynNA)$standard.errors[1,3] 
       ),2),
       ",",
       
       round(exp(summary(glm_nephropathy_ynNA)$coefficients[1,3] + 
                   1.96*summary(glm_nephropathy_ynNA)$standard.errors[1,3] 
       ),2)
       
)


# Results : Description of clusters (yn) --------------------






glm_nephropathy_ynNA2 = tab2_df %>% 
  mutate(d_retinopathy_ynNA = case_when(is.na(d_retinopathy_yn) ~ 0,
                                        d_retinopathy_yn == "0, No" ~ 0,
                                        TRUE ~ 1),
         d_nephropathy_ynNA = case_when(is.na(d_nephropathy_yn) ~ 0,
                                        d_nephropathy_yn == "0, No" ~ 0,
                                        TRUE ~ 1)) %>% 
  glm(d_nephropathy_ynNA ~ cluster + ageonset + d_duration + gender + hba,
                 family = binomial(),data=.)

broom::tidy(glm_retinopathy_ynNA2) %>% 
  mutate(display = paste0(round(exp(estimate),2), " (",
                          round(exp(estimate - 1.96*std.error),2),",",
                          round(exp(estimate + 1.96*std.error),2),")"
                          ))

broom::tidy(glm_nephropathy_ynNA2) %>% 
  mutate(display = paste0(round(exp(estimate),2), " (",
                          round(exp(estimate - 1.96*std.error),2),",",
                          round(exp(estimate + 1.96*std.error),2),")"
  ))


# Results - Prevalent Odds Sensitivity
glm_retinopathy_ynNA3 = tab2_df %>% 
  dplyr::filter(!is.na(d_retinopathy_yn)) %>% 
  mutate(d_retinopathy_yn = case_when(
                                        d_retinopathy_yn == "0, No" ~ 0,
                                        TRUE ~ 1)) %>% 
  glm(d_retinopathy_yn ~ cluster + ageonset + d_duration + gender + hba,
      family = binomial(),data=.)

glm_nephropathy_ynNA3 = tab2_df %>% 
  dplyr::filter(!is.na(d_nephropathy_yn)) %>% 
  mutate(d_nephropathy_yn = case_when(
    d_nephropathy_yn == "0, No" ~ 0,
    TRUE ~ 1)) %>% 
  glm(d_nephropathy_yn ~ cluster + ageonset + d_duration + gender + hba,
      family = binomial(),data=.)

broom::tidy(glm_retinopathy_ynNA3) %>% 
  mutate(display = paste0(round(exp(estimate),2), " (",
                          round(exp(estimate - 1.96*std.error),2),",",
                          round(exp(estimate + 1.96*std.error),2),")"
  ))

broom::tidy(glm_nephropathy_ynNA3) %>% 
  mutate(display = paste0(round(exp(estimate),2), " (",
                          round(exp(estimate - 1.96*std.error),2),",",
                          round(exp(estimate + 1.96*std.error),2),")"
  ))


# Risk of complications : surv_paper3.R -----
table(paper3_cluster_df$matched_retinopathy_yn,useNA="always")

with(cox_nephropathy_df,table(!is.na(yes), cluster))
with(cox_retinopathy_df,table(!is.na(yes), cluster))
