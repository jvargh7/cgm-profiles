library(NMF)

# TAB 2 DF ---------------------
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
  left_join(readRDS(paste0(path_cgm_working,"/cgm_output/mdrf_cgm_summary_2021-02-07.rds")) %>% 
              dplyr::select(file_name,starts_with("grade")),
            by = "file_name") %>% 
  mutate_at(vars(starts_with("grade")),~as.numeric(.))



# DISTRIBUTIONS -----------------------

dens_plot <- function(df,var_name = "hba",title="A",xlab="HbA1c (%)"){
  cluster_colors = c("darkgreen","purple","red")
  names(cluster_colors) = c("TIR profile","Hypo profile","Hyper profile")
  plot_output <- ggplot(data=df,aes_string(x=var_name)) +
    geom_density(aes(fill=cluster),alpha=0.4) +
    theme_bw(base_size=12) +
    labs(x=xlab,y="Density",title = title) +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank())
  
  return(plot_output)
  
  
}

# SFIG 3 -----------------------
sfig3a <- dens_plot(tab2_df,"bmi","A","BMI (kg/m2)")
sfig3b <- dens_plot(tab2_df,"hba","B","HbA1c (%)")
sfig3c <- dens_plot(tab2_df,"fbs","C","Fasting glucose (mg/dL)")

sfig3 <- ggpubr::ggarrange(sfig3a,sfig3b,sfig3c,
                          ncol = 1,nrow=3,
                          legend="bottom",
                          common.legend=TRUE)

sfig3

# SFIG 4 -----------------------
sfig4a <- dens_plot(tab2_df,"percent_time_70_180","A","Time between 70-180 mg/dL (%)")
sfig4b <- dens_plot(tab2_df,"percent_time_under_70","B","Time below 70 mg/dL (%)")
sfig4c <- dens_plot(tab2_df,"percent_time_over_180","C","Time above 180 mg/dL")
sfig4d <- dens_plot(tab2_df,"cv","D","Coefficient of variation")
sfig4e <- dens_plot(tab2_df,"lbgi","E","Low Blood Glucose Index")
sfig4f <- dens_plot(tab2_df,"hbgi","F","High Blood Glucose Index")

sfig4 <- ggpubr::ggarrange(sfig4a,sfig4b,sfig4c,
                           sfig4d,sfig4e,sfig4f,
                           ncol = 3,nrow=2,
                           legend="bottom",
                           common.legend=TRUE)

sfig4


# SFIG 8 -----------------------
sfig8a <- dens_plot(tab2_df,"grade_hypo","A","%GRADE hypoglycemia")
sfig8b <- dens_plot(tab2_df,"grade_eu","B","%GRADE euglycemia")
sfig8c <- dens_plot(tab2_df,"grade_hyper","C","%GRADE hyperglycemia")

sfig8 <- ggpubr::ggarrange(sfig8a,sfig8b,sfig8c,
                           ncol = 1,nrow=3,
                           legend="bottom",
                           common.legend=TRUE)

sfig8










