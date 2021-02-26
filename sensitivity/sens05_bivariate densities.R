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
  left_join(readRDS(paste0(path_cgm_working,"/cgm_output/mdrf_cgm_summary_2021-02-22.rds")) %>% 
              dplyr::select(file_name,starts_with("grade")),
            by = "file_name") %>% 
  mutate_at(vars(starts_with("grade")),~as.numeric(.))



# DISTRIBUTIONS -----------------------

bivariate_plot <- function(df,xvar_name = "hba",yvar_name = "cv",title="A",xlab="HbA1c (%)",ylab="Coefficient of Variation"){
  cluster_colors = c("darkgreen","purple","red")
  names(cluster_colors) = c("TIR profile","Hypo profile","Hyper profile")
  plot_output <- ggplot(data=df,aes_string(x=xvar_name,y=yvar_name,fill="cluster",col="cluster")) +
    geom_point(size=0.4,alpha=0.2) +
    # stat_ellipse(geom="polygon",alpha=0.2,type="norm") +
    stat_density2d() +
    theme_bw(base_size=12) +
    labs(x=xlab,y=ylab,title = title) +
    scale_color_manual(values = cluster_colors) +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.title = element_blank())
  
  return(plot_output)
  
  
}


# SFIG 5 -----------------------
sfig5a <- bivariate_plot(tab2_df,"bmi","percent_time_70_180","A","Body Mass Index (kg/m2)","Time between 70-180 mg/dL (%)")
sfig5b <- bivariate_plot(tab2_df,"bmi","percent_time_under_70","B","Body Mass Index (kg/m2)","Time below 70 mg/dL (%)")
sfig5c <- bivariate_plot(tab2_df,"bmi","percent_time_over_180","C","Body Mass Index (kg/m2)","Time above 180 mg/dL")
sfig5d <- bivariate_plot(tab2_df,"bmi","cv","D","Body Mass Index (kg/m2)","Coefficient of variation")
sfig5e <- bivariate_plot(tab2_df,"bmi","lbgi","E","Body Mass Index (kg/m2)","Low Blood Glucose Index")
sfig5f <- bivariate_plot(tab2_df,"bmi","hbgi","F","Body Mass Index (kg/m2)","High Blood Glucose Index")

sfig5 <- ggpubr::ggarrange(sfig5a,sfig5b,sfig5c,
                           sfig5d,sfig5e,sfig5f,
                           ncol = 3,nrow=2,
                           legend="bottom",
                           common.legend=TRUE)

sfig5


# SFIG 6 -----------------------
sfig6a <- bivariate_plot(tab2_df,"hba","percent_time_70_180","A","HbA1c (%)","Time between 70-180 mg/dL (%)")
sfig6b <- bivariate_plot(tab2_df,"hba","percent_time_under_70","B","HbA1c (%)","Time below 70 mg/dL (%)")
sfig6c <- bivariate_plot(tab2_df,"hba","percent_time_over_180","C","HbA1c (%)","Time above 180 mg/dL")
sfig6d <- bivariate_plot(tab2_df,"hba","cv","D","HbA1c (%)","Coefficient of variation")
sfig6e <- bivariate_plot(tab2_df,"hba","lbgi","E","HbA1c (%)","Low Blood Glucose Index")
sfig6f <- bivariate_plot(tab2_df,"hba","hbgi","F","HbA1c (%)","High Blood Glucose Index")

sfig6 <- ggpubr::ggarrange(sfig6a,sfig6b,sfig6c,
                           sfig6d,sfig6e,sfig6f,
                           ncol = 3,nrow=2,
                           legend="bottom",
                           common.legend=TRUE)

sfig6

# SFIG 7 -----------------------
sfig7a <- bivariate_plot(tab2_df,"fbs","percent_time_70_180","A","Fasting glucose (mg/dL)","Time between 70-180 mg/dL (%)")
sfig7b <- bivariate_plot(tab2_df,"fbs","percent_time_under_70","B","Fasting glucose (mg/dL)","Time below 70 mg/dL (%)")
sfig7c <- bivariate_plot(tab2_df,"fbs","percent_time_over_180","C","Fasting glucose (mg/dL)","Time above 180 mg/dL")
sfig7d <- bivariate_plot(tab2_df,"fbs","cv","D","Fasting glucose (mg/dL)","Coefficient of variation")
sfig7e <- bivariate_plot(tab2_df,"fbs","lbgi","E","Fasting glucose (mg/dL)","Low Blood Glucose Index")
sfig7f <- bivariate_plot(tab2_df,"fbs","hbgi","F","Fasting glucose (mg/dL)","High Blood Glucose Index")

sfig7 <- ggpubr::ggarrange(sfig7a,sfig7b,sfig7c,
                           sfig7d,sfig7e,sfig7f,
                           ncol = 3,nrow=2,
                           legend="bottom",
                           common.legend=TRUE)

sfig7


# SFIG 7 -----------------------
sfig9a <- bivariate_plot(tab2_df,"percent_time_under_70","percent_time_over_180","A","Time below 70 mg/dL (%)","Time above 180 mg/dL")
sfig9b <- bivariate_plot(tab2_df,"percent_time_70_180","cv","B","Time between 70-180 mg/dL (%)","Coefficient of variation")
sfig9c <- bivariate_plot(tab2_df,"lbgi","hbgi","C","Low Blood Glucose Index","High Blood Glucose Index")
sfig9d <- bivariate_plot(tab2_df,"grade_hypo","grade_hyper","D","%GRADE hypo","%GRADE hyper")

sfig9 <- ggpubr::ggarrange(sfig9a,sfig9b,
                           sfig9c,sfig9d,
                           ncol = 2,nrow=2,
                           legend="bottom",
                           common.legend=TRUE)

sfig9
saveRDS(sfig9, paste0(path_cgm_repo,"/paper3/output/sfig9.RDS"))
if(Sys.info()["user"]=="JVARGH7"){
  library(GGally)
  sfig9 <- readRDS(paste0(path_cgm_repo,"/paper3/output/sfig9.RDS"))
  sfig9
  tiff(paste0(path_cgm_repo,"/paper3/output/sfig9.tif"),res=300,width = 12,height=8,units = "in")
  sfig9
  dev.off()
}