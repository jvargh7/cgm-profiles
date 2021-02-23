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
              dplyr::select(file_name,starts_with("grade"),G),
            by = "file_name") %>% 
  mutate_at(vars(starts_with("grade"),G),~as.numeric(.))


tab2_df %>% 
  compareGroups::compareGroups(cluster~G + hbgi + lbgi + grade_hyper + grade_hypo,data=.,
                               method = 2) %>% 
  compareGroups::createTable(.,show.all=TRUE,show.p.overall = FALSE,q.type=c(2,2)) %>% 
  compareGroups::export2csv(.,paste0(path_cgm_repo,"/paper 3 review/table1_addl.csv"))
