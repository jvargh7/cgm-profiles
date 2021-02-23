date = "2021-02-07"


tab1_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
  dplyr::filter(!is.na(hba)) %>% 
  group_by(patient_id) %>% 
  arrange(record_date) %>% 
  
  dplyr::mutate(nvisit = n()) %>% 
  dplyr::mutate(nvisit = as.character(nvisit)) %>% 
  dplyr::filter(record_date == min(record_date)) %>% 
  ungroup()

cgm_hourly_summary_t1df <- readRDS(paste0(path_cgm_working,"/cgm_output/mdrf_cgm_hourly_",date,".RDS")) %>% 
  dplyr::filter(file_name %in% tab1_df$file_name) %>% 
  left_join(tab2_df %>% 
              dplyr::select(file_name,cluster),
            by = "file_name") %>% 
  mutate(cluster_hypofirst = fct_relevel(cluster,"Hypo profile","TIR profile"))

sfig10a <- cgm_hourly_summary_t1df %>% 
  group_by(cluster_hypofirst,hours_range_min) %>% 
  rename(GRpercent70_180 = 'GRpercent70 - 180',
         GRpercent_under54 = 'GRpercent0 - 54',
         GRpercent54_69 = 'GRpercent54 - 69',
         GRpercent180_250 = 'GRpercent181 - 250',
         GRpercent_over250 = 'GRpercent250 - 1000'
  ) %>% 
  summarize_at(vars(contains("GRpercent")),~mean(.,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(hypo = GRpercent_under54 + GRpercent54_69,
         hyper = GRpercent180_250 + GRpercent_over250) %>% 
  ggplot(data=.,aes(x=factor(hours_range_min),fill=hypo,y=cluster_hypofirst)) +
  geom_tile() +
  scale_fill_continuous(low="white",high="red") +
  labs(x="24 hour format",y="",title = "A",fill = "")

sfig10b <- cgm_hourly_summary_t1df %>% 
  group_by(cluster_hypofirst,hours_range_min) %>% 
  rename(GRpercent70_180 = 'GRpercent70 - 180',
         GRpercent_under54 = 'GRpercent0 - 54',
         GRpercent54_69 = 'GRpercent54 - 69',
         GRpercent180_250 = 'GRpercent181 - 250',
         GRpercent_over250 = 'GRpercent250 - 1000'
  ) %>% 
  summarize_at(vars(contains("GRpercent")),~mean(.,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(hypo = GRpercent_under54 + GRpercent54_69,
         hyper = GRpercent180_250 + GRpercent_over250) %>% 
  ggplot(data=.,aes(x=factor(hours_range_min),fill=hyper,y=cluster_hypofirst)) +
  geom_tile() +
  scale_fill_continuous(low="white",high="red") +
  labs(x="24 hour format",y="",title = "B",fill = "")

sfig10 <- ggpubr::ggarrange(sfig10a,sfig10b,
                           ncol = 1,nrow=2,
                           legend="right",
                           common.legend=TRUE)

sfig10
