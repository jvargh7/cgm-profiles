nmf_output <- readRDS(paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))
coef <- coef(nmf_output$fit$`3`)

coef_df <- coef %>% 
  data.frame() %>% 
  mutate(cluster = rownames(.)) %>% 
  mutate(cluster = factor(cluster,levels=c(3,1,2), 
                          labels=c("TIR profile","Hypo profile","Hyper profile"),
                          ordered = TRUE)) %>% 
  
  mutate(cluster_hypofirst = fct_relevel(cluster,"Hypo profile","TIR profile")) %>% 
  pivot_longer(cols=-c(contains("cluster")),names_to="var",values_to = "values") %>% 
  mutate(var = case_when(var == "average_sensor" ~ 1,
                         var == "gmi" ~ 2,
                         var == "cv" ~ 3,
                         var == "percent_time_over_250" ~ 4,
                         var == "percent_time_180_250" ~ 5,
                         var == "percent_time_70_180" ~ 6,
                         var == "percent_time_54_70" ~ 7,
                         var == "percent_time_under_54" ~ 8,
  ),
  var = factor(var, levels=c(1:8),labels=c("Average of sensor values",
                                           "GMI",
                                           "Coefficient of variation",
                                           "Percentage of time over 250 mg/dL",
                                           "Percentage of time between 181-250 mg/dL",
                                           "Percentage of time between 70-180 mg/dL",
                                           "Percentage of time between 54-69 mg/dL",
                                           "Percentage of time under 54 mg/dL"
  )),
  values =  round(values,2),
  # values = case_when(values > 1 ~ 1,
  #                    TRUE ~ round(values,2)),
  values_2d = cut(values,breaks=seq(0,1.1,by=0.1),
                  # labels = c("0-0.10","")
                  include.lowest = TRUE))
saveRDS(coef_df,paste0(path_cgm_repo,"/paper3/coef_df.RDS"))

sfig2 <- ggplot(data = coef_df,aes(x= var,y=cluster_hypofirst,fill=values_2d,label=values)) +
  geom_tile() +
  geom_text() + 
  scale_fill_grey("Loadings",start = 0.9,end=0.5)+
  theme_bw() +
  xlab("")+
  ylab("") +
  theme(axis.text.x = element_text(angle = 60,hjust =0.95))

sfig2

