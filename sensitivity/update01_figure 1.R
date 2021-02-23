
paper3_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) %>% 
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

library(GGally)

# https://stackoverflow.com/questions/57450913/removing-background-color-for-column-labels-while-keeping-plot-background-color/57455317#57455317
cor_func <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"), 
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  
  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  ) + #removed theme_void()
    # theme(panel.background = element_rect(fill = fill))
    theme(panel.background = element_rect(fill = NA))
}

agp_fig_revised <- paper3_df %>% 
  dplyr::select(average_sensor,gmi,cv,
                percent_time_over_250,
                percent_time_180_250,
                percent_time_70_180,
                percent_time_54_70,
                percent_time_under_54,
                percent_time_over_180,
                percent_time_under_70) %>%
  ggpairs(.,
          columnLabels = c("Average \n(mg/dL)", "GMI", "CV",
                           "Time > 250 \nmg/dL (%)",
                           "Time 180-250 \nmg/dL (%)",
                           "Time 70-180 \nmg/dL (%)",
                           "Time 54-69 \nmg/dL (%)",
                           "Time < 54 \nmg/dL (%)",
                           "Time > 180 \nmg/dL (%)",
                           "Time < 70 \nmg/dL (%)"
          ), 
          upper = list(continuous = wrap(cor_func,
                                         method = 'pearson', symbol = expression('\u03C1 ='))),
          lower = list(continuous = function(data, mapping, ...) {
            # ggally_smooth(data = data, mapping = mapping) +
            #   geom_point(size=0.2) +
            #   theme(panel.background = element_blank())
            ggplot(data=data,mapping=mapping) +
              geom_point(color="grey",alpha=0.3,size=0.5) +
              geom_smooth(color="black",method="loess",size=0.8) +
              theme(panel.background = element_blank())
            
          }),
          diag = list(continuous = function(data, mapping, ...) {
            ggally_barDiag(data = data, mapping = mapping) + 
              theme(panel.background = element_blank())}
          ))

mytheme = theme(strip.background = element_rect(fill = "white"),
                strip.text = element_text(color = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())


agp_fig_revised <- agp_fig_revised + mytheme
agp_fig_revised
saveRDS(agp_fig_revised,paste0(path_cgm_repo,"/paper3/output/agp_fig_revised.RDS"))


if(Sys.info()["user"]=="JVARGH7"){
  library(GGally)
  agp_fig_revised <- readRDS(paste0(path_cgm_repo,"/paper3/output/agp_fig_revised.RDS"))
  agp_fig_revised
  tiff(paste0(path_cgm_repo,"/paper3/output/agp_fig_revised.tif"),res=300,width = 12,height=8,units = "in")
  agp_fig_revised
  dev.off()
}


# agp_fig <- paper3_df %>%
#   dplyr::select(average_sensor,gmi,cv,
#                 percent_time_over_250,
#                 percent_time_180_250,
#                 percent_time_70_180,
#                 percent_time_54_70,
#                 percent_time_under_54) %>%
#   # # http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs
#   # pairs.panels(.,
#   #              method = "pearson", # correlation method
#   #              hist.col = "#00AFBB",
#   #              density = TRUE,  # show density plots
#   #              ellipses = TRUE # show correlation ellipses
#   # )
#   ggpairs(.,
#           # columns = c("average_sensor", "gmi", "cv",
#           #               "percent_time_over_250",
#           #               "percent_time_180_250",
#           #               "percent_time_70_180",
#           #               "percent_time_54_70",
#           #               "percent_time_under_5"),
#           columnLabels = c("Average (mg/dL)", "GMI", "CV",
#                            "Time > 250 mg/dL (%)",
#                            "Time 180-250 mg/dL (%)",
#                            "Time 70-180 mg/dL (%)",
#                            "Time 54-69 mg/dL (%)",
#                            "Time < 54 mg/dL (%)"
#                            ))


