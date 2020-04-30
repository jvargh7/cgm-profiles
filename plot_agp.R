
plot_agp <- function(cgm_demo,cgm_df,
                     ts_freq = (24*60/15),
                     low_range_cutoff = 70,
                     high_range_cutoff = 180,
                     hypo_cutoff = 54,
                     hyper_cutoff = 250,
                     aggregate_h = NULL,
                     additional_vars = NULL){
  
  cgm_df <- cgm_df %>% 
    dplyr::select(timestamp_id,file_name,sensorglucose) %>% 
    dplyr::filter(file_name %in% c(cgm_demo$file_name))
  
  if(!is.null(additional_vars)) {
    
    cgm_df <- cgm_df %>% 
      left_join(cgm_demo %>% 
                  dplyr::select(additional_vars,file_name),
                by = "file_name")
  }
  
  p <- cgm_df %>% 
    
    # spread(key=file_name,value=sensorglucose) %>% 
    # mutate_at(vars(-contains("timestamp_id")),
    #           .f=function(x) ts(x,frequency = ts_freq)) %>% 
    
    # group_by(file_name) %>%
    # mutate(sensorglucose = ts(sensorglucose,
    #                           frequency = ts_freq)) %>%
    
    ggplot(data=.,aes(x=timestamp_id,
                      y=sensorglucose,
                      group=file_name,
                      col=file_name)) +
    geom_path()
  
  if(!is.null(aggregate_h)){
    
    
    agg_colors = c("grey25","grey50","grey75")
    agg_linetypes = c("dotted","dashed","twodash")
    for(i in 1:length(aggregate_h)){
      aggregated_df <- cgm_df %>% 
        dplyr::select(timestamp_id,file_name,sensorglucose) %>% 
        dplyr::filter(file_name %in% c(cgm_demo$file_name))  
      aggregated_df <- aggregate_cgm(id = aggregated_df$file_name,
                                     timestamp_id=aggregated_df$timestamp_id,
                                     value=aggregated_df$sensorglucose,
                                     ts_freq = ts_freq,
                                     period_hours = aggregate_h[[i]]) %>% 
        gather(key="file_name",value="sensorglucose",-timestamp_id,-period_id)
      
      p <- p + geom_path(data=aggregated_df,aes(x=timestamp_id,
                                                y=sensorglucose,
                                                group=file_name
      ),col=agg_colors[[i]],size=1.2,linetype = agg_linetypes[[i]])
      
    }
                        
                       
  }
  
  
  p <- p + facet_grid(file_name~.) +
    theme_bw() +
    xlab("Timestamp") +
    ylab("Sensor Glucose (15-min readings)") +
    scale_x_continuous(limits=c(0,1400),
                       breaks = seq(0,1400,by=96)) + 
    annotate("rect", xmin = 0, xmax = 1400, 
             ymin = low_range_cutoff, ymax = high_range_cutoff, 
             fill = "lightgreen", alpha = 0.1, color = "darkgreen") +
    geom_hline(yintercept = c(hypo_cutoff,hyper_cutoff),
               color="orange") + 
    geom_vline(xintercept=seq(0,1400,by=96),col="grey",linetype="dashed")
  
  if(sum(cgm_df$sensorglucose>hyper_cutoff) > 0){
    p <- p + geom_ribbon(aes(ymin = hyper_cutoff, 
                             ymax = ifelse(sensorglucose>hyper_cutoff,
                                           sensorglucose,
                                           NA)), 
                         fill = "red", alpha = .1, color = "orange")
    
  }
  
  if(sum(cgm_df$sensorglucose < hypo_cutoff) > 0){
    p <- p + geom_ribbon(aes(ymax = hypo_cutoff, 
                             ymin = ifelse(sensorglucose < hypo_cutoff,
                                           sensorglucose,
                                           NA)), 
                         fill = "red", alpha = .1, color = "orange")
    
  }
     
    
  
  
  
  return(p)
  
  
}

## Aggregation of CGM data -----

aggregate_cgm <- function(id,timestamp_id,value,ts_freq = 96,period_hours = 3,
                          df_format = "wide",col_names=NULL){
  
  if(period_hours >= 0.5){
    start_index = period_hours*(60/15)
    by_index = start_index
    
    
    
    df <- data.frame(id = id,
                     timestamp_id = timestamp_id,
                     value = value) %>% 
      ungroup() %>% 
      dplyr::arrange(timestamp_id) %>% 
      spread(key=id,value=value) %>% 
      map_df(.,.f=function(x){
        x <- ts(x,frequency = ts_freq);
        period.apply(x, INDEX=seq(start_index,length(x),by=by_index), FUN=function(x) mean(x))
      }) %>% 
      mutate(period_id = 1:nrow(.)) %>% 
      mutate_at(vars(-contains("id")),.f=function(x) ts(x,frequency = 24/period_hours))
    
    if(df_format == "long"){
      df <- df %>% 
        dplyr::select(-timestamp_id) %>% 
        gather(key = "id",value="value",-period_id)
      
      if(!is.null(col_names)){
        names(df) = col_names
      }
      
    }
    
    
    return(df)
    
  } else(return(NA))
  
  
  
  
}
