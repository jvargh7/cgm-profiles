
library(NMF )
paper3_df <- readRDS(paste0(path_cgm_working,"/paper3/paper3_df.RDS")) 
nmf_output <- readRDS(paste0(path_cgm_working,"/paper3/nmf_comparison_output.RDS"))
nmf_output2 <- readRDS(paste0(path_cgm_working,"/paper3/nmf_output.RDS"))
source(paste0(path_cgm_repo,"/paper3/paper3_table2_nmf strata.R"))

mdrf_anonymized = readRDS(paste0(path_cgm_datasets,"/mdrf_anonymized.RDS")) %>% 
  dplyr::filter(file_name %in% paper3_df$file_name)
cluster_colors = c("grey","black","black")
cluster_fills = c("black","grey","white")

# sfig5_df <- mdrf_anonymized %>% 
#   left_join(tab2_df %>% 
#               dplyr::select(cluster,file_name),
#             by = "file_name") %>% 
#   group_by(cluster,d_hour) %>% 
#   summarize(median = median(sensorglucose),
#             q25 = quantile(sensorglucose,probs=0.25),
#             q75 = quantile(sensorglucose,probs=0.75)
#             )
# 
# sfig5_df %>% 
#   ungroup() %>% 
#   dplyr::filter(!is.na(cluster)) %>% 
#   ggplot(data=.,aes(x=d_hour,y=median,ymin=q25,ymax=q75,group=cluster,color=cluster)) +
#   geom_point() +
#   geom_errorbar() +
#   theme_bw()

sfig5_df <- mdrf_anonymized %>%
  left_join(tab2_df %>%
              dplyr::select(cluster,file_name),
            by = "file_name") %>%
  group_by(cluster,file_name,d_hour) %>%
  summarize(median = median(sensorglucose),
            q25 = quantile(sensorglucose,probs=0.25),
            q75 = quantile(sensorglucose,probs=0.75)
            )


source(paste0(path_cgm_repo,"/package/plot_agp.R"))
set.seed(1980)
#Cluster 1
r3_c1 <- tab2_df %>%
  dplyr::filter(cluster == "TIR profile"
  ) %>% 
  sample_n(.,size=1) %>% 
  ungroup()
r3_c2 <- tab2_df %>%
  dplyr::filter(cluster == "Hypo profile"
  ) %>% 
  sample_n(.,size=1) %>% 
  ungroup()
r3_c3 <- tab2_df %>%
  dplyr::filter(cluster == "Hyper profile"
  ) %>% 
  sample_n(.,size=1) %>% 
  ungroup()

sfig5a <- plot_agp_bw(r3_c1,mdrf_anonymized,facet_names = "") +
  ggtitle("A") +
  theme(
    # text = element_text(size=15),
    legend.position = "none") +
  scale_y_continuous(breaks=seq(0,300,by=50),
                     limits=c(0,300))

sfig5b <-  plot_agp_bw(r3_c2,mdrf_anonymized,facet_names = "") +
  ggtitle("B") +
  theme(
        # text = element_text(size=15),
        legend.position = "none") +
  scale_y_continuous(breaks=seq(0,300,by=50),
                     limits=c(0,300))

sfig5c <-  plot_agp_bw(r3_c3,mdrf_anonymized,facet_names = "") +
  ggtitle("C") +
  theme(
    # text = element_text(size=15),
    legend.position = "none") +
  scale_y_continuous(breaks=seq(0,300,by=50),
                     limits=c(0,300))

sfig5d = sfig5_df %>%
  ungroup() %>%
  dplyr::filter(!is.na(cluster)) %>%
  ggplot(data=.,aes(x=factor(d_hour),y=median)) +
  geom_boxplot(aes(fill=cluster, color = cluster),outlier.colour = NA) +
  # geom_errorbar() +
  theme_bw(base_size = 12) +
  xlab("Hour (24 hour format)") +
  ylab("Median (mg/dL)") +
  theme(legend.position = "bottom",
        # legend.text = element_text(size = 12),
        # axis.title = element_text(size = 15),
        # axis.text = element_text(size=12)
        ) +
  scale_fill_manual(name="", values = cluster_fills) +
  scale_color_manual(name="", values = cluster_colors) +
  scale_y_continuous(breaks=seq(0,300,by=50),
                     limits=c(0,300)) +
  ggtitle("D")

sfig5 <- ggpubr::ggarrange(sfig5a,
                           sfig5b,
                           sfig5c,
                           sfig5d,
                          ncol = 2,nrow=2)
sfig5
saveRDS(sfig5,paste0(path_cgm_repo,"/paper3/output/sfig5.RDS"))

ggsave(paste0(path_cgm_repo,"/paper3/output/sfig5_mdrf.png"),height=8,width=12,units="in",dpi = 300)


if(Sys.info()["user"]=="JVARGH7"){
sfig5 <- readRDS(paste0(path_cgm_repo,"/paper3/output/sfig5.RDS"))
sfig5
tiff(paste0(path_cgm_repo,"/paper3/output/sfig5.tif"),res=300,width = 12,height=8,units = "in")
sfig5
dev.off()
}
