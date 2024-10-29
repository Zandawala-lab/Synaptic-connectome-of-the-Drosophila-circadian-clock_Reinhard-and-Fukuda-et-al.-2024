#----Figure_S13-----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S13 for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
  source("./scripts/01_setup.R")
}
synapses_curated = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_supervoxel_id = col_character(),
                                             pre_pt_root_id = col_character(),
                                             post_pt_supervoxel_id = col_character(),
                                             post_pt_root_id = col_character()))

if (!paste0("classification_v",v,".csv") %in% input_files) {
  stop("please go to https://codex.flywire.ai/api/download n\ 
    and dowload the clssification file for the current version and save it in './input'.")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
#-------------------------------------------------------------------------------
clk_input_sum = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
clk_input_sum = clk_input_sum[clk_input_sum$n_synapses >=  5,]
clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_input_sum = left_join(clk_input_sum,clk_join,by = "post_pt_root_id")
colnames(clk_join) = c("col","pre_pt_root_id","hemisphere_pre")
clk_input_sum = left_join(clk_input_sum,clk_join[c("col","pre_pt_root_id")],by = "pre_pt_root_id")
cluster = unique(clk$clk_name)

for (i in cluster) {
  clk_input = clk_input_sum[clk_input_sum$post_pt_root_id%in%clk[clk$clk_name %in% i,]$clk_id,]
  clk_input_ordered = clk_input%>%
    group_by(pre_pt_root_id)%>%
    summarise(n_total = sum(n_synapses))
  clk_input_ordered = clk_input_ordered[order(clk_input_ordered$n_total,decreasing = T),]
  clk_input_ordered$position = c(paste0("000",c(1:9)),paste0("00",c(10:99)),
                                 paste0("0",c(100:999)),
                                 1000:9999)[1:length(clk_input_ordered$pre_pt_root_id)]
  clk_input = left_join(clk_input,clk_input_ordered[,c("pre_pt_root_id","position"),
                                                    by="pre_pt_root_id"])
  clk_input$input_id_pos = paste(clk_input$position,clk_input$pre_pt_root_id,sep="_")
  
  ggplot(clk_input,aes(x=input_id_pos, y=n_synapses,fill= col))+ 
    geom_col()+
    scale_fill_manual(values = c("DN1a"=DN1a,"DN1pA"=DN1p_A,"DN1pB"=DN1p_B,
                                 "DN1pC"=DN1p_C,"DN1pD"=DN1p_D,"DN1pE"=DN1p_E,
                               "DN2"=DN2,"APDN3"=APDN3,"l-CPDN3"=l_CPDN3,
                               "s-CPDN3A"=s_CPDN3A,"s-CPDN3B"=s_CPDN3B,
                               "s-CPDN3C"=s_CPDN3C,"s-CPDN3D"=s_CPDN3D,
                               "s-CPDN3E"=s_CPDN3E,"LPN"=LPN,"LNd_CRYp"=LNd_CRYp,
                               "LNd_CRYn"=LNd_CRYn,"LN_ITP"=LN_ITP,"s-LNv"=s_LNv,
                               "l-LNv"=l_LNv),na.value = "lightgrey")+
    geom_hline(yintercept = 0)+
    theme(strip.text.x = element_text(size = 6),
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, vjust = 0, hjust=1,size = 5),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      strip.background =element_rect(fill="white"),
      panel.background = element_rect(fill = 'white'),
      panel.spacing = unit(20, "points")
    )
  ggsave(filename = paste0(PATH_output,"Figure_S13/",i,"_input.pdf"),width = 40,
         height = 20,units = "cm")
}