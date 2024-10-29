#----Figure_04D_S14-------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 04 Panel D and S14 for the connectivity analysis 
# for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_output_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_preparations.py")
  source("./scripts/01_setup.R")
}
synapses_curated =  read_csv(paste0(PATH_input,"tmp/clk_neuron_output_filtered_v",v,".csv"),
                      col_types =  cols(pre_pt_supervoxel_id = col_character(),
                                        pre_pt_root_id = col_character(),
                                        post_pt_supervoxel_id = col_character(),
                                        post_pt_root_id = col_character()))
if (!paste0("classification_v",v,".csv") %in% input_files) {
  stop("please go to https://codex.flywire.ai/api/download n\ 
       and dowload the classification file and save it in './input'.")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
#-------------------------------------------------------------------------------
synapses_sum_output = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
clk_output_sum = synapses_sum_output[synapses_sum_output$n_synapses >= 5,]
                       
clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
clk_output_sum = left_join(clk_output_sum,clk_join,by = "pre_pt_root_id")
colnames(clk_join) = c("col","post_pt_root_id","hemisphere_post")
clk_output_sum = left_join(clk_output_sum,clk_join,by = "post_pt_root_id")
cluster = unique(clk$clk_name)

for (i in cluster) {
  clk_output = clk_output_sum[clk_output_sum$pre_pt_root_id%in% 
                          clk[clk$clk_name %in% i,]$clk_id,]
  clk_output_ordered = clk_output%>%
    group_by(post_pt_root_id)%>%
    summarise(n_total = sum(n_synapses))
  clk_output_ordered = clk_output_ordered[order(
    clk_output_ordered$n_total,decreasing = T),]
  clk_output_ordered$position = c(paste0("000",c(1:9)),paste0("00",c(10:99)),
                                  paste0("0",c(100:999)),
                                  1000:9999)[1:length(clk_output_ordered$post_pt_root_id)]
  clk_output = left_join(clk_output,clk_output_ordered[,c("post_pt_root_id","position"),
                                                       by="post_pt_root_id"])
  clk_output$output_id_pos = paste(clk_output$position,clk_output$post_pt_root_id,sep="_")
  ggplot(clk_output,aes(x=output_id_pos, y=n_synapses,fill= col))+ 
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
  if (i %in% c("DN1pA","s-CPDN3C","s-CPDN3D")) {
    ggsave(filename = paste0(PATH_output,"Figure_04/Figure_04D_",i,"_output.pdf"),
           width = 40,height = 20,units = "cm")
  }else{
    ggsave(filename = paste0(PATH_output,"Figure_S14/",i,"_output.pdf"),width = 40,
           height = 20,units = "cm")
  }
}

