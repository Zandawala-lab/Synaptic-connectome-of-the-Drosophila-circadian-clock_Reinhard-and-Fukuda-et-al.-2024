#----Figure_S04_S05-------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S04/ S05 for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
 reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
}

clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
          col_types = cols(clk_id = col_character()),delim = ",")

synapses_curated = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
               col_types = cols(pre_pt_supervoxel_id = col_character(),
                         pre_pt_root_id = col_character(),
                         post_pt_supervoxel_id = col_character(),
                         post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
synapses_sum = synapses_curated %>% 
 group_by(pre_pt_root_id, post_pt_root_id) %>% 
 summarise(n_synapses = length(post_pt_root_id))
synapses_sum_clk = synapses_sum[synapses_sum$pre_pt_root_id%in%clk$clk_id,]

clk_join = clk
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "pre_pt_root_id")

colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "post_pt_root_id")

clk_order = c("s-LNv","l-LNv","LN_ITP","LNd_CRYp",
       "LNd_CRYn","LPN","l-CPDN3","APDN3","s-CPDN3A",
       "s-CPDN3B","s-CPDN3C","s-CPDN3D","s-CPDN3E",
       "DN2","DN1pA","DN1pB","DN1pC","DN1pD","DN1pE","DN1a")
clk_level = c(1:length(clk_order)) 

synapses_sum_clk$name_pre_order = synapses_sum_clk$name_pre
synapses_sum_clk$name_post_order = synapses_sum_clk$name_post

for (i in clk_level) {
 synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i],]$name_pre_order = 
  paste(letters[i],
        synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i],]$name_pre,sep = "_")
 synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i],]$name_post_order = 
  paste(letters[i],
        synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i],]$name_post,sep = "_")
}

synapses_sum_clk$pre_name_hemisphere = paste(synapses_sum_clk$hemisphere_pre,
                      synapses_sum_clk$name_pre_order,
                      synapses_sum_clk$pre_pt_root_id,sep = "_")
synapses_sum_clk$post_name_hemisphere = paste(synapses_sum_clk$hemisphere_post,
                       synapses_sum_clk$name_post_order,
                       synapses_sum_clk$post_pt_root_id,sep = "_")

hemispheres = c("left","right")
filenames = c("Figure_S04", "Figure_S05")

for (i in c(1:2)) {
  synapses_sum_clk_grouped = synapses_sum_clk[synapses_sum_clk$hemisphere_pre == hemispheres[i],] %>% 
    group_by(post_name_hemisphere,pre_name_hemisphere,hemisphere_post) %>% 
    summarise(n_synapses_sum = sum(n_synapses,na.rm = T),
              n_pre_partners = length(unique(pre_pt_root_id)),
              n_post_partners = length(unique(post_pt_root_id)))
  
  synapses_sum_clk_grouped$col = ifelse(synapses_sum_clk_grouped$n_synapses_sum >= 5,
                                        "white","black")

  ggplot(synapses_sum_clk_grouped[synapses_sum_clk_grouped$n_synapses_sum >= 1,],
         aes(x = pre_name_hemisphere,y = post_name_hemisphere,fill = n_synapses_sum),
         color = "black")+
    geom_tile()+
    geom_text(aes(label = n_synapses_sum,colour = col),show.legend = FALSE,
              size = (3.5/(1/0.35)))+ # 1/0.35 is the ratio between points used by theme() and mm used by geom_text()
    scale_colour_manual(values = c("black","white"))+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    labs(y = "postsynaptic",x = "presynaptic")+
    scale_fill_gradientn(colors = grey_scale(length(synapses_sum_clk_grouped$n_synapses_sum)),
                         trans = "log",na.value = "white")+
    geom_vline(xintercept = c(0,
        seq(1.5,(length(unique(synapses_sum_clk_grouped$pre_name_hemisphere))+0.5),by = 1)),
        linewidth = 0.01)+
    geom_hline(yintercept = c(108.5),linewidth = 0.5)+
    geom_hline(yintercept = c(0,
        seq(1.5,(length(unique(synapses_sum_clk_grouped$post_name_hemisphere))+0.5),by = 1)),
        linewidth = 0.01)+
    theme(axis.text.x = element_text( angle = 90,hjust = 0,vjust = 0.5, colour = "black"),
          axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
          text = element_text(size = 5),
          legend.position = "",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",colour = "black")
    )
  ggsave(paste0(PATH_output,"Figure_S04_05/",filenames[i],"_v",v,".pdf"), 
         width = 184,height = 244,units = c("mm"))
  
}
