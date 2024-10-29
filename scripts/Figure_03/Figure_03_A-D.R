#----Figure_03_B and analysis as basis for A,C,D--------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 03 Panel B for the connectivity analysis for 
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
  # the use of other classification files might change the results slightly 
  # depending on the changes in classification. 
  # The classification file used in the paper was updated last 26-12-2023
 stop("please go to https://codex.flywire.ai/api/download n\ 
    and dowload the classification file for the current version and save it in './input'.")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
          col_types = cols(clk_id = col_character()),delim = ",")

classification = read_delim(paste0(PATH_input,"classification_v",v,".csv"), 
               delim = ",",
               escape_double = FALSE,
               col_types = cols(root_id = col_character(),
                       flow = col_character()),
               trim_ws = TRUE)
classification = as.data.frame(classification)
#-------------------------------------------------------------------------------
#updating classification--------------------------------------------------------
classification[classification$hemibrain_type%in%c("HBeyelet"),]$super_class = "sensory"
#-------------------------------------------------------------------------------
#Figure_03_B--------------------------------------------------------------------
clk_input_sum = synapses_curated %>%
 group_by(pre_pt_root_id, post_pt_root_id) %>%
 summarise(n_synapses = length(post_pt_root_id))
write.csv(clk_input_sum, paste0(PATH_input,"tmp/clk_input_sum_filtered_v",v,".csv"))

classification_join = classification
colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
clk_input = left_join(clk_input_sum[clk_input_sum$n_synapses >= 5,],
           classification_join,by = "pre_pt_root_id")
clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_input = left_join(clk_input,clk_join,by = "post_pt_root_id")
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
clk_input = left_join(clk_input,clk_join,by = "pre_pt_root_id")
clk_input_total = clk_input %>%
 group_by(name_post) %>%
 summarise(n_synapses_total = sum(n_synapses,na.rm = T),
      n_pre_partners_total = length(unique(pre_pt_root_id)))
#-------------------------------------------------------------------------------
clk_input$clk_marker = "no_clock"
clk_input[clk_input$pre_pt_root_id%in%clk$clk_id,]$clk_marker = "clock"

clk_input$post_pt_root_id_name = paste(clk_input$post_pt_root_id,
                    clk_input$name_post,sep = "_")

clk_input_sum_grouped = clk_input %>%
 group_by(
  clk_marker,
  super_class,
  name_post
 )%>%
 summarise(n_synapses_sum = sum(n_synapses,na.rm = T),
      avrg_synapses = mean(n_synapses,na.rm = T),
      n_pre_partners = length(unique(pre_pt_root_id)),
      n_post_partners = length(unique(post_pt_root_id)))

clk_input_sum_grouped = left_join(clk_input_sum_grouped,clk_input_total,
                 by = "name_post")
clk_input_sum_grouped$perc_of_input = clk_input_sum_grouped$n_synapses_sum/
                   clk_input_sum_grouped$n_synapses_total

clk_input_sum_grouped$name_post = factor(clk_input_sum_grouped$name_post,
                   levels = rev(c("s-LNv","l-LNv","LN_ITP",
                          "LNd_CRYp","LNd_CRYn","LPN",
                          "l-CPDN3","APDN3","s-CPDN3A",
                          "s-CPDN3B","s-CPDN3C",
                          "s-CPDN3D","s-CPDN3E","DN2",
                          "DN1pE","DN1pD","DN1pC",
                          "DN1pB","DN1pA","DN1a"))
)
clk_input_sum_grouped$fill = paste(clk_input_sum_grouped$super_class,
                clk_input_sum_grouped$clk_marker,sep = "_")
pattern_density = unique(clk_input_sum_grouped$fill)
clk_input_sum_grouped$super_class_corr = clk_input_sum_grouped$super_class
clk_input_sum_grouped$super_class_corr[clk_input_sum_grouped$perc_of_input<0.025] = "others"
clk_input_sum_grouped$super_class_corr[is.na(clk_input_sum_grouped$super_class)] = "undefined"
clk_input_sum_grouped$fill_corr = paste(clk_input_sum_grouped$super_class_corr,
                  clk_input_sum_grouped$clk_marker,sep = "_")
# in case the script is used for updated data there might be the need to update
# the factor level accordingly (the same is true for " scale_pattern_density_manual" 
# and "scale_fill_manual" in ggplot)
clk_input_sum_grouped = clk_input_sum_grouped%>%
 group_by(name_post,fill_corr)%>%
 summarize(perc_of_input = sum(perc_of_input))
clk_input_sum_grouped$fill_corr = factor(clk_input_sum_grouped$fill_corr,
                  levels = rev(c("others_clock",
                          "others_no_clock",
                          "undefined_no_clock",
                          "optic_no_clock",
                          "descending_no_clock",
                          "visual_centrifugal_no_clock",
                          "central_clock",
                          "central_no_clock",
                          "visual_projection_clock",
                          "visual_projection_no_clock",
                          "sensory_no_clock"
                          )))

ggplot(clk_input_sum_grouped,aes(x = "",y = perc_of_input,fill = fill_corr))+
 geom_bar_pattern(stat = "identity",col = "black",linewidth = 0.1,pattern = "stripe",
          aes(pattern_density = fill_corr),pattern_fill = "white",
          pattern_colour = NA,pattern_angle = 0,pattern_size = 0.1,
          pattern_alpha = 1,pattern_key_scale_factor = 0.3)+
 coord_polar("y", start = 0)+
 facet_wrap(name_post~.,ncol = 7)+
 scale_pattern_density_manual(values = 
                 c("others_clock" = 0.1,
                  "undefined_no_clock" = 0,
                  "central_clock" = 0.1,
                  "others_no_clock" = 0,
                  "visual_projection_clock" = 0.1,
                  "central_no_clock" = 0,
                  "descending_no_clock" = 0,
                  "optic_no_clock" = 0,
                  "visual_centrifugal_no_clock" = 0,
                  "visual_projection_no_clock" = 0,
                  #"ascending_no_clock" = 0,
                  "sensory_no_clock" = 0
                 ))+
 
 scale_x_discrete(expand = c(0,0))+
 scale_fill_manual(values = c("central_clock" = "#70B04F",
                "visual_projection_clock" = "#FFC800",
                "others_clock" = "black",
                "undefined_no_clock" = "white",
                "central_no_clock" = "#70B04F",
                "descending_no_clock" = "cyan",
                "optic_no_clock" = "#999999",
                "visual_centrifugal_no_clock" = "#1B458F",
                "visual_projection_no_clock" = "#FFC800",
                "others_no_clock" = "black",
                "sensory_no_clock" = "magenta"))+
 ylab(label = "")+
 xlab(label = "proportion of input synapses")+
 theme(panel.background = element_rect(fill = NA,color = NA),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0,0,0,0), "mm"),
    panel.border = element_rect(colour = "black", fill = NA,
                   linewidth = NA),
    legend.title = element_blank(),
    text = element_text(size = 8, colour = "black"),
    axis.text = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    panel.spacing = unit(0, "line"),
    legend.position = "bottom",
    axis.ticks = element_blank()
    
 )
ggsave(paste0(PATH_output,"Figure_03/Figure_03_B_v",v,".pdf"), width = 21,
    height = 8, units = "cm")
#-------------------------------------------------------------------------------
#analysis for Figure_03_A-------------------------------------------------------
cell_count = clk_input%>%
 group_by(super_class)%>%
 summarize(n_superclass = length(unique(pre_pt_root_id)),
      n_clock = length(unique(post_pt_root_id)))
#-------------------------------------------------------------------------------
#analysis for Figure_03_C,D-----------------------------------------------------
col_to_keep =c("pre_pt_root_id","cell_type","class","hemibrain_type","name_pre",
        "name_post","hemisphere_post","post_pt_root_id","n_synapses","super_class")
clk_input = clk_input[,col_to_keep]
clk_input_visual_proj = clk_input[clk_input$super_class %in% c("visual_projection"),]
clk_input_visual_centri = clk_input[clk_input$super_class %in% c("visual_centrifugal"),]
clk_input_sensory = clk_input[clk_input$super_class %in% ("sensory"),]
clk_input_optic = clk_input[clk_input$super_class %in% c("optic"),]
clk_input_endocrine = clk_input[clk_input$super_class %in% c("endocrine"),]
clk_input_descending = clk_input[clk_input$super_class %in% c("descending"),]
clk_input_central = clk_input[clk_input$super_class %in% c("central"),]
clk_input_ascending = clk_input[clk_input$super_class %in% c("ascending"),]
clk_input_NA = clk_input[is.na(clk_input$super_class),]

write.csv(clk_input_visual_proj,paste0(PATH_output,
    "Figure_03/Figure_03_C_input_to_clk_from_visual_projection_",v,".csv"))
write.csv(clk_input_visual_centri,paste0(PATH_output,
    "Figure_03/Figure_03_C_input_to_clk_from_visual_centrifugal_",v,".csv"))
write.csv(clk_input_optic,paste0(PATH_output,
    "Figure_03/Figure_03_C_input_to_clk_from_optic_",v,".csv"))
write.csv(clk_input_central,paste0(PATH_output,
    "Figure_03/Figure_03_C_input_to_clk_from_central_",v,".csv"))
write.csv(clk_input_sensory,paste0(PATH_output,
    "Figure_03/Figure_03_D_input_to_clk_from_sensory_",v,".csv"))
write.csv(clk_input_descending,paste0(PATH_output,
    "Figure_03/input_to_clk_from_descending_",v,".csv"))
write.csv(clk_input_ascending,paste0(PATH_output,
    "Figure_03/input_to_clk_from_ascending_",v,".csv"))
write.csv(clk_input_endocrine,paste0(PATH_output,
    "Figure_03/input_to_clk_from_endocrine_",v,".csv"))
write.csv(clk_input_NA,paste0(PATH_output,
    "Figure_03/input_to_clk_from_orphans_",v,".csv"))
#-------------------------------------------------------------------------------
# filtering for the top 10 strongest inputs to show-----------------------------
clk_input_central$hemibrain_type = ifelse(is.na(clk_input_central$hemibrain_type),
                     yes = clk_input_central$cell_type,
                     no = clk_input_central$hemibrain_type)
clk_input_central = clk_input_central[!clk_input_central$pre_pt_root_id %in% clk$clk_id,]
top_tmp = clk_input_central%>%
 group_by(hemibrain_type)%>%
 summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_input_central_top = clk_input_central[clk_input_central$hemibrain_type %in% top10_tmp,]
write.csv(clk_input_central_top,paste0(PATH_output,
     "Figure_03/Figure_03_C_input_to_clk_from_central_top_10_",v,".csv"))

clk_input_visual_proj$hemibrain_type = ifelse(is.na(clk_input_visual_proj$hemibrain_type),
                     yes = clk_input_visual_proj$cell_type,
                     no = clk_input_visual_proj$hemibrain_type)
clk_input_visual_proj = clk_input_visual_proj[!clk_input_visual_proj$pre_pt_root_id %in% clk$clk_id,]
top_tmp = clk_input_visual_proj%>%
 group_by(hemibrain_type)%>%
 summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_input_visual_proj_top = clk_input_visual_proj[clk_input_visual_proj$hemibrain_type %in% top10_tmp,]
write.csv(clk_input_visual_proj_top,paste0(PATH_output,
     "Figure_03/Figure_03_C_input_to_clk_from_visual_proj_top_10_",v,".csv"))

clk_input_visual_centri$hemibrain_type = ifelse(is.na(clk_input_visual_centri$hemibrain_type),
                       yes = clk_input_visual_centri$cell_type,
                       no = clk_input_visual_centri$hemibrain_type)
clk_input_visual_centri = clk_input_visual_centri[!clk_input_visual_centri$pre_pt_root_id %in% clk$clk_id,]
top_tmp = clk_input_visual_centri%>%
 group_by(hemibrain_type)%>%
 summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_input_visual_centri_top = clk_input_visual_centri[clk_input_visual_centri$hemibrain_type %in% top10_tmp,]
write.csv(clk_input_visual_centri_top,paste0(PATH_output,
     "Figure_03/Figure_03_C_input_to_clk_from_visual_centri_top_10_",v,".csv"))

clk_input_optic$hemibrain_type = ifelse(is.na(clk_input_optic$hemibrain_type),
                       yes = clk_input_optic$cell_type,
                       no = clk_input_optic$hemibrain_type)
clk_input_optic = clk_input_optic[!clk_input_optic$pre_pt_root_id %in% clk$clk_id,]
top_tmp = clk_input_optic%>%
 group_by(hemibrain_type)%>%
 summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_input_optic_top = clk_input_optic[clk_input_optic$hemibrain_type %in% top10_tmp,]
write.csv(clk_input_optic_top,paste0(PATH_output,
     "Figure_03/Figure_03_C_input_to_clk_from_optic_top_10_",v,".csv"))
#-------------------------------------------------------------------------------