#----Figure_04_A_B_C and analysis as basis for A,F-H----------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 04 Panel A,B,C for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
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
  # the use of other classification files might change the results slightly 
  # depending on the changes in classification. 
  # The classification file used in the paper was updated last 26-12-2023
  stop("please go to https://codex.flywire.ai/api/download n\ 
       and dowload the classification file and save it in './input'.")
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
#Figure_04_B and analysis for D-F-----------------------------------------------
synapses_sum_output = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
colnames(classification) = c("post_pt_root_id",colnames(classification)[-1])
clk_output = left_join(synapses_sum_output[synapses_sum_output$n_synapses >= 5,],
                       classification,by = "post_pt_root_id")
clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
clk_output = left_join(clk_output,clk_join,by = "pre_pt_root_id")
write.csv(clk_output,paste0(PATH_input,"tmp/clk_output_sum_filtered_v",v,".csv"))

#output for Figure D-F----------------------------------------------------------
clk_output_strong = clk_output[clk_output$n_synapses >= 50 & !clk_output$post_pt_root_id%in%clk$clk_id,]
write.csv(clk_output_strong,paste0(PATH_output,"Figure_04/Figure_04_F-H_clk_output_50_v",v,".csv"))
#-----------------------------------------------------------------------------

SLP316 = flywire_latestid(rootid = c("720575940457052743", "720575940515316116",
          "720575940561308162","720575940605868054", "720575940609931128", 
          "720575940613608190","720575940622007846", "720575940629287724",
          "720575940631925133","720575940633984804", "720575940636136311",
          "720575940653941153"), version = v)
clk_output[clk_output$post_pt_root_id %in% SLP316,]$super_class = "central"

colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_output = left_join(clk_output,clk_join,by = "post_pt_root_id")
clk_output_total = clk_output %>%
  group_by(name_pre) %>%
  summarise(n_synapses_total = sum(n_synapses,na.rm = T),
            n_pre_partners_total = length(unique(pre_pt_root_id)))
clk_output$clk_marker = "no_clock"
clk_output[clk_output$post_pt_root_id %in% clk$clk_id,]$clk_marker = "clock"
clk_output$pre_pt_root_id_name = paste(clk_output$pre_pt_root_id,clk_output$name_pre,sep = "_")

n_post_partner = clk_output %>%
  group_by(super_class) %>%
  summarise(n_total = length(unique(post_pt_root_id)))

clk_output_grouped = clk_output %>%
  group_by(
    clk_marker,
    super_class,
    name_pre) %>%
  summarise(n_synapses_sum = sum(n_synapses,na.rm = T),
            avrg_synapses = mean(n_synapses,na.rm = T),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)))

clk_output_grouped = left_join(clk_output_grouped,clk_output_total,by = "name_pre")
clk_output_grouped$perc_of_output = clk_output_grouped$n_synapses_sum/
                                    clk_output_grouped$n_synapses_total

clk_output_grouped$name_pre = factor(clk_output_grouped$name_pre,
                                     levels = rev(c("s-LNv","l-LNv","LN_ITP",
                                                    "LNd_CRYp","LNd_CRYn","LPN",
                                                    "l-CPDN3","APDN3","s-CPDN3A",
                                                    "s-CPDN3B","s-CPDN3C",
                                                    "s-CPDN3D","s-CPDN3E","DN2",
                                                    "DN1pE","DN1pD","DN1pC",
                                                    "DN1pB","DN1pA","DN1a"))
)
clk_output_grouped$fill = paste(clk_output_grouped$super_class,
                                clk_output_grouped$clk_marker,sep = "_")
pattern_density = unique(clk_output_grouped$fill)
clk_output_grouped$super_class_corr = clk_output_grouped$super_class
clk_output_grouped$super_class_corr[clk_output_grouped$perc_of_output<0.025] = "others"
clk_output_grouped$super_class_corr[is.na(clk_output_grouped$super_class)] = "undefined"
clk_output_grouped$fill_corr = paste(clk_output_grouped$super_class_corr,
                                     clk_output_grouped$clk_marker,sep = "_")

clk_output_grouped = clk_output_grouped %>%
  group_by(name_pre,fill_corr) %>%
  summarize(perc_of_output = sum(perc_of_output))
clk_output_grouped$fill_corr = factor(clk_output_grouped$fill_corr,
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
                                                   "endocrine_no_clock"
                                                   )))

ggplot(clk_output_grouped,aes(x = "",y = perc_of_output,fill = fill_corr))+
  geom_bar_pattern(stat = "identity",col = "black",linewidth = 0.1,
                   pattern = "stripe",aes(pattern_density = fill_corr),
                   pattern_fill = "white", pattern_colour = NA,pattern_angle = 0,
                   pattern_size = 0.1, pattern_alpha = 1,
                   pattern_key_scale_factor = 0.3)+
  coord_polar("y", start = 0)+
  facet_wrap(name_pre~.,ncol = 7)+
  scale_pattern_density_manual(values =
                                 c("others_clock" = 0.1,
                                   "central_clock" = 0.1,
                                   "undefined_no_clock" = 0,
                                   "others_no_clock" = 0,
                                   "visual_projection_clock" = 0.1,
                                   "central_no_clock" = 0,
                                   "descending_no_clock" = 0,
                                   "optic_no_clock" = 0,
                                   "visual_centrifugal_no_clock" = 0,
                                   "visual_projection_no_clock" = 0,
                                   "endocrine_no_clock" = 0

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
                               "endocrine_no_clock" = "#b10071ff"))+
  ylab(label = "")+
  xlab(label = "proportion of output synapses")+
  theme(panel.background = element_rect(fill = NA,color  = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = NA),
        legend.title = element_blank(),
        text =  element_text(size = 8, colour = "black"),
        axis.text = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        panel.spacing = unit(0, "line"),
        legend.position = "bottom",
        axis.ticks = element_blank()
        
  )
ggsave(paste0(PATH_output,"Figure_04/Figure_04_B_v",v,".pdf"),
       width = 21, height =  8, units = "cm")
#-------------------------------------------------------------------------------
#analysis for Figure_04_A-------------------------------------------------------
cell_count = clk_output%>%
  group_by(super_class)%>%
  summarize(n_superclass = length(unique(post_pt_root_id)),
            n_clock = length(unique(pre_pt_root_id)))
#-------------------------------------------------------------------------------
#analysis for Figure_04_C-------------------------------------------------------
col_to_keep = c("post_pt_root_id","cell_type","class","hemibrain_type",
                "name_pre","hemisphere_pre","pre_pt_root_id","n_synapses","super_class")
clk_output = clk_output[,col_to_keep]
clk_output_visual_proj = clk_output[clk_output$super_class %in% c("visual_projection"),]
clk_output_visual_centri = clk_output[clk_output$super_class %in% c("visual_centrifugal"),]
clk_output_sensory = clk_output[clk_output$super_class %in% c("sensory"),]
clk_output_optic = clk_output[clk_output$super_class %in% c("optic"),]
clk_output_endocrine = clk_output[clk_output$super_class %in% c("endocrine"),]
clk_output_descending = clk_output[clk_output$super_class %in% c("descending"),]
clk_output_central = clk_output[clk_output$super_class %in% c("central"),]
clk_output_ascending = clk_output[clk_output$super_class %in% c("ascending"),]
clk_output_NA = clk_output[is.na(clk_output$super_class),]

write.csv(clk_output_visual_proj,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_visual_projection_v",v,".csv"))
write.csv(clk_output_visual_centri,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_visual_centrifugal_v",v,".csv"))
write.csv(clk_output_optic,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_optic_v",v,".csv"))
write.csv(clk_output_endocrine,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_endocrine_v",v,".csv"))
write.csv(clk_output_descending,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_descending_v",v,".csv"))
write.csv(clk_output_central,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_central_v",v,".csv"))
write.csv(clk_output_ascending,paste0(PATH_output,
  "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_ascending_v",v,".csv"))
# write.csv(clk_output_NA,paste0(PATH_output,
# "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_NA_v",v,".csv"))
# write.csv(clk_output_sensory,paste0(PATH_output,
# "Figure_04/Figure_04_C_clk_output_filtered_5_syn_threshold_sensory_v",v,".csv"))
#-------------------------------------------------------------------------------
# filter for top 10 output to show in Figure 04---------------------------------
clk_output_central$hemibrain_type = ifelse(is.na(clk_output_central$hemibrain_type),
                                          yes = clk_output_central$cell_type,
                                          no = clk_output_central$hemibrain_type)
clk_output_central = clk_output_central[!clk_output_central$post_pt_root_id %in% clk$clk_id,]
top_tmp = clk_output_central%>%
  group_by(hemibrain_type)%>%
  summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_output_central_top = clk_output_central[clk_output_central$hemibrain_type %in% top10_tmp,]
write.csv(clk_output_central_top,paste0(PATH_output,
  "Figure_04/Figure_04_C_output_to_clk_from_central_top_10_",v,".csv"))

clk_output_visual_proj$hemibrain_type = ifelse(is.na(clk_output_visual_proj$hemibrain_type),
                                          yes = clk_output_visual_proj$cell_type,
                                          no = clk_output_visual_proj$hemibrain_type)
clk_output_visual_proj = clk_output_visual_proj[!clk_output_visual_proj$post_pt_root_id %in% clk$clk_id,]
top_tmp = clk_output_visual_proj%>%
  group_by(hemibrain_type)%>%
  summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_output_visual_proj_top = clk_output_visual_proj[clk_output_visual_proj$hemibrain_type %in% top10_tmp,]
write.csv(clk_output_visual_proj_top,paste0(PATH_output,
  "Figure_04/Figure_04_C_output_to_clk_from_visual_proj_top_10_",v,".csv"))

clk_output_visual_centri$hemibrain_type = ifelse(is.na(clk_output_visual_centri$hemibrain_type),
                                            yes = clk_output_visual_centri$cell_type,
                                            no = clk_output_visual_centri$hemibrain_type)
clk_output_visual_centri = clk_output_visual_centri[!clk_output_visual_centri$post_pt_root_id %in% clk$clk_id,]
top_tmp = clk_output_visual_centri%>%
  group_by(hemibrain_type)%>%
  summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_output_visual_centri_top = clk_output_visual_centri[clk_output_visual_centri$hemibrain_type %in% top10_tmp,]
write.csv(clk_output_visual_centri_top,paste0(PATH_output,
  "Figure_04/Figure_04_C_output_to_clk_from_visual_centri_top_10_",v,".csv"))

clk_output_optic$hemibrain_type = ifelse(is.na(clk_output_optic$hemibrain_type),
                                        yes = clk_output_optic$cell_type,
                                        no = clk_output_optic$hemibrain_type)
clk_output_optic = clk_output_optic[!clk_output_optic$post_pt_root_id %in% clk$clk_id,]
top_tmp = clk_output_optic%>%
  group_by(hemibrain_type)%>%
  summarize(total_syn = sum(n_synapses))
top10_tmp = top_tmp[order(-top_tmp$total_syn),]$hemibrain_type[1:10]
clk_output_optic_top = clk_output_optic[clk_output_optic$hemibrain_type %in% top10_tmp,]
write.csv(clk_output_optic_top,paste0(PATH_output,
  "Figure_04/Figure_04_C_output_to_clk_from_optic_top_10_",v,".csv"))
