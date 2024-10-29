#----Figure_03_G,H -------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 03 Panel G and H for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# The first run might take a bit longer as all synapses have to be downloaded.
#-------------------------------------------------------------------------------
# needs script 'Figure_03_A-D.R' first to be run!
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
  source("./scripts/01_setup.R")
}
snyapses_curated = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
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
#updating classification--------------------------------------------------------
classification[classification$hemibrain_type%in%c("HBeyelet"),]$super_class = "sensory"
#-------------------------------------------------------------------------------
#preparation--------------------------------------------------------------------
photoreceptors = c("R1-6")
for (i in photoreceptors) {
  ids = classification[classification$cell_type == i & classification$side %in% c("right") ,]
  ids = ids[!is.na(ids$root_id),]$root_id
  write.csv(ids,paste0(PATH_input,"tmp/Figure_03_",i,"_ids_v",v,".csv"))
}
photoreceptors = c("ocellar_retinula_cell","OCG02c","R7","R8")
for (i in photoreceptors) {
  ids = classification[classification$cell_type == i,]
  ids = ids[!is.na(ids$root_id),]$root_id
  write.csv(ids,paste0(PATH_input,"tmp/Figure_03_",i,"_ids_v",v,".csv"))
}
ids = classification[classification$hemibrain_type == "HBeyelet",]
ids = ids[!is.na(ids$root_id),]$root_id
write.csv(ids,paste0(PATH_input,"tmp/Figure_03_HB_ids_v",v,".csv"))
#-------------------------------------------------------------------------------
if (0<sum(!c(
             paste0("Figure_03_GH_R1_6_output_filtered_v",v,".csv"),
             paste0("Figure_03_GH_R7_output_filtered_v",v,".csv"),
             paste0("Figure_03_GH_R8_output_filtered_v",v,".csv"),
             paste0("Figure_03_GH_OC_output_filtered_v",v,".csv"),
             paste0("Figure_03_GH_OCG_input_filtered_v",v,".csv"),
             paste0("Figure_03_GH_OCG_output_filtered_v",v,".csv"),
             paste0("Figure_03_GH_HB_output_filtered_v",v,".csv")) %in% input_files)) {
  reticulate::source_python("./scripts/Preparations/Figure_03_G_H_preparation.py")
  source("./scripts/01_setup.R")
}
synapses_R1_6 = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_R1_6_output_filtered_v",v,".csv"),
                        col_types = cols(pre_pt_root_id = col_character(),
                                         post_pt_root_id = col_character()))
synapses_R7 = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_R7_output_filtered_v",v,".csv"),
                        col_types = cols(pre_pt_root_id = col_character(),
                                         post_pt_root_id = col_character()))
synapses_R8 = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_R8_output_filtered_v",v,".csv"),
                       col_types = cols(pre_pt_root_id = col_character(),
                                        post_pt_root_id = col_character()))
synapses_OC = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_OC_output_filtered_v",v,".csv"),
                      col_types = cols(pre_pt_root_id = col_character(),
                                       post_pt_root_id = col_character()))

input_synapses_OCG = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_OCG_input_filtered_v",v,".csv"),
                      col_types = cols(pre_pt_root_id = col_character(),
                                       post_pt_root_id = col_character()))
output_synapses_OCG = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_OCG_output_filtered_v",v,".csv"),
                             col_types = cols(pre_pt_root_id = col_character(),
                                              post_pt_root_id = col_character()))

synapses_HB = read_csv(paste0(PATH_input,"tmp/Figure_03_GH_HB_output_filtered_v",v,".csv"),
                      col_types = cols(pre_pt_root_id = col_character(),
                                       post_pt_root_id = col_character()))

clk_synapses_sum_input = snyapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))

classification_join = classification
colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
clk_synapses_sum_input = left_join(clk_synapses_sum_input[clk_synapses_sum_input$n_synapses >= 5,],
                     classification_join,by = "pre_pt_root_id")

clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_synapses_sum_input = left_join(clk_synapses_sum_input,clk_join,by = "post_pt_root_id")
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
clk_synapses_sum_input = left_join(clk_synapses_sum_input,clk_join,by = "pre_pt_root_id")


classification_join = classification
colnames(classification_join) = c("post_pt_root_id",colnames(classification)[-1])
synapses_R1_6_sum = synapses_R1_6%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_R1_6_sum = left_join(synapses_R1_6_sum,
                             classification_join,by = "post_pt_root_id")

synapses_R7_sum = synapses_R7%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_R7_sum = left_join(synapses_R7_sum,
                                  classification_join,by = "post_pt_root_id")

synapses_R8_sum = synapses_R8%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_R8_sum = left_join(synapses_R8_sum,
                            classification_join,by = "post_pt_root_id")

synapses_OCG_output_sum = output_synapses_OCG%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_OCG_output_sum = left_join(synapses_OCG_output_sum,
                             classification_join,by = "post_pt_root_id")

synapses_OC_sum = synapses_OC%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_OC_sum = left_join(synapses_OC_sum,
                             classification_join,by = "post_pt_root_id")

synapses_OC_total = synapses_OC_sum%>%
  group_by(pre_pt_root_id)%>%
  summarize(n_total = sum(n_synapses))
synapses_OC_sum = left_join(synapses_OC_sum,synapses_OC_total,by = "pre_pt_root_id")
synapses_OC_sum$output_proportion = synapses_OC_sum$n_synapses/synapses_OC_sum$n_total

synapses_OCG_sum = input_synapses_OCG%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(post_pt_root_id))

classification_join = classification
colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
synapses_OCG_sum = left_join(synapses_OCG_sum,classification_join,by = "pre_pt_root_id")
 
synapses_OCG_total = synapses_OCG_sum%>%
  group_by(post_pt_root_id)%>%
  summarize(n_total = sum(n_synapses))
synapses_OCG_sum = left_join(synapses_OCG_sum,synapses_OCG_total,by = "post_pt_root_id")
synapses_OCG_sum$input_proportion = synapses_OCG_sum$n_synapses/synapses_OCG_sum$n_total

colnames(classification_join) = c("post_pt_root_id",colnames(classification)[-1])
synapses_HB_sum = synapses_HB%>%
  group_by(pre_pt_root_id,post_pt_root_id)%>%
  summarize(n_synapses = length(pre_pt_root_id))

synapses_HB_sum = left_join(synapses_HB_sum,
                           classification_join,by = "post_pt_root_id")
# Figure_03_G:------------------------------------------------------------------

bridge_neurons_R1_6_input = synapses_R1_6_sum[synapses_R1_6_sum$post_pt_root_id%in%
                                               clk_synapses_sum_input$pre_pt_root_id&
                                               synapses_R1_6_sum$n_synapses >= 5,]
bridge_neurons_R1_6_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                     bridge_neurons_R1_6_input$post_pt_root_id&
                                                     clk_synapses_sum_input$n_synapses >= 5,]

bridge_neurons_R7_input = synapses_R7_sum[synapses_R7_sum$post_pt_root_id%in%
                                               clk_synapses_sum_input$pre_pt_root_id&
                                               synapses_R7_sum$n_synapses >= 5,]                                         
bridge_neurons_R7_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                     bridge_neurons_R7_input$post_pt_root_id&
                                                     clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_R7_input,paste0(PATH_output,"Figure_03/Figure_03_G_R7_input_to_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))
write.csv(bridge_neurons_R7_output,paste0(PATH_output,"Figure_03/Figure_03_G_R7_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))


bridge_neurons_R8_input = synapses_R8_sum[synapses_R8_sum$post_pt_root_id%in%
                                            clk_synapses_sum_input$pre_pt_root_id&
                                            synapses_R8_sum$n_synapses >= 5,]                                         
bridge_neurons_R8_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                      bridge_neurons_R8_input$post_pt_root_id&
                                                      clk_synapses_sum_input$n_synapses >= 5,]


write.csv(bridge_neurons_R8_input,paste0(PATH_output,"Figure_03/Figure_03_G_R8_input_to_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))
write.csv(bridge_neurons_R8_output,paste0(PATH_output,"Figure_03/Figure_03_G_R8_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))

bridge_neurons_OC_input = synapses_OC_sum[synapses_OC_sum$post_pt_root_id%in%
                                               clk_synapses_sum_input$pre_pt_root_id&
                                               synapses_OC_sum$n_synapses >= 5,]                                         
bridge_neurons_OC_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                     bridge_neurons_OC_input$post_pt_root_id&
                                                     clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_OC_output,paste0(PATH_output,"Figure_03/Figure_03_G_OC_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))
write.csv(bridge_neurons_OC_input,paste0(PATH_output,"Figure_03/Figure_03_G_OC_input_to_bridge_neurons_for_disyn_to_clk_5_threshold",v,".csv"))

bridge_neurons_HB_input = synapses_HB_sum[synapses_HB_sum$post_pt_root_id%in%
                                           clk_synapses_sum_input$pre_pt_root_id&
                                           synapses_HB_sum$n_synapses >= 5,] 
write.csv(bridge_neurons_HB_input,paste0(PATH_output,"Figure_03/Figure_03_G_HB_bridge_neurons_input_5_threshold_v",v,".csv"))
bridge_neurons_HB_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                   bridge_neurons_HB_input$post_pt_root_id&
                                                   clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_HB_output,paste0(PATH_output,"Figure_03/Figure_03_G_HB_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))

#Figure_03_H:-------------------------------------------------------------------

bridge_neurons_R1_6_input = synapses_R1_6_sum[synapses_R1_6_sum$post_pt_root_id%in%
                                               clk_synapses_sum_input$pre_pt_root_id&
                                               synapses_R1_6_sum$n_synapses >= 3,]
bridge_neurons_R1_6_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                     bridge_neurons_R1_6_input$post_pt_root_id&
                                                     clk_synapses_sum_input$n_synapses >= 5,]

bridge_neurons_R7_input = synapses_R7_sum[synapses_R7_sum$post_pt_root_id%in%
                                            clk_synapses_sum_input$pre_pt_root_id&
                                            synapses_R7_sum$n_synapses >= 3,]                                         
bridge_neurons_R7_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                    bridge_neurons_R7_input$post_pt_root_id&
                                                    clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_R7_input,paste0(PATH_output,"Figure_03/Figure_03_H_R7_input_to_bridge_neurons_for_disyn_to_clk_3_threshold_v",v,".csv"))
write.csv(bridge_neurons_R7_output,paste0(PATH_output,"Figure_03/Figure_03_H_R7_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))


bridge_neurons_R8_input = synapses_R8_sum[synapses_R8_sum$post_pt_root_id%in%
                                            clk_synapses_sum_input$pre_pt_root_id&
                                            synapses_R8_sum$n_synapses >= 3,]                                         
bridge_neurons_R8_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                    bridge_neurons_R8_input$post_pt_root_id&
                                                    clk_synapses_sum_input$n_synapses >= 5,]


write.csv(bridge_neurons_R8_input,paste0(PATH_output,"Figure_03/Figure_03_H_R8_input_to_bridge_neurons_for_disyn_to_clk_3_threshold_v",v,".csv"))
write.csv(bridge_neurons_R8_output,paste0(PATH_output,"Figure_03/Figure_03_H_R8_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))


bridge_neurons_OC_input = synapses_OC_sum[synapses_OC_sum$post_pt_root_id%in%
                                            clk_synapses_sum_input$pre_pt_root_id&
                                            synapses_OC_sum$n_synapses >= 3,]                                         
bridge_neurons_OC_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                    bridge_neurons_OC_input$post_pt_root_id&
                                                    clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_OC_input,paste0(PATH_output,"Figure_03/Figure_03_H_OC_input_to_bridge_neurons_for_disyn_to_clk_3_threshold_v",v,".csv"))
write.csv(bridge_neurons_OC_output,paste0(PATH_output,"Figure_03/Figure_03_H_OC_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))


bridge_neurons_HB_input = synapses_HB_sum[synapses_HB_sum$post_pt_root_id%in%
                                            clk_synapses_sum_input$pre_pt_root_id&
                                            synapses_HB_sum$n_synapses >= 3,]
write.csv(bridge_neurons_HB_input,paste0(PATH_output,"Figure_03/Figure_03_H_HB_bridge_neurons_input_3_threshold_v",v,".csv"))
bridge_neurons_HB_output = clk_synapses_sum_input[clk_synapses_sum_input$pre_pt_root_id%in%
                                                    bridge_neurons_HB_input$post_pt_root_id&
                                                    clk_synapses_sum_input$n_synapses >= 5,]
write.csv(bridge_neurons_HB_output,paste0(PATH_output,"Figure_03/Figure_03_H_HB_bridge_neurons_for_disyn_to_clk_5_threshold_v",v,".csv"))
