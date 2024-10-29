#----Figure_04E ----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 04 Panel E for the connectivity analysis 
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
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_output_sum = left_join(clk_output_sum,clk_join,by = "post_pt_root_id")

clk_output_sum = clk_output_sum[clk_output_sum$pre_pt_root_id %in%
                                  clk[clk$clk_name %in% c("DN1pA",
                                                          "s-CPDN3C",
                                                          "s-CPDN3D"),]$clk_id,]
clk_output_sum = clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                  clk[clk$clk_name %in% c("DN1pC","DN1pD","DN1pE",
                                                          "s-CPDN3C","s-CPDN3D"),]$clk_id,]
clk_output_sum$pre_post = paste0(clk_output_sum$name_pre,clk_output_sum$name_post)
clk_output_sum = clk_output_sum[!clk_output_sum$pre_post%in%c("s-CPDN3Cs-CPDN3D",
                                                              "s-CPDN3Ds-CPDN3C",
                                                              "s-CPDN3Ds-CPDN3D",
                                                              "s-CPDN3Cs-CPDN3C"),]
clk_output_sum$pre_post = NULL
clk_output_sum$connection_status_pre = NA
clk_output_sum$connection_status_post = NA
clk_output_sum$connection_status_pre = ifelse(clk_output_sum$pre_pt_root_id %in%
                                            clk_output_sum$post_pt_root_id,
                                          yes= "pre_post",no = "pre")
clk_output_sum$connection_status_post = ifelse(clk_output_sum$post_pt_root_id %in%
                                                clk_output_sum$pre_pt_root_id,
                                              yes= "pre_post",no = "post")
clk_output_sum[clk_output_sum$name_pre %in% c("s-CPDN3C","s-CPDN3D"),]$name_pre = "s-CPDN3C_D"
clk_output_sum[clk_output_sum$name_post %in% c("s-CPDN3C","s-CPDN3D"),]$name_post = "s-CPDN3C_D"

clk_output_sum[clk_output_sum$name_post %in% c("DN1pC","DN1pD","DN1pE"),]$name_post = "DN1pC-E"

clk_output_grouped = clk_output_sum %>%
  group_by(name_pre,name_post,hemisphere_pre,hemisphere_post,
           connection_status_pre,connection_status_post)%>%
  summarize(n_syn_total = sum(n_synapses))
cell_count_pre = clk_output_sum %>%
  group_by(name_pre,hemisphere_pre,
           connection_status_pre)%>%
  summarize(n_pre = length(unique(pre_pt_root_id)))
colnames(cell_count_pre) = c("name","hemisphere","connection_status","n")
cell_count_post = clk_output_sum %>%
  group_by(name_post,hemisphere_post,
           connection_status_post)%>%
  summarize(n_pre = length(unique(post_pt_root_id)))
colnames(cell_count_post) = c("name","hemisphere","connection_status","n")

cell_count = rbind(cell_count_pre,cell_count_post[cell_count_post$connection_status == "post",])

write.csv(clk_output_grouped,paste0(PATH_output,"Figure_04/Figure_04_E_connections.csv"))
write.csv(cell_count,paste0(PATH_output,"Figure_04/Figure_04_E_cell_count.csv"))
