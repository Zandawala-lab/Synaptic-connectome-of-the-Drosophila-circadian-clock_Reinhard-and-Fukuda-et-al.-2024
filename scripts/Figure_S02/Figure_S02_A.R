#----Figure_S02_A --------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S02 Panel A for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# set variables:----------------------------------------------------------------
con_threshold = 1
hemisphere = "right"

if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
}
clk_input = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
                     col_types = cols(pre_pt_supervoxel_id  = col_character(),
                                      pre_pt_root_id = col_character(),
                                      post_pt_supervoxel_id = col_character(),
                                      post_pt_root_id = col_character()))

if (!paste0("clk_neuron_output_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_preparations.py")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")

clk_output = read_csv(paste0(PATH_input,"/tmp/clk_neuron_output_filtered_v",v,".csv"),
                      col_types = cols(pre_pt_supervoxel_id = col_character(),
                                       pre_pt_root_id = col_character(),
                                       post_pt_supervoxel_id = col_character(),
                                       post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
clk_input_sum = clk_input %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))

clk_join = clk
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_input_sum = left_join(clk_input_sum,clk_join,by = "post_pt_root_id")
clk_input_sum$post_name_hemisphere_id = paste(
  clk_input_sum$hemisphere_post,
  clk_input_sum$name_post,
  clk_input_sum$post_pt_root_id,sep="_"
)

clk_input_sum_heat = clk_input_sum[
  clk_input_sum$hemisphere_post==hemisphere&
    clk_input_sum$n_synapses>=con_threshold, c("n_synapses", "pre_pt_root_id",
                                               "post_name_hemisphere_id")]
clk_input_sum_heat_w = spread(clk_input_sum_heat,                                  
                              key =post_name_hemisphere_id ,
                              value = as.numeric(n_synapses))
row_names = clk_input_sum_heat_w$pre_pt_root_id
clk_input_sum_heat_w = as.matrix(clk_input_sum_heat_w[,-1])
row.names(clk_input_sum_heat_w) = row_names
clk_input_sum_heat_w[is.na(clk_input_sum_heat_w)] = 0
input_similarity = cosine_sim(clk_input_sum_heat_w,transpose=F)

clk_output_sum = clk_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(pre_pt_root_id))

colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
clk_output_sum = left_join(clk_output_sum,clk_join,by = "pre_pt_root_id")
clk_output_sum$pre_name_hemisphere_id = paste(
  clk_output_sum$hemisphere_pre,
  clk_output_sum$name_pre,
  clk_output_sum$pre_pt_root_id,sep="_")

clk_output_sum_heat = clk_output_sum[
  clk_output_sum$hemisphere_pre==hemisphere&
    clk_output_sum$n_synapses>=con_threshold, c("n_synapses", "post_pt_root_id",
                                                "pre_name_hemisphere_id")]
clk_output_sum_heat_w = spread(clk_output_sum_heat,                                  
                               key =pre_name_hemisphere_id ,
                               value = as.numeric(n_synapses))
row_names = clk_output_sum_heat_w$post_pt_root_id
clk_output_sum_heat_w = as.matrix(clk_output_sum_heat_w[,-1])
row.names(clk_output_sum_heat_w) = row_names
clk_output_sum_heat_w[is.na(clk_output_sum_heat_w)] = 0
output_similarity = cosine_sim(clk_output_sum_heat_w,transpose=F)

similarity_matrix = (input_similarity + output_similarity)/2
hc_similarity = nhclust(scoremat = (similarity_matrix))

d_hc = colour_clusters(hc_similarity, k=17)
con_clust = as.data.frame(cutree(hc_similarity,k =20)) #defining clusters

pdf(file = paste0(PATH_output,"Figure_S02/Figure_S02_A_",hemisphere,
                  "_threshold_",con_threshold,"_v",v,".pdf"),width = 20,height = 20)

heatmap(similarity_matrix,
        Rowv = d_hc, Colv = d_hc,
        scale = "none",
        margins = c(20,20),
        col=grey_scale(length(clk$clk_id))
)
dev.off() 


