#----Figure_05_A----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 05 Panel A for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# The first run might take  long as all synapses have to be downloaded for
# disynaptic connections.
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_output_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_preparations.py")
  source("./scripts/01_setup.R")
}
clk_output = read_csv(paste0(PATH_input,"tmp/clk_neuron_output_filtered_v",v,".csv"), 
                      col_types = cols(pre_pt_root_id = col_character(), 
                                       post_pt_root_id = col_character()))

if (!paste0("classification_v",v,".csv") %in% input_files) {
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
DN_ids<-na.omit(classification[classification$super_class == "descending",]$root_id)
write.csv(DN_ids,paste0(PATH_input,"./tmp/Figure_05_DN_ids_v",v,".csv"),row.names = FALSE)

if (!paste0("Figure_05_A_DN_input_from_clk_bridge_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_05_A_preparations.py")
  source("./scripts/01_setup.R")
}

DN_input = read_csv(paste0(PATH_input,"tmp/Figure_05_A_DN_input_from_clk_bridge_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_root_id = col_character(),
                                             post_pt_root_id = col_character()))

DN_input_sum = DN_input %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))

DN_input_sum = DN_input_sum[DN_input_sum$n_synapses>=5,]

clk_output_sum = clk_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))
clk_output_sum = clk_output_sum[clk_output_sum$n_synapses>=5,]

clk_to_DN_direct_con = clk_output_sum[clk_output_sum$post_pt_root_id %in% DN_ids,]
clk_join = clk
colnames(clk_join) = c("clk_name","pre_pt_root_id","hemisphere_clk")

clk_to_DN_direct_con = left_join(clk_to_DN_direct_con, clk_join,
                                 by = "pre_pt_root_id")

classification_join = classification
colnames(classification_join) = c("post_pt_root_id",colnames(classification_join)[-1])

clk_to_DN_direct_con = left_join(clk_to_DN_direct_con,classification_join,
                                 by = "post_pt_root_id")
write.csv(clk_to_DN_direct_con,
          paste0(PATH_output,"Figure_05/Figure_05_A_clk_to_DN_direct_v",v,".csv"))


print(paste("number of clk neurons contacting directly DNs:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           DN_ids,]$pre_pt_root_id))))
print(paste("number of DNs reached directly by clk neurons:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           DN_ids,]$post_pt_root_id))))

print(paste("number of clk neurons connected disynaptic to DN:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           DN_input_sum$pre_pt_root_id,]$pre_pt_root_id))))
print(paste("number of bridge neurons to DN:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           DN_input_sum$pre_pt_root_id,]$post_pt_root_id))))
print(paste("number of DNs reached by bridge neurons:",
            length(unique(DN_input_sum[DN_input_sum$pre_pt_root_id %in%
                                         clk_output_sum$post_pt_root_id,]$post_pt_root_id))))

colnames(DN_input_sum)[colnames(DN_input_sum) == "pre_pt_root_id"] = "bridge_pt_root_id"
colnames(DN_input_sum)[colnames(DN_input_sum) == "post_pt_root_id"] = "DN_post_pt_root_id"
colnames(DN_input_sum)[colnames(DN_input_sum) == "n_synapses"] = "n_synapses_to_DN"

colnames(clk_output_sum)[colnames(clk_output_sum) == "pre_pt_root_id"] = "clk_pre_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "post_pt_root_id"] = "bridge_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "n_synapses"] = "n_synapses_from_clk"

classification_join = classification
classification_join = classification_join[classification_join$root_id %in%
                                            DN_input_sum$bridge_pt_root_id,]
colnames(classification_join) = c("bridge_pt_root_id",
                                  paste0(colnames(classification)[-1],"_bridge"))

connection_clk_to_DN = left_join(classification_join,DN_input_sum,
                                 by = "bridge_pt_root_id")

classification_join = classification
classification_join = classification_join[classification_join$root_id %in%
                                            DN_input_sum$DN_post_pt_root_id,]
colnames(classification_join) = c("DN_post_pt_root_id",
                                  paste0(colnames(classification)[-1],"_DN"))

connection_clk_to_DN = left_join(connection_clk_to_DN,classification_join,
                                 by = "DN_post_pt_root_id")

clk_join = clk[clk$clk_id %in% 
           clk_output_sum[clk_output_sum$bridge_pt_root_id %in%
           DN_input_sum$bridge_pt_root_id,]$clk_pre_pt_root_id,]

colnames(clk_join) = c("clk_name","clk_pre_pt_root_id","hemisphere_clk")

clk_output_sum = left_join(clk_output_sum,clk_join,by="clk_pre_pt_root_id")

connection_clk_to_DN = left_join(clk_output_sum,connection_clk_to_DN,
                           by="bridge_pt_root_id",relationship = "many-to-many")

connection_clk_to_DN = connection_clk_to_DN[!is.na(
  connection_clk_to_DN$DN_post_pt_root_id),]

write.csv(connection_clk_to_DN,
          paste0(PATH_output,"Figure_05/Figure_05_A_clk_to_DN_disynaptic_v",v,".csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
