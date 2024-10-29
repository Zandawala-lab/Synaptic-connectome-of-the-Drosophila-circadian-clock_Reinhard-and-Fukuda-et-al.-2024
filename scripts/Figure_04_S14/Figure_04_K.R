#----Figure_04 K------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 04 Panel K for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# The first run might take  long as all synapses have to be downloaded for
# disynaptic connections.
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
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

NSC = read_delim(paste0(PATH_input,"NSC_v",v,".csv"),
                 col_types = cols(NSC_id = col_character()),delim = ",")

if (!paste0("NSC_input_from_clk_bridge_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_K_preparations.py")
  source("./scripts/01_setup.R")
}
NSC_input = read_csv(paste0(PATH_input,"tmp/Figure_04_K_NSC_input_from_clk_bridge_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_root_id = col_character(),
                                             post_pt_root_id = col_character()))

if (!paste0("clk_neuron_output_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_K_preparations.py")
  source("./scripts/01_setup.R")
}
clk_output = read_csv(paste0(PATH_input,"tmp/clk_neuron_output_filtered_v",v,".csv"), 
                      col_types = cols(pre_pt_root_id = col_character(), 
                                       post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
NSC_input_sum = NSC_input %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))

NSC_input_sum = NSC_input_sum[NSC_input_sum$n_synapses>=5,]

clk_output_sum = clk_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))
clk_output_sum = clk_output_sum[clk_output_sum$n_synapses>=5,]

clk_to_NSC_direct_con = clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                         NSC$NSC_id,]
clk_join = clk
colnames(clk_join) = c("clk_name","pre_pt_root_id", "hemisphere_clk")

clk_to_NSC_direct_con = left_join(clk_to_NSC_direct_con, clk_join,
                                  by = "pre_pt_root_id")

NSC_join = NSC
colnames(NSC_join) = c("nc_name","post_pt_root_id","hemisphere")

clk_to_NSC_direct_con = left_join(clk_to_NSC_direct_con,NSC_join,
                                  by = "post_pt_root_id")
write.csv(clk_to_NSC_direct_con,paste0(PATH_output,
                        "Figure_04/Figure_04_K_clk_to_NSC_direct_v",v,".csv"))

print(paste("number of clk neurons contacting directly NSCs:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           NSC$NSC_id,]$pre_pt_root_id))))
print(paste("number of NSCs reached directly by clk neurons:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           NSC$NSC_id,]$post_pt_root_id))))

print(paste("number of clk neurons connected disynaptic to NSC:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                               NSC_input_sum$pre_pt_root_id,]$pre_pt_root_id))))
print(paste("number of bridge neurons to NSC:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                               NSC_input_sum$pre_pt_root_id,]$post_pt_root_id))))
print(paste("number of NSCs reached by bridge neurons:",
            length(unique(NSC_input_sum[NSC_input_sum$pre_pt_root_id %in%
                               clk_output_sum$post_pt_root_id,]$post_pt_root_id))))

colnames(NSC_input_sum)[colnames(NSC_input_sum) == "pre_pt_root_id"] = "bridge_pt_root_id"
colnames(NSC_input_sum)[colnames(NSC_input_sum) == "post_pt_root_id"] = "NSC_post_pt_root_id"
colnames(NSC_input_sum)[colnames(NSC_input_sum) == "n_synapses"] = "n_synapses_to_NSC"

colnames(clk_output_sum)[colnames(clk_output_sum) == "pre_pt_root_id"] = "clk_pre_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "post_pt_root_id"] = "bridge_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "n_synapses"] = "n_synapses_from_clk"

classification_join = classification
classification_join = classification_join[classification_join$root_id %in%
                                            NSC_input_sum$bridge_pt_root_id,]
colnames(classification_join) = c("bridge_pt_root_id",
                                  paste0(colnames(classification)[-1],"_bridge"))

connection_clk_to_NSC = left_join(classification_join,
                                  NSC_input_sum,by = "bridge_pt_root_id")

classification_join = classification
classification_join = classification_join[classification_join$root_id %in%
                                            NSC_input_sum$NSC_post_pt_root_id,]
colnames(classification_join) = c("NSC_post_pt_root_id",
                                  paste0(colnames(classification)[-1],"_NSC"))

connection_clk_to_NSC = left_join(connection_clk_to_NSC,
                                  classification_join,by = "NSC_post_pt_root_id")

clk_join = clk[clk$clk_id %in% clk_output_sum[clk_output_sum$bridge_pt_root_id %in%
                          NSC_input_sum$bridge_pt_root_id,]$clk_pre_pt_root_id,]

colnames(clk_join) = c("clk_name","clk_pre_pt_root_id","hemisphere_clk")

clk_output_sum = left_join(clk_output_sum,clk_join,by="clk_pre_pt_root_id")

connection_clk_to_NSC = left_join(clk_output_sum,
                                  connection_clk_to_NSC,by="bridge_pt_root_id",
                                  relationship = "many-to-many")

connection_clk_to_NSC = 
  connection_clk_to_NSC[!is.na(connection_clk_to_NSC$NSC_post_pt_root_id),]

NSC_join = NSC[,-3]
colnames(NSC_join) = c("NSC_name", "NSC_post_pt_root_id")
connection_clk_to_NSC = left_join(connection_clk_to_NSC, NSC_join,
                                  by = "NSC_post_pt_root_id")

write.csv(connection_clk_to_NSC,
          paste0(PATH_output,"Figure_04/Figure_04_K_clk_to_NSC_disynaptic_v",v,".csv"),
          row.names = FALSE)

connection_clk_to_NSC_grouped_NSC = connection_clk_to_NSC %>%
  group_by(NSC_name) %>%
  summarize(n = length(unique(NSC_post_pt_root_id)))
colnames(connection_clk_to_NSC_grouped_NSC) = c("NSC_name", "n")

connection_clk_to_NSC_grouped_NSC$NSC_name = 
  factor(connection_clk_to_NSC_grouped_NSC$NSC_name,
         levels = c("l-NSC","m-NSC","SEZ-NSC"))

ggplot(connection_clk_to_NSC_grouped_NSC,aes(x="",y=n,fill=NSC_name,group=NSC_name))+
  geom_bar(stat="identity",col="white",linewidth=0.1)+
  coord_polar("y", start=0)+
  scale_x_discrete(expand=c(0,0))+
  ylab(label="")+
  xlab(label="proportion of output synapses")+
  theme(panel.background=element_rect(fill=NA,color =NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "black", fill = NA,linewidth=NA),
        legend.title = element_blank(),
        text= element_text(size = 8, colour = "black"),
        axis.text = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        panel.spacing = unit(0, "line"),
        legend.position = "bottom",
        axis.ticks = element_blank()
        
  )
ggsave(paste0(PATH_output,"Figure_04/Figure_04_K_NSC_v",v,".pdf"),
       width = 21, height= 8, units = "cm")
#-------------------------------------------------------------------------------