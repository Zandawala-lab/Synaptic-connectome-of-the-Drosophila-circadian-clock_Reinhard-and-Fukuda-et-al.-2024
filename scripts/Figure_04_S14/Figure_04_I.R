#----Figure_04_I----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 04 Panel I for the connectivity analysis for 
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
CX_ids<-na.omit(classification[classification$class == "CX",]$root_id)
write.csv(CX_ids,paste0(PATH_input,"tmp/Figure_04_CX_ids_v",v,".csv"),row.names = FALSE)

if (!paste0("Figure_04_I_CX_input_from_clk_bridge_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_I_preparations.py")
  source("./scripts/01_setup.R")
}

CX_input = read_csv(paste0(PATH_input,
                           "tmp/Figure_04_I_CX_input_from_clk_bridge_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_root_id = col_character(),
                                             post_pt_root_id = col_character()))

CX_input_sum = CX_input %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))
CX_input_sum = CX_input_sum[CX_input_sum$n_synapses>=5,]

clk_output_sum = clk_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))
clk_output_sum = clk_output_sum[clk_output_sum$n_synapses>=5,]
# calculating direct connections------------------------------------------------
clk_to_CX_direct_con = clk_output_sum[clk_output_sum$post_pt_root_id %in% CX_ids,]

clk_join = clk
colnames(clk_join) = c("clk_name","pre_pt_root_id", "hemisphere_clk")

clk_to_CX_direct_con = left_join(clk_to_CX_direct_con, clk_join, by = "pre_pt_root_id")

classification_join = classification
colnames(classification_join) = c("post_pt_root_id",colnames(classification_join)[-1])

clk_to_CX_direct_con = left_join(clk_to_CX_direct_con,classification_join,
                                 by = "post_pt_root_id")
write.csv(clk_to_CX_direct_con,
          paste0(PATH_output,"Figure_04/Figure_04_I_clk_to_CX_direct_v",v,".csv"))

print(paste("number of clk neurons contacting directly CX neurons:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           CX_ids,]$pre_pt_root_id))))
print(paste("number of CX neurons reached directly by clk neurons:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           CX_ids,]$post_pt_root_id))))
#-------------------------------------------------------------------------------
#disynaptic connections:--------------------------------------------------------
print(paste("number of clk neurons connected disynaptic to CX:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           CX_input_sum$pre_pt_root_id,]$pre_pt_root_id))))
print(paste("number of bridge neurons to CX:",
            length(unique(clk_output_sum[clk_output_sum$post_pt_root_id %in%
                                           CX_input_sum$pre_pt_root_id,]$post_pt_root_id))))
print(paste("number of CX neurons reached by bridge neurons:",
            length(unique(CX_input_sum[CX_input_sum$pre_pt_root_id %in%
                                         clk_output_sum$post_pt_root_id,]$post_pt_root_id))))

colnames(CX_input_sum)[colnames(CX_input_sum) == "pre_pt_root_id"] = "bridge_pt_root_id"
colnames(CX_input_sum)[colnames(CX_input_sum) == "post_pt_root_id"] = "CX_post_pt_root_id"
colnames(CX_input_sum)[colnames(CX_input_sum) == "n_synapses"] = "n_synapses_to_CX"

colnames(clk_output_sum)[colnames(clk_output_sum) == "pre_pt_root_id"] = "clk_pre_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "post_pt_root_id"] = "bridge_pt_root_id"
colnames(clk_output_sum)[colnames(clk_output_sum) == "n_synapses"] = "n_synapses_from_clk"

classification_join = classification
classification_join = classification_join[classification_join$root_id %in% 
                                            CX_input_sum$bridge_pt_root_id,]
colnames(classification_join) = c("bridge_pt_root_id",
                                  paste0(colnames(classification)[-1],"_bridge"))

connection_clk_to_cx = left_join(classification_join,CX_input_sum,
                                 by = "bridge_pt_root_id")

classification_join = classification
classification_join = classification_join[classification_join$root_id %in%
                                            CX_input_sum$CX_post_pt_root_id,]
colnames(classification_join) = c("CX_post_pt_root_id",
                                  paste0(colnames(classification)[-1],"_CX"))

connection_clk_to_cx = left_join(connection_clk_to_cx,
                                 classification_join,by = "CX_post_pt_root_id")

clk_join = clk[clk$clk_id %in% clk_output_sum[clk_output_sum$bridge_pt_root_id %in%
                               CX_input_sum$bridge_pt_root_id,]$clk_pre_pt_root_id,]

colnames(clk_join) = c("clk_name","clk_pre_pt_root_id","hemisphere_clk")

clk_output_sum = left_join(clk_output_sum,clk_join,by="clk_pre_pt_root_id")

connection_clk_to_cx = left_join(clk_output_sum,
                                 connection_clk_to_cx,by="bridge_pt_root_id",
                                 relationship = "many-to-many")

connection_clk_to_cx = connection_clk_to_cx[!is.na(connection_clk_to_cx$CX_post_pt_root_id),]
#-------------------------------------------------------------------------------
#plotting:----------------------------------------------------------------------
connection_clk_to_cx$CX_type = NA
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
  c("FB7A,FB7K","FB5Y",	"FB6F", "FB5I",	"FB6H","FB7B",	
    "FB6S,FB6T,FB7C,FB7F,FB7M,FB8D,FB8F","FB7L", "hDeltaD",	
    "vDeltaA,vDeltaB,vDeltaC,vDeltaD,vDeltaE,vDeltaJ,vDeltaK,vDeltaL,vDeltaM",
    "hDeltaE",	"FB8A,FB8H",	"hDeltaL","FB8B",	"FB6K",	"FS1A,FS1B",	"FB4N",
    "FB1H", "FB5C,FB5D,FB5E,FB5F,FB5G,FB5O,FB5P,FB5Q,FB5X,FB5Z","FB4Y", "FB1G",
    "FB8C","FB5J,FB5K,FB5L,FB5M,FB5N,FB5W",	"FB5H","FB6Z",	"hDeltaF", "FB5A",
    "hDeltaK", "FB2A",	"FB6M",	"FB2B",	"hDeltaC","FB6V",	"FB1C",	"FC2C",
    "FB4O,FB4P,FB4Q,FB4R","FC2B",	"FB4X",	"FB4A,FB4D,FB4E,FB4F,FB4G,FB4H,FB4I,FB4J",
    "FB4K",	"FB4M",	"FB5AB", "FS3",	"FB9A,FB9C",	"FS4A,FS4B,FS4C", "FB5V",
    "FB6O",	"FB4B",	"FB4C","FB8I","FB8G","FB8E", "FB8D","FB7J", "FB7G","FB7F,FB8F",
    "FB7E","FB6W", "FB6T","FB6R","FB6E","FB5Z","FB5X", "FB5O","FB5C,FB5F,FB5G,FB5P,FB5Q",
    "FB4P_a","FB4P,FB4Q,FB4R","FB4O"),]$CX_type = "FB"
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
     c("ExR1",	"ExR2", "ER4m",	"ExR3", "ExR5",	"ER3w", "ER3d",	"EL",	"ER3m",
       "ER2", "ER4d",	"ER3p",	"ExR7",	"ER3a", "ExR6",	"ExR4",	"ER5",	"ER1"),
     ]$CX_type = "EB"
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
                       c("SA1", "SA2", "SA3"),]$CX_type = "others" #asymetric_body
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
                       c("LNO1", "LCNOpm", "LNOa", "LNO2"),]$CX_type = "NO"
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
                       c("EPG"),]$CX_type = "EPG"
connection_clk_to_cx[connection_clk_to_cx$hemibrain_type_CX %in% 
                       c("PFGs"),]$CX_type = "PFGs"
connection_clk_to_cx[is.na(connection_clk_to_cx$hemibrain_type_CX),]$CX_type = "others"

write.csv(connection_clk_to_cx,
          paste0(PATH_output,"Figure_04/Figure_04_I_clk_to_CX_disynaptic_v",v,".csv"),
          row.names = FALSE)

connection_clk_to_cx_grouped = connection_clk_to_cx%>%
  group_by(CX_type)%>%
  summarize(n = length(unique(CX_post_pt_root_id)))

connection_clk_to_cx_grouped$CX_type = factor(
  connection_clk_to_cx_grouped$CX_type,
  levels = c("PFGs","EPG","NO","asymetric_body", "others","EB","FB"))

ggplot(connection_clk_to_cx_grouped,aes(x = "",y = n,fill = CX_type))+
  geom_bar(stat = "identity",col = "black",linewidth = 0.1)+
  coord_polar("y", start = 0)+
  scale_x_discrete(expand = c(0,0))+
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
ggsave(paste0(PATH_output,"Figure_04/Figure_04_I_v",v,".pdf"),
       width = 21, height =  8, units = "cm")
#-------------------------------------------------------------------------------