#----Figure_03_E----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 03 Panel E for the connectivity analysis for 
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
synapses_sum_input = synapses_curated %>%
 group_by(pre_pt_root_id, post_pt_root_id) %>%
 summarise(n_synapses = length(post_pt_root_id))
classification_join = classification
colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
clk_input = left_join(synapses_sum_input[synapses_sum_input$n_synapses >= 5,],
      classification_join,by = "pre_pt_root_id")

clk_input_strong = clk_input[clk_input$n_synapses >= 80 & !clk_input$pre_pt_root_id%in%clk$clk_id,]
write.csv(clk_input_strong,paste0(PATH_output,"Figure_03/Figure_03_E_clk_input_80_v",v,".csv"))
#--------------------------------------------------------------------------------

if (!paste0("Figure_03_E_pre_clk_input_filtered_v",v,".csv") %in% input_files) {
 reticulate::source_python("./scripts/Preparations/Figure_03_E_preparation.py")
  source("./scripts/01_setup.R")
} 
pre_clk_input = read_csv(paste0(PATH_input,"tmp/Figure_03_E_pre_clk_input_filtered_v",v,".csv"),
   col_types = cols(pre_pt_root_id = col_character(),
       post_pt_root_id = col_character()))

pre_clk_input_sum = pre_clk_input%>%
 group_by(pre_pt_root_id,post_pt_root_id)%>%
 summarize(n_synapses = length(pre_pt_root_id))
pre_clk_input_sum = pre_clk_input_sum[pre_clk_input_sum$n_synapses>= 5,]

colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
pre_clk_input_sum = left_join(pre_clk_input_sum,classification_join,by = "pre_pt_root_id")
colnames(classification_join) = c("post_pt_root_id",paste(colnames(classification)[-1],"post",sep = "_"))
pre_clk_input_sum = left_join(pre_clk_input_sum,classification_join,by = "post_pt_root_id")
write.csv(pre_clk_input_sum,paste0(PATH_output,"Figure_03/Figure_03_E_pre_clk_input_sum_v",v,".csv"))

if (!paste0("Figure_03_E_pre_OCG02_input_filtered_v",v,".csv") %in% input_files) {
 reticulate::source_python("./scripts/Preparations/Figure_03_E_preparation.py")
 source("./scripts/01_setup.R")
} 
pre_OCG_input = read_csv(paste0(PATH_input,"tmp/Figure_03_E_pre_OCG02_input_filtered_v",v,".csv"),
      col_types = cols(pre_pt_root_id = col_character(),
           post_pt_root_id = col_character()))
pre_OCG_input_sum = pre_OCG_input%>%
 group_by(pre_pt_root_id,post_pt_root_id)%>%
 summarize(n_synapses = length(pre_pt_root_id))
pre_OCG_input_sum = pre_OCG_input_sum[pre_OCG_input_sum$n_synapses>= 1,]
 
classification_join = as.data.frame(classification)

colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
pre_OCG_input_sum = left_join(pre_OCG_input_sum,classification_join,by = "pre_pt_root_id")
colnames(classification_join) = c("post_pt_root_id",paste(colnames(classification)[-1],"post",sep = "_"))
pre_OCG_input_sum = left_join(pre_OCG_input_sum,classification_join,by = "post_pt_root_id")
write.csv(pre_OCG_input_sum,paste0(PATH_output,"Figure_03/Figure_03_E_pre_OCG_input_sum_v",v,".csv"))

