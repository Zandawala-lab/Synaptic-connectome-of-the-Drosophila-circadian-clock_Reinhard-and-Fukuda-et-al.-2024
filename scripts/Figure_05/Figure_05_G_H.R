#----Figure_05G_H---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 05 Panel G,H for the connectivity analysis for 
# Reinhard and Fukuda et al. 
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
clk  =  read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                   col_types  =  cols(clk_id  =  col_character()),delim  =  ",")

if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
}

if (!paste0("classification_v",v,".csv") %in% input_files) {
  stop("please go to https://codex.flywire.ai/api/download n\ 
       and dowload the classification file for the current version and save it in './input'.")
}

classification  =  read_delim(paste0(PATH_input,"classification_v",v,".csv"), 
                              delim  =  ",",
                              escape_double  =  FALSE,
                              col_types  =  cols(root_id  =  col_character(),
                                                 flow  =  col_character()),
                              trim_ws  =  TRUE)
classification = as.data.frame(classification)

#-------------------------------------------------------------------------------
#Figure_05_G-------------------------------------------------------------------
clk_to_DN = read_delim(paste0(PATH_output,"Figure_05/Figure_05_A_clk_to_DN_direct_v783.csv"),
              col_types  =  cols(pre_pt_root_id  =  col_character(),
                              post_pt_root_id  =  col_character()),delim  =  ",")[,-1]
clk_to_DN_grouped = clk_to_DN%>%
  group_by(clk_name,cell_type,hemibrain_type,hemisphere_clk,side)%>%
  summarize(n_syn_total = sum(n_synapses))
clk_to_DN_grouped$cell_type = ifelse(is.na(clk_to_DN_grouped$cell_type),yes=clk_to_DN_grouped$hemibrain_type,no = clk_to_DN_grouped$cell_type)
clk_to_DN_grouped$hemibrain_type=NULL
clk_to_DN_grouped$both_hemispheres = FALSE
clk_to_DN_grouped[clk_to_DN_grouped$cell_type%in%unique(clk_to_DN_grouped[clk_to_DN_grouped$side=="left",]$cell_type)& 
                  clk_to_DN_grouped$cell_type%in%unique(clk_to_DN_grouped[clk_to_DN_grouped$side=="right",]$cell_type),]$both_hemispheres = TRUE
clk_to_DN_grouped$both_hemispheres_clk = FALSE                  
clk_to_DN_grouped[clk_to_DN_grouped$clk_name%in%unique(clk_to_DN_grouped[clk_to_DN_grouped$hemisphere_clk=="left",]$clk_name)& 
                    clk_to_DN_grouped$clk_name%in%unique(clk_to_DN_grouped[clk_to_DN_grouped$hemisphere_clk=="right",]$clk_name),]$both_hemispheres_clk = TRUE


clk_to_DN_grouped = clk_to_DN_grouped%>%
  group_by(clk_name,cell_type,both_hemispheres,both_hemispheres_clk)%>%
  summarize(n_syn_total = sum(n_syn_total))

clk_to_DN_grouped$clk_name = ifelse(clk_to_DN_grouped$both_hemispheres_clk==F,yes="others",no=clk_to_DN_grouped$clk_name)
clk_to_DN_grouped$cell_type = ifelse(clk_to_DN_grouped$both_hemispheres==F,yes="others",no=clk_to_DN_grouped$cell_type)


clk_to_DN_grouped = clk_to_DN_grouped%>%
  group_by(clk_name,cell_type)%>%
  summarize(n_syn_total = sum(n_syn_total))

DN_top3 = unique(clk_to_DN[clk_to_DN$cell_type %in% c("DNpe048","DNpe033","DNpe041"),]$post_pt_root_id)
write.csv(DN_top3,paste0(PATH_input,"/tmp/Figure_05_DN_top3.csv"))
#-------------------------------------------------------------------------------
#Figure_05_H--------------------------------------------------------------------
synapses_curated = read_csv(paste0(PATH_input,"/tmp/Figure_05_DN_top3_neuron_input_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_supervoxel_id = col_character(),
                                             pre_pt_root_id = col_character(),
                                             post_pt_supervoxel_id = col_character(),
                                             post_pt_root_id = col_character()))
synapses_sum = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
synapses_sum = synapses_sum[synapses_sum$n_synapses >= 5,]

clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")

synapses_sum = left_join(synapses_sum,clk_join,by = "pre_pt_root_id")

classification_join = classification
colnames(classification_join) = c("post_pt_root_id",colnames(classification)[-1])

synapses_sum = left_join(synapses_sum,classification_join,by = "post_pt_root_id")
synapses_sum$clk = NA
synapses_sum$clk = ifelse(is.na(synapses_sum$name_pre),yes="no_clk",no="clk")

synapses_sum_grouped = synapses_sum%>%
  group_by(clk,cell_type)%>%
  summarize(n_syn_total = sum(n_synapses))

syn_total = synapses_sum%>%
  group_by(cell_type)%>%
  summarize(n_syn_total_DN = sum(n_synapses))
synapses_sum_grouped = left_join(synapses_sum_grouped,syn_total, by="cell_type")
synapses_sum_grouped$perc = synapses_sum_grouped$n_syn_total/synapses_sum_grouped$n_syn_total_DN

plot = ggplot(synapses_sum_grouped,aes(
   x = "",y = perc, fill=clk))+
  geom_col()+
  facet_wrap(cell_type~.,ncol = 7)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0(PATH_output,"Figure_05/Figure_05_H.pdf"),width=10,height=15,units = "cm")
source_data =plot$data
write.csv(source_data,"./output/Figure_05/Figure_05_H_source-Data.csv")
#-------------------------------------------------------------------------------
