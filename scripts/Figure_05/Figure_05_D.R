#----Figure_05------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 05 Panel D for the connectivity analysis for 
# Reinhard and Fukuda et al. 
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
clk  =  read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                   col_types  =  cols(clk_id  =  col_character()),delim  =  ",")

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

#nt_neurons_v_2024-02-26_v783_codex
nt = read_csv(paste0(PATH_input,"nt_neurons_v783.csv"),
              col_types  =   cols(root_id  =  col_character()))
#-------------------------------------------------------------------------------
clk_join = clk
colnames(clk_join) = c("clk_name","root_id","hemisphere_clk")
clk_nt = left_join(clk_join,nt, by = "root_id")
clk_nt$group=NULL
classification_join = classification
clk_nt = left_join(clk_nt,classification_join[,c("root_id","hemilineage",
                                                 "cell_type","hemibrain_type")],
                   by = "root_id")
clk_nt$cell_type = ifelse(is.na(clk_nt$cell_type),yes = clk_nt$hemibrain_type,
                          no =clk_nt$cell_type )
clk_nt$hemibrain_type=NULL
#confidence threshold for FAFB: ~0.62
clk_nt[clk_nt$nt_type_score<=0.62,]$nt_type = NA
clk_nt = as.data.frame(clk_nt)
#correct for nt_content based on T2A Gal4 expression of AChT, VGlut (Fig. 05F)
clk_nt[clk_nt$root_id %in% clk[clk$clk_name =="DN1pA",]$clk_id,]$nt_type = "GLUT"
clk_nt[clk_nt$root_id %in% clk[clk$clk_name =="LN_ITP",]$clk_id,]$nt_type = "ACH"

for (i in unique(clk_nt$root_id)) {
  no_nt = is.na(clk_nt[clk_nt$root_id ==i,]$nt_type)
  if (no_nt == T) {
  tmp_hemilineage = ifelse(clk_nt[clk_nt$root_id==i,]$hemilineage[1]%in% c("putative_primary","primary"),
                           yes=NA,no=clk_nt[clk_nt$root_id==i,]$hemilineage[1])
  tmp_cell_type = clk_nt[clk_nt$root_id==i,]$cell_type[1] 
  tmp_root_id = i
  
  if (!is.na(tmp_hemilineage)) {
    all_same_nt = length(na.omit(unique(clk_nt[clk_nt$hemilineage==tmp_hemilineage &
                                         clk_nt$root_id != tmp_root_id,]$nt_type)))==1
    if (all_same_nt == T) {
      tmp_nt = na.omit(clk_nt[clk_nt$hemilineage==tmp_hemilineage &
                        clk_nt$root_id != tmp_root_id ,]$nt_type)[1]
    }else{
      tmp_nt = NA
    }
  }
  if (is.na(tmp_hemilineage)|is.na(tmp_nt) & !is.na(tmp_cell_type)) {
    all_same_nt = length(na.omit(unique(clk_nt[clk_nt$cell_type==tmp_cell_type &
                                         clk_nt$root_id != tmp_root_id,]$nt_type)))==1
    if (all_same_nt == T) {
      tmp_nt = na.omit(clk_nt[clk_nt$cell_type==tmp_cell_type &
                        clk_nt$root_id != tmp_root_id ,]$nt_type)[1]
    }else{
      tmp_nt = NA
    }
  }
  clk_nt[clk_nt$root_id==tmp_root_id,]$nt_type = tmp_nt
  }
}
clk_nt$nt_type = ifelse(clk_nt$nt_type%in%c("SER","OCT","DA"),
                        yes=NA,no = clk_nt$nt_type)
write.csv(clk_nt,"./output/Figure_05/Figure_05_clk_singel_source-Data.csv")
clk_number = clk%>%
  group_by(clk_name)%>%
  summarize(n_total = length(clk_id))
clk_nt_grouped = clk_nt %>%
  group_by(clk_name,nt_type)%>%
  summarize(n= length(root_id))
clk_nt_grouped = left_join(clk_nt_grouped,clk_number,by= "clk_name")
clk_nt_grouped = as.data.frame(clk_nt_grouped)
clk_nt_grouped$perc = clk_nt_grouped$n/clk_nt_grouped$n_total

clk_nt_grouped$clk_name = factor(clk_nt_grouped$clk_name,
                                        levels = c("s-LNv","l-LNv","LN_ITP",
                                                   "LNd_CRYp","LNd_CRYn","LPN",
                                                   "l-CPDN3","APDN3","s-CPDN3A",
                                                   "s-CPDN3B","s-CPDN3C","s-CPDN3D",
                                                   "s-CPDN3E","DN2","DN1pE",
                                                   "DN1pD","DN1pC","DN1pB",
                                                   "DN1pA","DN1a"))

plot = ggplot(clk_nt_grouped,aes(x=perc,y=clk_name, fill=nt_type))+
  geom_col()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_manual(values  =  c("ACH" = "#eb3036ff",
                                 "GLUT" = "#53ab89ff",
                                 "NA" = "lightgrey"))+
  theme(legend.position = "",
        axis.text.x = element_blank())
ggsave(paste0(PATH_output,"Figure_05/Figure_05_D_clk.pdf"),width = 3,height=5,units = "cm")
source_data =plot$data
write.csv(source_data,"./output/Figure_05/Figure_05_clk_source-Data.csv")


if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
  source("./scripts/01_setup.R")
}
synapses_curated  =   read_csv(paste0(PATH_input,"/tmp/clk_neuron_input_filtered_v",v,".csv"),
                               col_types  =   cols(post_pt_supervoxel_id  =  col_character(),
                                                   post_pt_root_id  =  col_character(),
                                                   pre_pt_supervoxel_id  =  col_character(),
                                                   pre_pt_root_id  =  col_character()))

clk_input_sum = synapses_curated %>%
  group_by(post_pt_root_id, pre_pt_root_id) %>%
  summarise(n_synapses = length(pre_pt_root_id))
classification_join = classification
colnames(classification_join) = c("pre_pt_root_id",colnames(classification)[-1])
clk_input = left_join(clk_input_sum[clk_input_sum$n_synapses >=  5,],
                      classification_join[,c("pre_pt_root_id","hemilineage",
                                             "cell_type","hemibrain_type")],
                      by  = "pre_pt_root_id")
clk_input$cell_type = ifelse(is.na(clk_input$cell_type),yes = clk_input$hemibrain_type,no =clk_input$cell_type )
clk_input$hemibrain_type=NULL

colnames(nt) = c("pre_pt_root_id",colnames(nt)[-1])
clk_input_nt = left_join(clk_input,nt, by = "pre_pt_root_id")
clk_input_nt$group=NULL

clk_join = clk
colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
clk_input_nt = left_join(clk_input_nt,clk_join,by = "post_pt_root_id")

#confidence threshold for FAFB: ~0.62
clk_input_nt$nt_type = ifelse(clk_input_nt$nt_type_score<=0.62,yes = NA,
                              no = clk_input_nt$nt_type)
#correct for nt_content based on T2A Gal4 expression of AChT, VGlut
clk_input_nt[clk_input_nt$pre_pt_root_id %in% clk[clk$clk_name =="DN1pA",]$clk_id,]$nt_type = "GLUT"
clk_input_nt[clk_input_nt$pre_pt_root_id %in% clk[clk$clk_name =="LN_ITP",]$clk_id,]$nt_type = "ACH"

for (i in unique(clk_input_nt$pre_pt_root_id)) {
  no_nt = is.na(unique(clk_input_nt[clk_input_nt$pre_pt_root_id==i,]$nt_type))
  if (no_nt == T) {
    tmp_hemilineage = ifelse(clk_input_nt[clk_input_nt$pre_pt_root_id==i,]$hemilineage[1]%in% c("putative_primary","primary"),
                             yes=NA,no=clk_input_nt[clk_input_nt$pre_pt_root_id==i,]$hemilineage[1])
    tmp_cell_type = clk_input_nt[clk_input_nt$pre_pt_root_id==i,]$cell_type[1] 
    tmp_pre_pt_root_id =i
    
    if (!is.na(tmp_hemilineage)) {
      all_same_nt = length(na.omit(unique(clk_input_nt[clk_input_nt$hemilineage==tmp_hemilineage &
                                                   clk_input_nt$pre_pt_root_id != tmp_pre_pt_root_id,]$nt_type)))==1
      if (all_same_nt == T) {
        tmp_nt = na.omit(clk_input_nt[clk_input_nt$hemilineage==tmp_hemilineage &
                                  clk_input_nt$pre_pt_root_id != tmp_pre_pt_root_id ,]$nt_type)[1]
      }else{
        tmp_nt = NA
      }
    }
    if ((is.na(tmp_hemilineage)|is.na(tmp_nt))==T & !is.na(tmp_cell_type)) {
      all_same_nt = length(na.omit(unique(clk_input_nt[clk_input_nt$cell_type==tmp_cell_type &
                                                   clk_input_nt$pre_pt_root_id != tmp_pre_pt_root_id,]$nt_type)))==1
      if (all_same_nt == T) {
        tmp_nt = na.omit(clk_input_nt[clk_input_nt$cell_type==tmp_cell_type &
                                  clk_input_nt$pre_pt_root_id != tmp_pre_pt_root_id ,]$nt_type)[1]
      }else{
        tmp_nt = NA
      }
    }
    clk_input_nt[clk_input_nt$pre_pt_root_id==tmp_pre_pt_root_id,]$nt_type = tmp_nt
  }
}
clk_input_nt$nt_type = ifelse(clk_input_nt$nt_type%in%c("SER","OCT","DA"),yes=NA,no = clk_input_nt$nt_type)


clk_input_nt$clk = ifelse(clk_input_nt$pre_pt_root_id %in%clk$clk_id,yes="clk",no="no_clk")
clk_input_nt_grouped = clk_input_nt%>%
  group_by(name_post,nt_type,clk)%>%
  summarize(n_syn = sum(n_synapses),n_neurons = length(pre_pt_root_id))
clK_input_total = clk_input_nt%>%
  group_by(name_post)%>%
  summarize(n_syn_total = sum(n_synapses),n_neurons_total = length(pre_pt_root_id))
clk_input_nt_grouped = left_join(clk_input_nt_grouped,clK_input_total,by="name_post")
clk_input_nt_grouped$perc_syn = clk_input_nt_grouped$n_syn/clk_input_nt_grouped$n_syn_total
clk_input_nt_grouped$perc_neuron = clk_input_nt_grouped$n_neurons/clk_input_nt_grouped$n_neurons_total

clk_input_nt_grouped$name_post = factor(clk_input_nt_grouped$name_post,
                                        levels = c("s-LNv","l-LNv","LN_ITP",
                                                   "LNd_CRYp","LNd_CRYn","LPN",
                                                   "l-CPDN3","APDN3","s-CPDN3A",
                                                   "s-CPDN3B","s-CPDN3C","s-CPDN3D",
                                                   "s-CPDN3E","DN2","DN1pE",
                                                   "DN1pD","DN1pC","DN1pB",
                                                   "DN1pA","DN1a"))
clk_input_nt_grouped$nt_type_clk_id = paste(clk_input_nt_grouped$nt_type,clk_input_nt_grouped$clk,sep="_")
clk_input_nt_grouped$nt_type_clk_id = factor(clk_input_nt_grouped$nt_type_clk_id,
                                             levels= c("NA_no_clk","NA_clk","GABA_no_clk",
                                                       "GLUT_no_clk","GLUT_clk",
                                                       "ACH_clk","ACH_no_clk"))
plot = ggplot(clk_input_nt_grouped,aes(x=perc_syn,y=name_post, fill=nt_type_clk_id))+
  geom_bar_pattern(stat = "identity",col = "black",linewidth = 0.1,pattern = "stripe",
                   aes(pattern_density = nt_type_clk_id),pattern_fill = "white",
                   pattern_colour = NA,pattern_angle = 80,pattern_size = 0.01,
                   pattern_alpha = 1,pattern_key_scale_factor  =  1, pattern_spacing=0.02)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_pattern_density_manual(values  =
                                 c("GLUT_clk" = 0.1,
                                   "GLUT_no_clk" = 0,
                                   "ACH_clk" = 0.1,
                                   "ACH_no_clk" = 0,
                                   "GABA_no_clk" = 0,
                                   "NA_no_clk" = 0,
                                   "NA_clk" = 0.1
                                 ))+
  scale_fill_manual(values  =  c("ACH_clk" = "#eb3036ff",
                                 "ACH_no_clk" = "#eb3036ff",
                                 "GLUT_clk" = "#53ab89ff",
                                 "GLUT_no_clk" = "#53ab89ff",
                                 "GABA_no_clk" = "#4275b6ff",
                                 "NA_no_clk" = "black",
                                 "NA_clk" = "black"))+
  theme(legend.position = "",
        axis.text.x = element_blank())
ggsave(paste0(PATH_output,"Figure_05/Figure_05_D_pre_partner.pdf"),width = 7,height=5,units = "cm")
source_data =plot$data
write.csv(source_data,"./output/Figure_05/Figure_05_input_source-Data.csv")


if (!paste0("clk_neuron_output_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_preparation.py")
  source("./scripts/01_setup.R")
}
synapses_curated  =   read_csv(paste0(PATH_input,"/tmp/clk_neuron_output_filtered_v",v,".csv"),
                               col_types  =   cols(post_pt_supervoxel_id  =  col_character(),
                                                   post_pt_root_id  =  col_character(),
                                                   pre_pt_supervoxel_id  =  col_character(),
                                                   pre_pt_root_id  =  col_character()))

clk_output_sum = synapses_curated %>%
  group_by(post_pt_root_id, pre_pt_root_id) %>%
  summarise(n_synapses = length(pre_pt_root_id))
classification_join = classification
colnames(classification_join) = c("post_pt_root_id",colnames(classification)[-1])
clk_output = left_join(clk_output_sum[clk_output_sum$n_synapses >=  5,],
                      classification_join[,c("post_pt_root_id","hemilineage",
                                             "cell_type","hemibrain_type")],
                      by  = "post_pt_root_id")
clk_output$cell_type = ifelse(is.na(clk_output$cell_type),yes = clk_output$hemibrain_type,no =clk_output$cell_type )
clk_output$hemibrain_type=NULL

colnames(nt) = c("post_pt_root_id",colnames(nt)[-1])
clk_output_nt = left_join(clk_output,nt, by = "post_pt_root_id")
clk_output_nt$group=NULL

clk_join = clk
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_post")
clk_output_nt = left_join(clk_output_nt,clk_join,by = "pre_pt_root_id")

#confidence threshold for FAFB: ~0.62
clk_output_nt$nt_type = ifelse(clk_output_nt$nt_type_score<=0.62,yes = NA,no = clk_output_nt$nt_type)
#correct for nt_content based on T2A Gal4 expression of AChT, VGlut
clk_output_nt[clk_output_nt$post_pt_root_id %in% clk[clk$clk_name =="DN1pA",]$clk_id,]$nt_type = "GLUT"
clk_output_nt[clk_output_nt$post_pt_root_id %in% clk[clk$clk_name =="LN_ITP",]$clk_id,]$nt_type = "ACH"


for (i in unique(clk_output_nt$post_pt_root_id)) {
  no_nt = is.na(clk_output_nt[clk_output_nt$post_pt_root_id==i,]$nt_type[1])
  if (no_nt == T) {
    tmp_hemilineage = ifelse(clk_output_nt[clk_output_nt$post_pt_root_id==i,]$hemilineage[1]%in% c("putative_primary","primary"),yes=NA,no=clk_output_nt[clk_output_nt$post_pt_root_id==i,]$hemilineage[1])
    tmp_cell_type = clk_output_nt[clk_output_nt$post_pt_root_id==i,]$cell_type[1] 
    tmp_post_pt_root_id = i
    
    if (!is.na(tmp_hemilineage)) {
      all_same_nt = length(na.omit(unique(clk_output_nt[clk_output_nt$hemilineage==tmp_hemilineage &
                                                         clk_output_nt$post_pt_root_id != tmp_post_pt_root_id,]$nt_type)))==1
      if (all_same_nt == T) {
        tmp_nt = na.omit(clk_output_nt[clk_output_nt$hemilineage==tmp_hemilineage &
                                        clk_output_nt$post_pt_root_id != tmp_post_pt_root_id ,]$nt_type)[1]
      }else{
        tmp_nt = NA
      }
    }
    if (is.na(tmp_hemilineage)|is.na(tmp_nt) & !is.na(tmp_cell_type)) {
      all_same_nt = length(na.omit(unique(clk_output_nt[clk_output_nt$cell_type==tmp_cell_type &
                                                         clk_output_nt$post_pt_root_id != tmp_post_pt_root_id,]$nt_type)))==1
      if (all_same_nt == T) {
        tmp_nt = na.omit(clk_output_nt[clk_output_nt$cell_type==tmp_cell_type &
                                        clk_output_nt$post_pt_root_id != tmp_post_pt_root_id ,]$nt_type)[1]
      }else{
        tmp_nt = NA
      }
    }
    clk_output_nt[clk_output_nt$post_pt_root_id==tmp_post_pt_root_id,]$nt_type = tmp_nt
  }
}
clk_output_nt$nt_type = ifelse(clk_output_nt$nt_type%in%c("SER","OCT","DA"),
                               yes=NA,no = clk_output_nt$nt_type)


clk_output_nt$clk = ifelse(clk_output_nt$post_pt_root_id %in%clk$clk_id,yes="clk",no="no_clk")
clk_output_nt_grouped = clk_output_nt%>%
  group_by(name_pre,nt_type,clk)%>%
  summarize(n_syn = sum(n_synapses),n_neurons = length(post_pt_root_id))
clK_output_total = clk_output_nt%>%
  group_by(name_pre)%>%
  summarize(n_syn_total = sum(n_synapses),n_neurons_total = length(post_pt_root_id))
clk_output_nt_grouped = left_join(clk_output_nt_grouped,clK_output_total,by="name_pre")
clk_output_nt_grouped$perc_syn = clk_output_nt_grouped$n_syn/clk_output_nt_grouped$n_syn_total
clk_output_nt_grouped$perc_neuron = clk_output_nt_grouped$n_neurons/clk_output_nt_grouped$n_neurons_total

clk_output_nt_grouped$name_pre = factor(clk_output_nt_grouped$name_pre,
       levels = c("s-LNv","l-LNv","LN_ITP",
                      "LNd_CRYp","LNd_CRYn","LPN",
                      "l-CPDN3","APDN3","s-CPDN3A",
                      "s-CPDN3B","s-CPDN3C","s-CPDN3D",
                      "s-CPDN3E","DN2","DN1pE",
                      "DN1pD","DN1pC","DN1pB",
                      "DN1pA","DN1a"))
clk_output_nt_grouped$nt_type_clk_id = paste(clk_output_nt_grouped$nt_type,clk_output_nt_grouped$clk,sep="_")
clk_output_nt_grouped$nt_type_clk_id = paste(clk_output_nt_grouped$nt_type,clk_output_nt_grouped$clk,sep="_")
clk_output_nt_grouped$nt_type_clk_id = factor(clk_output_nt_grouped$nt_type_clk_id,
                                             levels= c("NA_no_clk","NA_clk","GABA_no_clk",
                                                       "GLUT_no_clk","GLUT_clk",
                                                       "ACH_clk","ACH_no_clk"))

plot = ggplot(clk_output_nt_grouped,aes(x=perc_syn,y=name_pre, fill=nt_type_clk_id))+
  geom_bar_pattern(stat = "identity",col = "black",linewidth = 0.1,pattern = "stripe",
                   aes(pattern_density = nt_type_clk_id),pattern_fill = "white",
                   pattern_colour = NA,pattern_angle = 80,pattern_size = 0.01,
                   pattern_alpha = 1,pattern_key_scale_factor  =  1, pattern_spacing=0.02)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_pattern_density_manual(values  =
                                 c("GLUT_clk" = 0.1,
                                   "GLUT_no_clk" = 0,
                                   "ACH_clk" = 0.1,
                                   "ACH_no_clk" = 0,
                                   "GABA_no_clk" = 0,
                                   "NA_no_clk" = 0,
                                   "NA_clk" = 0.1
                                 ))+
  scale_fill_manual(values  =  c("ACH_clk" = "#eb3036ff",
                                 "ACH_no_clk" = "#eb3036ff",
                                 "GLUT_clk" = "#53ab89ff",
                                 "GLUT_no_clk" = "#53ab89ff",
                                 "GABA_no_clk" = "#4275b6ff",
                                 "NA_no_clk" = "black",
                                 "NA_clk" = "black"))+
  theme(legend.position = "",
        axis.text.x = element_blank())
ggsave(paste0(PATH_output,"Figure_05/Figure_05_D_post_partner.pdf"),width = 7,height=5,units = "cm")
source_data =plot$data
write.csv(source_data,"./output/Figure_05/Figure_05_output_source-Data.csv")

