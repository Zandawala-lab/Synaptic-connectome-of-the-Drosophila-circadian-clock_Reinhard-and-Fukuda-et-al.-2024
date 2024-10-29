#----Figure_S03-----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S03 for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
}

clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")

synapses_curated = read_csv(paste0(PATH_input,"/tmp/clk_neuron_input_filtered_v",v,".csv"),
                            col_types = cols(pre_pt_supervoxel_id = col_character(),
                                             pre_pt_root_id = col_character(),
                                             post_pt_supervoxel_id = col_character(),
                                             post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
synapses_sum<-synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses=length(post_pt_root_id))
synapses_sum_clk<-synapses_sum[synapses_sum$pre_pt_root_id%in%
                                 clk$clk_id &synapses_sum$n_synapses>=1,]

clk_join<-clk
colnames(clk_join)<-c("name_pre","pre_pt_root_id","hemisphere_pre")
synapses_sum_clk<-left_join(synapses_sum_clk,clk_join,by = "pre_pt_root_id")

colnames(clk_join)<-c("name_post","post_pt_root_id","hemisphere_post")
synapses_sum_clk<-left_join(synapses_sum_clk,clk_join,by = "post_pt_root_id")

synapses_sum_clk$pre_name_hemisphere<-paste(synapses_sum_clk$hemisphere_pre,
                                            synapses_sum_clk$name_pre,sep="_")
synapses_sum_clk$post_name_hemisphere<-paste(synapses_sum_clk$hemisphere_post,
                                             synapses_sum_clk$name_post,sep="_")

synapses_sum_clk_grouped<-synapses_sum_clk %>%
  group_by(post_name_hemisphere,pre_name_hemisphere,hemisphere_post)%>%
  summarise(n_synapses_sum=sum(n_synapses,na.rm = T),
            avrg_synapses=mean(n_synapses,na.rm = T),
            n_pre_partners=length(unique(pre_pt_root_id)),
            n_post_partners=length(unique(post_pt_root_id)))

synapses_sum_clk_grouped$col<-ifelse(synapses_sum_clk_grouped$avrg_synapses>5,
                                     "white","black")

synapses_sum_clk_grouped$post_name_hemisphere<-factor(
  synapses_sum_clk_grouped$post_name_hemisphere,
  levels = c("left_s-LNv","left_l-LNv","left_LN_ITP","left_LNd_CRYp",
             "left_LNd_CRYn","left_LPN","left_l-CPDN3","left_APDN3",
             "left_s-CPDN3A","left_s-CPDN3B","left_s-CPDN3C","left_s-CPDN3D",
             "left_s-CPDN3E","left_DN2","left_DN1pA","left_DN1pB","left_DN1pC",
             "left_DN1pD","left_DN1pE","left_DN1a",
             "right_s-LNv","right_l-LNv","right_LN_ITP","right_LNd_CRYp",
             "right_LNd_CRYn","right_LPN","right_l-CPDN3","right_APDN3",
             "right_s-CPDN3A","right_s-CPDN3B","right_s-CPDN3C",
             "right_s-CPDN3D","right_s-CPDN3E","right_DN2","right_DN1pA",
             "right_DN1pB","right_DN1pC","right_DN1pD","right_DN1pE",
             "right_DN1a") )

synapses_sum_clk_grouped$pre_name_hemisphere<-factor(
  synapses_sum_clk_grouped$pre_name_hemisphere,
  levels = c("left_s-LNv","left_l-LNv","left_LN_ITP","left_LNd_CRYp",
             "left_LNd_CRYn","left_LPN","left_l-CPDN3","left_APDN3",
             "left_s-CPDN3A","left_s-CPDN3B","left_s-CPDN3C","left_s-CPDN3D",
             "left_s-CPDN3E","left_DN2","left_DN1pA","left_DN1pB","left_DN1pC",
             "left_DN1pD","left_DN1pE","left_DN1a",
             "right_s-LNv","right_l-LNv","right_LN_ITP","right_LNd_CRYp",
             "right_LNd_CRYn","right_LPN","right_l-CPDN3","right_APDN3",
             "right_s-CPDN3A","right_s-CPDN3B","right_s-CPDN3C","right_s-CPDN3D",
             "right_s-CPDN3E","right_DN2","right_DN1pA","right_DN1pB",
             "right_DN1pC","right_DN1pD","right_DN1pE","right_DN1a") )

ggplot(synapses_sum_clk_grouped,aes(x=pre_name_hemisphere,y=post_name_hemisphere,
                                    fill=avrg_synapses),color="black")+
  geom_raster()+
  geom_text(aes(label =n_synapses_sum,colour=col),show.legend=FALSE,size=(8/(1/0.35)))+ # 1/0.35 is the ratio between points used by theme() and mm used by geom_text()
  scale_colour_manual(values=c("black","white"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  labs(y="postsynaptic",x="presynaptic")+
  scale_fill_gradientn(colors = grey_scale(length(synapses_sum_clk_grouped$n_synapses_sum)),trans="log",na.value = "white")+
  geom_vline(xintercept = c(20.5),linewidth=1)+
  geom_vline(xintercept = c(0,seq(1.5,(length(unique(synapses_sum_clk_grouped$pre_name_hemisphere))+0.5),by=1)),linewidth=0.01)+
  geom_hline(yintercept = c(20.5),linewidth=1)+
  geom_hline(yintercept = c(0,seq(1.5,(length(unique(synapses_sum_clk_grouped$post_name_hemisphere))+0.5),by=1)),linewidth=0.01)+
  theme(axis.text.x = element_text( angle=90,hjust = 0,vjust = 0.5, colour = "black"),
        axis.text.y = element_text( angle=0,hjust = 0,colour = "black"),
        text = element_text(size=8),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S03/Figure_S03_v",v,".pdf"), width=8.5,height=6.5,units=c("in"))
