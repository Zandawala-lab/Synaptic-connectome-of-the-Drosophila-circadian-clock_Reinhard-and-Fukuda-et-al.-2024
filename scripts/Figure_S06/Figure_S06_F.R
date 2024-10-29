#----Figure_S06_F---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S06 Panel F for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
clk = read.csv(paste0(PATH_input,"clk_neurons_hemibrain_v1-2.csv"))
#-------------------------------------------------------------------------------
connectivity = as.data.frame(neuprint_simple_connectivity(bodyids = clk$clk_id,
                                                         prepost = "PRE"))

colnames(connectivity) = c("input","name","type",clk$clk_id)
connectivity_long = melt(as.data.table(connectivity),
                        id.vars = c("input","name","type"))

colnames(connectivity_long) = c("input","name","type","clk_id","nr_synapses")
connectivity_long$input = as.character(connectivity_long$input)

connectivity_merged = merge.data.frame(connectivity_long,clk,by = "clk_id")
connectivity_merged$name_clk_id = paste(connectivity_merged$clk_name,
                                       connectivity_merged$clk_id,sep="_")

#correction of the names
for (i in unique(clk$clk_name)) {
  connectivity_merged[connectivity_merged$input %in%
                        clk[clk$clk_name == i,]$clk_id,]$name = i
}

connectivity_merged = connectivity_merged[!is.na(
  connectivity_merged$nr_synapses),]

connectivity_summarized = connectivity_merged %>%
  group_by(name,clk_name) %>%
  summarise(avrg_synapses = mean(nr_synapses),
            nr_synapses = sum(nr_synapses,na.rm = T))

colnames(connectivity_summarized) = c("name","clk_name",
                                     "avrg_synapses","nr_synapses_name")
connectivity_merged_sum = merge.data.frame(connectivity_merged,
                                          connectivity_summarized,
                                          by =c("clk_name","name"))
connectivity_merged_sum$nr_synapses_name[
  connectivity_merged_sum$nr_synapses_name == 0] = NA
connectivity_merged_sum$col = ifelse(connectivity_merged_sum$avrg_synapses < 10,
                                    "black","white")

connectivity_merged_sum$clk_name = factor(connectivity_merged_sum$clk_name,
                                         levels = rev(c("DN1a","DN1pE","DN1pD",
                                                    "DN1pC","DN1pB","DN1pA",
                                                    "DN2","s-CPDN3D",
                                                    "s-CPDN3B","s-CPDN3A","l-CPDN3",
                                                    "APDN3","LPN","LNd_CRYn",
                                                    "LNd_CRYp","LN_ITP","l-LNv",
                                                    "s-LNv")))
connectivity_merged_sum$name = factor(connectivity_merged_sum$name,
                                     levels = rev(c("DN1a","DN1pE","DN1pD",
                                                    "DN1pC","DN1pB","DN1pA",
                                                    "DN2","s-CPDN3D",
                                                    "s-CPDN3B","s-CPDN3A","l-CPDN3",
                                                    "APDN3","LPN","LNd_CRYn",
                                                    "LNd_CRYp","LN_ITP","l-LNv",
                                                    "s-LNv")))
connectivity_merged_sum$avrg_synapses = round(
  connectivity_merged_sum$avrg_synapses,digits = 1)

ggplot(connectivity_merged_sum[connectivity_merged_sum$input
                               %in% clk$clk_id,],
       aes(x = name,y = clk_name),color="black")+
  geom_raster(aes(fill = avrg_synapses))+
  geom_text(aes(label =nr_synapses_name,colour = col),
            show.legend = FALSE,size = (8/(1/0.35)))+
  scale_colour_manual(values = c("black","white"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  labs(y = "postsynaptic",x = "presynaptic")+
  scale_fill_gradientn(colors = grey_scale(length(
    connectivity_merged_sum$avrg_synapses)),trans = "log",na.value = "white")+
  geom_vline(xintercept = c(0,seq(1.5,(length(unique(
    connectivity_merged_sum$name))+0.5),by = 1)),linewidth = 0.01)+
  geom_hline(yintercept = c(0,seq(1.5,(length(unique(
    connectivity_merged_sum$clk_name))+0.5),by = 1)),linewidth = 0.01)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0,
                                   vjust = 0.5, colour = "black"),
        axis.text.y = element_text(angle = 0,hjust = 0,colour = "black"),
        panel.background = element_rect(fill="white",colour = "black"),
        axis.text = element_text(size = 8,color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = ""
        )
ggsave(paste0(PATH_output,"Figure_S06/Figure_S06F.pdf"),
       width = 12,height = 12,units = c("cm"))
