#----Figure_S06_D---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S06 Panel D for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables-----------------------------------------------------------------
clk = read.csv(paste0(PATH_input,"clk_neurons_hemibrain_v1-2.csv"))
#-------------------------------------------------------------------------------
clk_count = clk%>%
  group_by(clk_name)%>%
  summarise(count = length(unique(clk_id)))
clk_count$clk_name = factor(clk_count$clk_name,
                            levels = c("DN1a","DN1pA","DN1pB","DN1pC","DN1pD",
                                     "DN1pE","DN2","s-CPDN3A","s-CPDN3B",
                                     "s-CPDN3D","l-CPDN3","APDN3","LPN",
                                     "LNd_CRYn","LNd_CRYp","LN_ITP","l-LNv",
                                     "s-LNv"))
ggplot(clk_count,aes(y = clk_name,x = count))+
  geom_col(aes(fill = clk_name))+
  scale_fill_manual(values = grid.col,limits = rev)+
  scale_y_discrete(position = "left",limits = rev)+
  scale_x_continuous(expand = c(0,0))+
  theme(
    panel.background = element_rect(fill = "white",colour = "black"),
    panel.border = element_rect(color = "black",fill = NA),
    legend.position = "",
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0.5,margin = margin(t = 0,
                                                    r = 15, b = 0, l = 0)),
    axis.text = element_text(size = 8)
  )
ggsave(paste0(PATH_output,"/Figure_S06/Figure_S06_D.pdf"),width = 8,
       height = 5, unit="cm")
