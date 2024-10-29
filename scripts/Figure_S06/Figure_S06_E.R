#----Figure_S06_E---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S06 Panel E for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
ipsi_gap = 2
contra_gap = 180
grid.col = NULL
grid.col = c("s-LNv" = s_LNv,"l-LNv" = l_LNv,"LN_ITP" = LN_ITP,
             "LNd_CRYp" = LNd_CRYp,"LNd_CRYn" = LNd_CRYn,"LPN" = LPN,
             "APDN3" = APDN3,"s-CPDN3B" = s_CPDN3B,"s-CPDN3A" = s_CPDN3A,
             "s-CPDN3D" = s_CPDN3A,"l-CPDN3" = l_CPDN3, "DN2" = DN2,
             "DN1pA" = DN1p_A,"DN1pB" = DN1p_B,"DN1pC" = DN1p_A,"DN1pD" = DN1p_A,
             "DN1pE" = DN1p_A, "DN1a" = DN1a)
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
                                       connectivity_merged$clk_id,sep = "_")

#correction of the names
for (i in unique(clk$clk_name)) {
  connectivity_merged[connectivity_merged$input %in%
                        clk[clk$clk_name == i,]$clk_id,]$name = i
}

connectivity_merged = connectivity_merged[!is.na(
  connectivity_merged$nr_synapses),]

connectivity_summarized = connectivity_merged %>%
  group_by(name,clk_name) %>%
  summarise(avrg_synapses = mean(nr_synapses,na.rm = T),
            nr_synapses = sum(nr_synapses,na.rm = T))

connectivity_summarized = connectivity_summarized[connectivity_summarized$name %in%
                     clk$clk_name & !is.na(connectivity_summarized$nr_synapses),]

col = as.data.frame(grid.col)
col$names = row.names(col)
colnames(col) = c("color","name")

connectivity_summarized = left_join(
  connectivity_summarized,col,by = "name")
connectivity_summarized$color = adjustcolor(
  connectivity_summarized$color,alpha.f = 0.5)
connectivity_summarized$avrg_synapses = NULL

circos.par(start.degree = 90,
           gap.after = c(
             "s-LNv"=contra_gap,"l-LNv"=ipsi_gap,"LN_ITP"=ipsi_gap,"LNd_CRYp"=ipsi_gap,
             "LNd_CRYn"=ipsi_gap,"LPN"=ipsi_gap,"APDN3"=ipsi_gap,"l-CPDN3" = ipsi_gap,
             "s-CPDN3D"=ipsi_gap,"s-CPDN3B"=ipsi_gap,"s-CPDN3A"=ipsi_gap,"DN2"=ipsi_gap,
             "DN1pE"=ipsi_gap,"DN1pD"=ipsi_gap,"DN1pC"=ipsi_gap,"DN1pB"=ipsi_gap,
             "DN1pA"=ipsi_gap,"DN1a"=ipsi_gap)
)
pdf(paste0(PATH_output,"Figure_S06/Figure_S06_E.pdf"),width = 10, height = 5)
chordDiagram(connectivity_summarized[c(1,2,3)],
             order = c("DN1a","DN1pE","DN1pD","DN1pC","DN1pB","DN1pA","DN2",
                       "s-CPDN3D","s-CPDN3B","s-CPDN3A","l-CPDN3","APDN3","LPN",
                       "LNd_CRYn","LNd_CRYp","LN_ITP","l-LNv","s-LNv"),
             col = connectivity_summarized$color,
             link.visible = connectivity_summarized$nr_synapses >= 10,
             grid.col = grid.col,
             directional = 1,
             direction.type = c("arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(7),
             abline(v = 0, lty = 2, col = "#000000",lwd = 3)
)
circos.clear()     
dev.off()  
