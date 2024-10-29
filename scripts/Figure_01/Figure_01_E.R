#----Figure_01_E----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 01 Panel E for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
ipsi_gap = 2
contra_gap = 10
grid.col = NULL

if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
  source("./scripts/01_setup.R")
}
synapses_curated = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",
                                   v,".csv"),
                            col_types = cols(pre_pt_supervoxel_id = col_character(),
                                            pre_pt_root_id = col_character(),
                                            post_pt_supervoxel_id = col_character(),
                                            post_pt_root_id = col_character()))
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
#-------------------------------------------------------------------------------
#setting the color for the following plot---------------------------------------
grid.col = c("left_s-LNv" = s_LNv,"left_l-LNv" = l_LNv,"left_LN_ITP" = LN_ITP,
           "left_LNd_CRYp" = LNd_CRYp,"left_LNd_CRYn" = LNd_CRYn,"left_LPN" = LPN,
           "left_l-CPDN3" = l_CPDN3,"left_APDN3" = APDN3,"left_s-CPDN3A" = s_CPDN3A,
           "left_s-CPDN3B" = s_CPDN3B,"left_s-CPDN3C" = s_CPDN3C,
           "left_s-CPDN3D" = s_CPDN3D, "left_s-CPDN3E" = s_CPDN3E,
           "left_DN2" = DN2,"left_DN1pB" = DN1p_B, "left_DN1pA" = DN1p_A,
           "left_DN1pC" = DN1p_A,"left_DN1pD" = DN1p_A, "left_DN1pE" = DN1p_A,
           "left_DN1a" = DN1a,
           "right_DN1a" = DN1a,"right_DN1pE" = DN1p_A, "right_DN1pD" = DN1p_A,
           "right_DN1pC" = DN1p_A,"right_DN1pA" = DN1p_A,"right_DN1pB" = DN1p_B,
           "right_DN2"= DN2,"right_s-CPDN3E" = s_CPDN3E, 
           "right_s-CPDN3D" = s_CPDN3D,"right_s-CPDN3C" = s_CPDN3C,
           "right_s-CPDN3B" = s_CPDN3B, "right_s-CPDN3A" = s_CPDN3A,
           "right_APDN3" = APDN3,"right_l-CPDN3" = l_CPDN3, "right_LPN" = LPN,
           "right_LNd_CRYn" = LNd_CRYn,"right_LNd_CRYp" = LNd_CRYp,
           "right_LN_ITP" = LN_ITP,"right_l-LNv" = l_LNv, "right_s-LNv" = s_LNv)
#-------------------------------------------------------------------------------
synapses_sum = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
synapses_sum = synapses_sum[synapses_sum$n_synapses >= 5,]
synapses_sum_clk = synapses_sum[synapses_sum$pre_pt_root_id %in% clk$clk_id,]

clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")

synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "pre_pt_root_id")

colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")

synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "post_pt_root_id")

synapses_sum_clk$pre_name_hemisphere = paste(synapses_sum_clk$hemisphere_pre,
                                            synapses_sum_clk$name_pre,sep = "_")
synapses_sum_clk$post_name_hemisphere = paste(synapses_sum_clk$hemisphere_post,
                                             synapses_sum_clk$name_post,sep = "_")


synapses_sum_clk_grouped = synapses_sum_clk %>%
  group_by(post_name_hemisphere,pre_name_hemisphere,hemisphere_post) %>%
  summarise(n_synapses_sum = sum(n_synapses,na.rm = T),
            avrg_synapses = mean(n_synapses,na.rm = T),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)))


synapses_sum_clk_grouped$avrg_synapses = NULL
synapses_sum_clk_grouped$hemisphere_post = NULL
synapses_sum_clk_grouped$col = NULL

col = as.data.frame(grid.col)
col$names=row.names(col)
colnames(col) = c("color","pre_name_hemisphere")
synapses_sum_clk_grouped = left_join(synapses_sum_clk_grouped,col,
                                    by = "pre_name_hemisphere")

synapses_sum_clk_grouped$color = adjustcolor(synapses_sum_clk_grouped$color,
                                            alpha.f = 0.5)

synapses_sum_clk_grouped$color[synapses_sum_clk_grouped$pre_name_hemisphere %in%
                                 c("left_DN1pA","right_DN1pA")] = DN1p_A

circos.par(start.degree = -102.5,
           gap.after = c(
             "left_s-LNv" = ipsi_gap,"left_l-LNv" = ipsi_gap,
             "left_LN_ITP" = ipsi_gap,"left_LNd_CRYp" = ipsi_gap,
             "left_LNd_CRYn" = ipsi_gap,"left_LPN" = ipsi_gap,
             "left_l-CPDN3" = ipsi_gap, "left_APDN3" = ipsi_gap,
             "left_s-CPDN3A" = ipsi_gap, "left_s-CPDN3B" = ipsi_gap,
             "left_s-CPDN3C" = ipsi_gap, "left_s-CPDN3D" = ipsi_gap,
             "left_s-CPDN3E" = ipsi_gap, "left_DN2" = ipsi_gap, 
             "left_DN1pB" = ipsi_gap, "left_DN1pA" = ipsi_gap, 
             "left_DN1pC" = ipsi_gap, "left_DN1pD" = ipsi_gap,
             "left_DN1pE" = ipsi_gap, "left_DN1a" = contra_gap,
             
             "right_DN1a" = ipsi_gap, "right_DN1pE" = ipsi_gap,
             "right_DN1pD" = ipsi_gap,"right_DN1pC" = ipsi_gap,
             "right_DN1pA" = ipsi_gap, "right_DN1pB" = ipsi_gap,
             "right_DN2" = ipsi_gap, "right_s-CPDN3E" = ipsi_gap, 
             "right_s-CPDN3D" = ipsi_gap, "right_s-CPDN3C" = ipsi_gap,
             "right_s-CPDN3B" = ipsi_gap, "right_s-CPDN3A" = ipsi_gap, 
             "right_APDN3" = ipsi_gap, "right_l-CPDN3" = ipsi_gap,
             "right_LPN" = ipsi_gap, "right_LNd_CRYn" = ipsi_gap,
             "right_LNd_CRYp" = ipsi_gap, "right_LN_ITP" = ipsi_gap,
             "right_l-LNv" = ipsi_gap, "right_s-LNv" = contra_gap+5
           )
)
pdf(paste0(PATH_output,"Figure_01/Figure_01_E_",v,".pdf"),
    width = 5, height = 5)
chordDiagram(synapses_sum_clk_grouped[,c(2,1,3)],
             order = c("left_s-LNv","left_l-LNv","left_LN_ITP","left_LNd_CRYp",
                     "left_LNd_CRYn","left_LPN","left_l-CPDN3","left_APDN3",
                     "left_s-CPDN3A","left_s-CPDN3B","left_s-CPDN3C",
                     "left_s-CPDN3D", "left_s-CPDN3E","left_DN2", "left_DN1pA",
                     "left_DN1pB","left_DN1pC", "left_DN1pD","left_DN1pE",
                     "left_DN1a",
                     "right_DN1a","right_DN1pE","right_DN1pD","right_DN1pC",
                     "right_DN1pB","right_DN1pA","right_DN2","right_s-CPDN3E",
                     "right_s-CPDN3D","right_s-CPDN3C","right_s-CPDN3B",
                     "right_s-CPDN3A","right_APDN3","right_l-CPDN3","right_LPN",
                     "right_LNd_CRYn","right_LNd_CRYp","right_LN_ITP",
                     "right_l-LNv","right_s-LNv"),
             col = synapses_sum_clk_grouped$color,
             link.visible = synapses_sum_clk_grouped$n_synapses_sum >= 5,
             grid.col = grid.col,
             directional = 1,
             direction.type = c( "arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(7),
             abline(v = -0.07, lty = 2, col = "#000000",lwd = 3),
             preAllocateTracks = list()
)
circos.clear()
dev.off()  