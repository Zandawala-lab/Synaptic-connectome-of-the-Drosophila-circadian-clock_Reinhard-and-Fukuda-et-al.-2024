#----Figure_01_F----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 01 Panel F for the connectivity analysis for 
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
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
grey_scale_tmp = grey_scale
grey_scale = colorRampPalette(c("lightgrey","black"))
grid.col_pre = NULL
DN1p_col = c("#648fff","#185ef0","#1B458F","#159fff")
#-------------------------------------------------------------------------------
synapses_sum = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
synapses_sum_clk = synapses_sum[synapses_sum$pre_pt_root_id%in%clk$clk_id,]

clk_join = clk
clk_join$clk_id_old = NULL
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")

synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "pre_pt_root_id")

colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")

synapses_sum_clk = left_join(synapses_sum_clk,clk_join,by = "post_pt_root_id")

clk_order = c("s-LNv","l-LNv","LN_ITP","LNd_CRYp",
             "LNd_CRYn","LPN","APDN3","s-CPDN3A",
             "s-CPDN3B","s-CPDN3C","s-CPDN3D",
             "s-CPDN3E","DN2","DN1pB","DN1pA",
             "DN1pC","DN1pD","DN1pE","DN1a")
clk_level = c(1:length(clk_order)) 

synapses_sum_clk$name_pre_order = synapses_sum_clk$name_pre
synapses_sum_clk$name_post_order = synapses_sum_clk$name_post

for (i in clk_level) {
  synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i]&
                     synapses_sum_clk$hemisphere_post == "left",]$name_post_order = 
    paste(letters[(length(clk_level)-(i-1))],
          synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i]&
                             synapses_sum_clk$hemisphere_post == "left",]$name_post,
          sep = "_")
  synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i]&
                     synapses_sum_clk$hemisphere_pre == "left",]$name_pre_order = 
    paste(letters[(length(clk_level)-(i-1))],
          synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i]&
                             synapses_sum_clk$hemisphere_pre == "left",]$name_pre,
          sep = "_")

  synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i]&
                     synapses_sum_clk$hemisphere_post == "right",]$name_post_order = 
    paste(letters[i],
          synapses_sum_clk[synapses_sum_clk$name_post == clk_order[i]&
                             synapses_sum_clk$hemisphere_post == "right",]$name_post,
          sep = "_")
  synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i]&
                     synapses_sum_clk$hemisphere_pre == "right",]$name_pre_order = 
    paste(letters[i],
          synapses_sum_clk[synapses_sum_clk$name_pre == clk_order[i]&
                             synapses_sum_clk$hemisphere_pre == "right",]$name_pre,
          sep = "_")
  }

synapses_sum_clk$pre_name_hemisphere = paste(synapses_sum_clk$hemisphere_pre,
                                            synapses_sum_clk$name_pre_order,
                                            synapses_sum_clk$pre_pt_root_id,sep = "_")
synapses_sum_clk$post_name_hemisphere = paste(synapses_sum_clk$hemisphere_post,
                                             synapses_sum_clk$name_post_order,
                                             synapses_sum_clk$post_pt_root_id,sep = "_")

synapses_sum_clk_grouped = synapses_sum_clk[synapses_sum_clk$name_pre%in%c("DN1pA")&
                                             synapses_sum_clk$hemisphere_pre == "right",] %>%
  group_by(post_name_hemisphere,pre_name_hemisphere,hemisphere_post,name_post)%>%
  summarise(n_synapses_sum = sum(n_synapses,na.rm = T),
            avrg_synapses = mean(n_synapses,na.rm = T),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)))

synapses_sum_clk_grouped = synapses_sum_clk_grouped[synapses_sum_clk_grouped$n_synapses_sum >= 5,]

grid.col = structure(grey_scale(length(unique(synapses_sum_clk_grouped$name_post))),
                    names = unique(synapses_sum_clk_grouped$name_post))
col = as.data.frame(grid.col)
col$names = row.names(col)
colnames(col) = c("grid.col","name_post")
synapses_sum_clk_grouped = left_join(synapses_sum_clk_grouped,col,by = "name_post")

grid.col_post = structure(synapses_sum_clk_grouped$grid.col,
                         names = synapses_sum_clk_grouped$post_name_hemisphere)

DN1p_sc = unique(synapses_sum_clk_grouped$pre_name_hemisphere)
for (i in c(1:length(DN1p_sc))) {
  grid.col_pre_tmp = structure(rep(DN1p_col[i],
                      times = length(
                        synapses_sum_clk_grouped[synapses_sum_clk_grouped$
                          pre_name_hemisphere == DN1p_sc[i],]$pre_name_hemisphere)),
                            names = rep(DN1p_sc[i],times = length(
                              synapses_sum_clk_grouped[synapses_sum_clk_grouped$
                                pre_name_hemisphere == DN1p_sc[i],]$pre_name_hemisphere)))
  grid.col_pre = c(grid.col_pre,grid.col_pre_tmp)
}

grid.col = c(grid.col_post,grid.col_pre)
grid.col_pre_names = unique(names(grid.col_pre))
grid.col_pre = as.data.frame(unique(grid.col_pre))
grid.col_pre$pre_name_hemisphere = grid.col_pre_names
colnames(grid.col_pre) = c("link_col","pre_name_hemisphere")
synapses_sum_clk_grouped = left_join(synapses_sum_clk_grouped,grid.col_pre,
                                    by = "pre_name_hemisphere")
sector_names_small_gap = c(unique(synapses_sum_clk_grouped$post_name_hemisphere),
                          unique(synapses_sum_clk_grouped$pre_name_hemisphere))
sector_names_big_gap = c("right_k_s-CPDN3D_720575940633954808")
sector_names_medium_gap = c("right_o_DN1pA_720575940633976403")
sector_names_medium_gap_II = c("left_q_LN_ITP_720575940634984800")
sector_names_small_gap = sector_names_small_gap[!sector_names_small_gap%in%
                                                 c(sector_names_big_gap,
                                                   sector_names_medium_gap,
                                                   sector_names_medium_gap_II)]
gap.after = c(structure(rep(2,times = length(sector_names_small_gap)),
                       names = sector_names_small_gap),
             structure(rep(80,times = length(sector_names_big_gap)),
                       names = sector_names_big_gap), #115
             structure(rep(40,times = length(sector_names_medium_gap)),
                       names = sector_names_medium_gap), #65
             structure(rep(40,times = length(sector_names_medium_gap_II)),
                       names = sector_names_medium_gap_II)
             )

circos.par(start.degree = 90,gap.after = gap.after)
pdf(paste0(PATH_output,"/Figure_01/Figure_01_F_",v,".pdf"),
    width = 5, height = 5)
chordDiagram(synapses_sum_clk_grouped[,c(2,1,4)],
             scale = T,
             col = synapses_sum_clk_grouped$link_col,
             grid.col = grid.col,
             directional = 1,
             direction.type = c( "arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(7),
             abline(v = 0, lty = 2, col = "#000000",lwd = 3),
             preAllocateTracks = list(),transparency = 0
)
circos.clear()
dev.off()
grey_scale = grey_scale_tmp