#----Figure_S01-----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S01 for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# 'morphology_cluster.csv' was last updated 05.2023 and is no longer available 
#                             via codex.org. It is based on nblast similarity 
#                             scores for 124,988/127,978 neurons
# nucleus segmentation from Mu et al. 2023 10.1101/2021.11.04.467197:
# https://neuromancer-seung-import.appspot.com/?local_id=2e3abb29a592113120b564a7f718273b#
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
source("./scripts/Preparations/Figure_S01_preparations_for_nuclei.R")
grid.col = c("s-LNv" = s_LNv,"l-LNv" = l_LNv,"LN_ITP" = LN_ITP,
             "LNd_CRYp" = LNd_CRYp,"LNd_CRYn" = LNd_CRYn,"LPN" = LPN,
             "l-CPDN3" = l_CPDN3,"APDN3" = APDN3,"s-CPDN3A" = s_CPDN3A,
             "s-CPDN3B" = s_CPDN3B,"s-CPDN3C" = s_CPDN3C,"s-CPDN3D" = s_CPDN3D,
             "s-CPDN3E" = s_CPDN3E,"DN2" = DN2,"DN1pB" = DN1p_B,
             "DN1pA" = DN1p_A,"DN1pC" = DN1p_A,"DN1pD" = DN1p_A,
             "DN1pE" = DN1p_A,"DN1a" = DN1a)
brain_surf = readOBJ(con = paste0(PATH_input,"FlyWire_brain_based_on_syn_Schlegel_et_al.obj"))
#-------------------------------------------------------------------------------
# plot clk nuclei locations:----------------------------------------------------
clk_with_nuclei$col= NA
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("DN1a")] = DN1a
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("DN1pA",
                                    "DN1pC", "DN1pD","DN1pE")] = DN1p_A
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("DN1pB")] = DN1p_B
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("DN2")] = DN2
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("APDN3")] = APDN3
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("l-CPDN3")] = l_CPDN3
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-CPDN3A")] = s_CPDN3A
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-CPDN3B")] = s_CPDN3B
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-CPDN3C")] = s_CPDN3C
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-CPDN3D")] = s_CPDN3D
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-CPDN3E")] = s_CPDN3E
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("LNd_CRYn")] = LNd_CRYn
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("LNd_CRYp")] = LNd_CRYp
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("LPN")] = LPN
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("LN_ITP")] = LN_ITP
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("s-LNv")] = s_LNv
clk_with_nuclei$col[clk_with_nuclei$clk_names %in% c("l-LNv")] = l_LNv

open3d()
plot3d(clk_with_nuclei$pt_position_x, clk_with_nuclei$pt_position_y,
       clk_with_nuclei$pt_position_z, col=clk_with_nuclei$col, 
       size = 0.5, add=T,type = "s")
plot3d(brain_surf, alpha = 0.025, add=T, col="black")
view3d(userMatrix = rotationMatrix(180*pi/180, 1, 0, 0),zoom=0.405)
s_CPDN3_nuclei = clk_with_nuclei[clk_with_nuclei$clk_names %in% 
                    c("s-CPDN3A","s-CPDN3B","s-CPDN3C","s-CPDN3D","s-CPDN3E"),]
clk_with_nuclei$col = NULL
#-------------------------------------------------------------------------------
warning(c("the following clk neurons were not included in morphology clustering:", 
      unique(clk[!clk$clk_id %in% morph_cluster_clk$pt_root_id,]$clk_name)))
morph_clk_all_nuclei$col= NA
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("DN1a"),]$clk_id] = DN1a
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("DN1pA","DN1pC", "DN1pD","DN1pE"),]$clk_id] = DN1p_A
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("DN1pB"),]$clk_id] = DN1p_B
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("DN2"),]$clk_id] = DN2
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("APDN3"),]$clk_id] = APDN3
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("l-CPDN3"),]$clk_id] = l_CPDN3
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-CPDN3A"),]$clk_id] = s_CPDN3A
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-CPDN3B"),]$clk_id] = s_CPDN3B
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-CPDN3C"),]$clk_id] = s_CPDN3C
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-CPDN3D"),]$clk_id] = s_CPDN3D
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-CPDN3E"),]$clk_id] = s_CPDN3E
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("LNd_CRYn"),]$clk_id] = LNd_CRYn
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("LNd_CRYp"),]$clk_id] = LNd_CRYp
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("LPN"),]$clk_id] = LPN
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("LN_ITP"),]$clk_id] = LN_ITP
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("s-LNv"),]$clk_id] = s_LNv
morph_clk_all_nuclei$col[morph_clk_all_nuclei$pt_root_id %in% clk[clk$clk_name %in% c("l-LNv"),]$clk_id] = l_LNv
morph_clk_all_nuclei$col[is.na(morph_clk_all_nuclei$col)] = "red"

clk_join = clk
colnames(clk_join) = c("clk_names","pt_root_id", "hemisphere")
morph_clk_all_nuclei = left_join(morph_clk_all_nuclei,clk_join)

morph_clust_plot = morph_clk_all_nuclei %>%
  group_by(morphology_cluster, clk_names, col) %>%
  summarize(n_neurons = length(pt_root_id))
morph_clk_all_nuclei[,c("clk_names","hemisphere")] =NULL

plot = ggplot(morph_clust_plot, aes(x=morphology_cluster,y = n_neurons, fill= clk_names))+
  geom_col()+
  scale_fill_manual(values = grid.col,na.value = "grey")+
  theme(
    panel.background = element_rect(fill = "white",colour = "black"),
    panel.border = element_rect(color = "black",fill = NA),
    legend.position = "right",
    axis.title.y = element_blank(),
    axis.text.y = element_text(vjust = 0.5,margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text.x = element_text(size = 8,angle = 45, hjust = 1)
  )+
  guides(fill=guide_legend(ncol=2))
ggsave(paste0(PATH_output,"Figure_S01/Figure_S01_morph_clust.pdf"),
       height = 10,width = 20, units = "cm")

morph_clust_s_CPDN3 = morph_clk_all_nuclei[morph_clk_all_nuclei$morphology_cluster %in% c("MC_37.161","MC_401.91")&
                                             !morph_clk_all_nuclei$pt_root_id %in% clk$clk_id,]
open3d()
plot3d(s_CPDN3_nuclei$pt_position_x, s_CPDN3_nuclei$pt_position_y,
       s_CPDN3_nuclei$pt_position_z, col=s_CPDN3_nuclei$col, 
       size = 0.5, add=T,type = "s")
plot3d(morph_clust_s_CPDN3$pt_position_x, morph_clust_s_CPDN3$pt_position_y,
       morph_clust_s_CPDN3$pt_position_z, col="red", 
       size = 0.5, add=T,type = "s")
plot3d(brain_surf, alpha = 0.025, add=T, col="black")
view3d(userMatrix = rotationMatrix(180*pi/180, 1, 0, 0),zoom=0.405)

clk_cluster = c("DN1a","DN1pB","DN1pC","DN1pD","DN1pE","s-CPDN3A","s-CPDN3B",
                "s-CPDN3C","s-CPDN3D","s-CPDN3E","APDN3","LPN")
clk_cand = matrix(nrow = 0, ncol = 5)
colnames(clk_cand) = c("pt_root_id_cand", "dist", "pt_root_id_clk", "clk_name", "hemisphere")
for (i in clk_cluster) {
  tmp_clk_name = i
  tmp_clk_root_id = clk[clk$clk_name %in% i,]$clk_id
  tmp_morph = unique(morph_cluster_clk[morph_cluster_clk$pt_root_id %in%
                                     tmp_clk_root_id,]$morphology_cluster)
  tmp_morph_cluster = as.data.frame(morph_clk_all_nuclei[
      morph_clk_all_nuclei$morphology_cluster %in% tmp_morph &
      !is.na(morph_clk_all_nuclei$pt_position_x),
      c("pt_root_id","pt_position_x","pt_position_y","pt_position_z")])
  row_names = tmp_morph_cluster$pt_root_id
  row.names(tmp_morph_cluster) = row_names
  tmp_morph_cluster$pt_root_id = NULL
  dist_morph_clust_w = dist(tmp_morph_cluster)
  dist_morph_clust = gather(as.data.frame(as.matrix(dist_morph_clust_w)))
  dist_morph_clust$pt_root_id = row_names
  colnames(dist_morph_clust) = c("pt_root_id","dist","pt_root_id.y")
  colnames(clk_join) = c("clk_name","pt_root_id","hemisphere")
  dist_morph_clust = left_join(dist_morph_clust,clk_join,by= "pt_root_id")
  colnames(clk_join) = c("clk_name.y","pt_root_id.y","hemisphere.y")
  dist_morph_clust = left_join(dist_morph_clust,clk_join, by= "pt_root_id.y")
  dist_morph_clust = dist_morph_clust[dist_morph_clust$dist!=0,]
  avrg_dist_clk_left = mean(na.rm =T,dist_morph_clust[dist_morph_clust$clk_name %in% tmp_clk_name & 
                                               dist_morph_clust$hemisphere == "left" & 
                                               dist_morph_clust$clk_name.y %in% tmp_clk_name & 
                                               dist_morph_clust$hemisphere.y == "left",]$dist)
  avrg_dist_clk_right = mean(na.rm =T, dist_morph_clust[dist_morph_clust$clk_name %in% tmp_clk_name & 
                                                dist_morph_clust$hemisphere == "right" & 
                                                dist_morph_clust$clk_name.y %in% tmp_clk_name & 
                                                dist_morph_clust$hemisphere.y == "right",]$dist)
  avrg_dist_clk = mean(x = c(avrg_dist_clk_right,avrg_dist_clk_left))
  
  clk_cand_tmp = dist_morph_clust[!dist_morph_clust$pt_root_id %in% clk$clk_id,]
  clk_cand_tmp = clk_cand_tmp[clk_cand_tmp$pt_root_id.y %in% clk[clk$clk_name %in% tmp_clk_name,]$clk_id &
                        clk_cand_tmp$dist <= avrg_dist_clk*2,]
  clk_cand_tmp$clk_name=NULL
  clk_cand_tmp$hemisphere=NULL
  colnames(clk_cand_tmp) = c("pt_root_id_cand", "dist", "pt_root_id_clk","clk_name","hemisphere")
  
    if (length(clk_cand_tmp[!is.na(clk_cand_tmp$pt_root_id_cand),]$pt_root_id_cand)>0) {
    clk_cand = rbind(clk_cand,clk_cand_tmp)
    root_id_for_plot = c(clk_cand_tmp$pt_root_id_cand, clk[clk$clk_name %in% unique(clk_cand_tmp$clk_name),]$clk_id)
    morph_clk_all_nuclei_plot = morph_clk_all_nuclei[morph_clk_all_nuclei$pt_root_id %in% root_id_for_plot,]
    open3d()
    plot3d(morph_clk_all_nuclei_plot$pt_position_x, morph_clk_all_nuclei_plot$pt_position_y,
               morph_clk_all_nuclei_plot$pt_position_z, col=morph_clk_all_nuclei_plot$col, 
               size = 0.5, add=T,type = "s")
    plot3d(brain_surf, alpha = 0.025, add=T, col="black")
    view3d(userMatrix = rotationMatrix(180*pi/180, 1, 0, 0),zoom=0.405)
    }
}
write.csv(clk_cand,paste0(PATH_output,"clk_candidates_based_on_morph_and_nucleus_location.csv"))

# # code for selecting single neurons manualy:----------------------------------
# root_id_for_plot = clk[clk$clk_name %in% c("DN1pC"),]$clk_id
# morph_clk_all_nuclei_plot = morph_clk_all_nuclei[
#   morph_clk_all_nuclei$morphology_cluster %in% unique(
#     morph_clk_all_nuclei[morph_clk_all_nuclei$pt_root_id %in%
#                            root_id_for_plot,]$morphology_cluster),]
# open3d()
# y = plot3d(morph_clk_all_nuclei_plot$pt_position_x, morph_clk_all_nuclei_plot$pt_position_y,
#            morph_clk_all_nuclei_plot$pt_position_z, col=morph_clk_all_nuclei_plot$col, 
#            size = 0.5, add=T,type = "s")
# plot3d(FAFB14.surf, alpha = 0.025, add=T, col="black")
# candidates = selectpoints3d(object = y)
# candidate_pt_root_id = matrix(nrow = 0, ncol = 1)
# colnames(candidate_pt_root_id) = "pt_root_id"
# for (i in c(1:length(candidates[,1]))) {
#   tmp_x_index = match(candidates[i,1],morph_clk_all_nuclei_plot$pt_position_x)
#   if (length(tmp_x_index)>1) {
#     tmp_y_index = match(candidates[i,2],morph_clk_all_nuclei_plot$pt_position_y)
#     if (length(tmp_y_index)>1) {
#       tmp_z_index = match(candidates[i,3],morph_clk_all_nuclei_plot$pt_position_z)
#       index = Reduce(intersect, list(tmp_x_index,tmp_y_index,tmp_z_index))
#     }else{
#       index = Reduce(intersect, list(tmp_x_index,tmp_y_index))
#     }
#   }else{
#     index = tmp_x_index
#   }
#   if (!is.na(morph_clk_all_nuclei_plot[index,]$pt_root_id)) {
#     candidate_pt_root_id = rbind(candidate_pt_root_id, 
#                                  morph_clk_all_nuclei_plot[index,]$pt_root_id)
#   }
# }
# #-------------------------------------------------------------------------------