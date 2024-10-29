#----Figure_S01_preparations_for_nuclei------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure S01
# for the connectivity analysis for Reinhard and Fukuda et al. 
# nucleus segmentation from Mu et al. 2023 10.1101/2021.11.04.467197:
# https://neuromancer-seung-import.appspot.com/?local_id=2e3abb29a592113120b564a7f718273b#
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("nuclei_coords_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_S01_preparation.py")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
nuclei_coords = read_delim(paste0(PATH_input,"tmp/nuclei_coords_v",v,".csv"), 
                           delim = ",", escape_double = FALSE, 
                           col_types =cols(pt_root_id = col_character(), 
                                           id = col_character(),
                                           pt_supervoxel_id = col_character()),
                           trim_ws = TRUE)
# morphology_cluster are only available for v630!
if (!"morphology_clusters_v630.csv" %in% input_files) {
  stop("please make sure the file 'morphology_clusters_v630.csv' is located in './input'.")
}
morph_cluster<-read_delim(paste0(PATH_input,"morphology_clusters_v630.csv"), 
                          delim = ",", escape_double = FALSE,
                          col_types = cols(root_id = col_character(),
                                           morphology_cluster = col_character()), 
                          trim_ws = TRUE)
#-------------------------------------------------------------------------------
duplicated_root_id = nuclei_coords[duplicated(nuclei_coords$pt_root_id) &
                                     nuclei_coords$pt_root_id != 0,]$pt_root_id
nuclei_coords_duplicated = nuclei_coords[nuclei_coords$pt_root_id %in% 
                              duplicated_root_id,c("id","volume", "pt_root_id")]

valid_duplicates = nuclei_coords_duplicated %>%
  group_by(pt_root_id) %>%
  top_n(1,volume)
valid_id = c(valid_duplicates$id, nuclei_coords[!nuclei_coords$pt_root_id %in%
                                                  duplicated_root_id,]$id)
nuclei_coords = nuclei_coords[nuclei_coords$id %in% valid_id,]

nuclei_coords$pt_position_x = xyzmatrix(nuclei_coords$pt_position)[,1]
nuclei_coords$pt_position_y = xyzmatrix(nuclei_coords$pt_position)[,2]
nuclei_coords$pt_position_z = xyzmatrix(nuclei_coords$pt_position)[,3]

clk_join = clk
colnames(clk_join) = c("clk_names","pt_root_id", "hemisphere")
clk_with_nuclei = left_join(clk_join,nuclei_coords, by= "pt_root_id")

# determine which clk neurons did not get a nucleus assigned
# and assign it manually:-------------------------------------------------------
clk_without_nuclei = clk_with_nuclei[is.na(clk_with_nuclei$id),]$pt_root_id
clk_without_nuclei_v630 = flywire_latestid(clk_without_nuclei,version = 630)
# write.csv(clk_without_nuclei_v630, paste0(PATH_output,"clk_neurons_without_assigned_nucleus_v630.csv"))

if (!"clk_neurons_manual_assigned_nucleus.csv" %in% input_files) {
  stop("please make sure the file 'clk_neurons_manual_assigned_nucleus.csv' is located in './input'.")
}

clk_manual_nucleus = read_delim(paste0(PATH_input,"clk_neurons_manual_assigned_nucleus.csv"),
                                col_types = cols(pt_root_id = col_character(),
                                                 id = col_character()),delim = ";")
# clk_without_nuclei_v630 = clk_without_nuclei_v630[!clk_without_nuclei_v630 %in% clk_manual_nucleus$pt_root_id]
#write.csv(clk_without_nuclei_v630, paste0(PATH_output,"clk_neurons_without_assigned_nucleus_v630.csv"))
clk_manual_nucleus$pt_root_id = flywire_latestid(clk_manual_nucleus$pt_root_id,version = v)
clk_manual_nucleus = left_join(clk_manual_nucleus,nuclei_coords,by = "id")
clk_manual_nucleus$pt_root_id.y = NULL
colnames(clk_manual_nucleus) = c("pt_root_id",colnames(clk_manual_nucleus)[-1])
clk_manual_nucleus = left_join(clk_manual_nucleus,clk_join,by = "pt_root_id")

clk_with_nuclei = rbind(clk_with_nuclei[!is.na(clk_with_nuclei$id),],clk_manual_nucleus)
write.csv(paste0(PATH_output,"clk_with_nuclei_v",v,".csv"))
#-------------------------------------------------------------------------------

# determin which neurons of the needed morphology clusters did not get a nucleus 
# assigned and assign it manually:----------------------------------------------
colnames(morph_cluster) = c("pt_root_id", "morphology_cluster")
clk_id_v630 = flywire_latestid(clk$clk_id,version = 630)
tmp_morph_cluster_clk = unique(morph_cluster[morph_cluster$pt_root_id %in% 
                                               clk_id_v630,]$morphology_cluster)

morph_cluster_clk = morph_cluster[morph_cluster$morphology_cluster %in% tmp_morph_cluster_clk,]
morph_cluster_clk$pt_root_id = flywire_latestid(morph_cluster_clk$pt_root_id,version = v)

morph_clk_all_nuclei = left_join(morph_cluster_clk, nuclei_coords)

clk_nuclei_morph = left_join(clk_with_nuclei[,-c(1,3)], morph_cluster_clk)

morph_clk_all_nuclei = rbind(morph_clk_all_nuclei [!morph_clk_all_nuclei$pt_root_id %in%
                                                     clk$clk_id, ], clk_nuclei_morph)

morph_clk_all_without_nuclei = morph_clk_all_nuclei[is.na(morph_clk_all_nuclei$id),]$pt_root_id
# write.csv(morph_clk_all_without_nuclei, paste0(PATH_output,"morph_clk_all_without_assigned_nucleus.csv"))

if (!"morph_clk_all_manual_assigned_nucleus.csv" %in% input_files) {
  stop("please make sure the file 'morph_clk_all_manual_assigned_nucleus.csv' is located in './input'.")
}

morph_clk_all_manual_nuclei = read_delim(paste0(PATH_input,"morph_clk_all_manual_assigned_nucleus.csv"),
                                  col_types = cols(pt_root_id = col_character(),
                                  id = col_character()),delim = ";")

# morph_clk_all_without_nuclei_v630 = flywire_latestid(morph_clk_all_without_nuclei,version = 630)
# morph_clk_all_without_nuclei_v630 = morph_clk_all_without_nuclei[!morph_clk_all_without_nuclei %in% morph_clk_all_manual_nuclei$pt_root_id]
#problem in mapping these ("720575940635878042" "720575940615062447") to version. They are always twice included...
#write.csv(morph_clk_all_without_nuclei_v630, paste0(PATH_output,"morph_clk_all_without_assigned_nucleus_v630.csv"))

morph_clk_all_manual_nuclei$pt_root_id = flywire_latestid(morph_clk_all_manual_nuclei$pt_root_id,version = v)
#neuron fragments without cell body and accordingly no nucleus will be excluded for further analysis -> this removes single fragments from the morphology clusters... 
morph_clk_all_manual_nuclei = left_join(morph_clk_all_manual_nuclei[!is.na(morph_clk_all_manual_nuclei$id),],nuclei_coords,by = "id")
morph_clk_all_manual_nuclei$pt_root_id.y = NULL
colnames(morph_clk_all_manual_nuclei) = c("pt_root_id",colnames(morph_clk_all_manual_nuclei)[-1])
morph_clk_all_manual_nuclei = left_join(morph_clk_all_manual_nuclei,morph_cluster_clk,by = "pt_root_id")

morph_clk_all_nuclei = unique(rbind(morph_clk_all_nuclei[!is.na(morph_clk_all_nuclei$id),],
                                    morph_clk_all_manual_nuclei))

#-------------------------------------------------------------------------------

