#----02_update_root_ids---------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the file for updating root ids of clk neurons and NSC used in the 
# connectivity analysis for Reinhard and Fukuda et al. 2024 
# This script has to be run if the analysis should be repeated with a different 
# version.
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
clk  =  read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                   col_types  =  cols(clk_id  =  col_character()),delim  =  ",")
NSC = read_delim(paste0(PATH_input,"NSC_v",v,".csv"),
                 col_types = cols(NSC_id = col_character()),delim = ",")
#-------------------------------------------------------------------------------
clk$root_id = flywire_latestid(rootid = clk$clk_id, version = v)
clk$clk_id = clk$root_id
clk$root_id = NULL
write.csv(clk,paste0(PATH_input,"clk_neurons_v",v,".csv"),row.names = F)

NSC$root_id = flywire_latestid(rootid = NSC$NSC_id, version = v)
NSC$NSC_id = NSC$root_id
NSC$root_id = NULL
write.csv(NSC,paste0(PATH_input,"NSC_v",v,".csv"),row.names = F)

# check single root_ids:
#flywire_latestid(rootid = "720575940634984800", version = v)
