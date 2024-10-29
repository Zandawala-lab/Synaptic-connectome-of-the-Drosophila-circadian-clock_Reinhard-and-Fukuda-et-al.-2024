#----Figure_04_K_preparation-----------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure 04 K
# for the connectivity analysis for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

clk_output = pd.read_csv("./input/tmp/clk_output_sum_filtered_v"+v+".csv")
bridge_cand_id = np.unique(clk_output["post_pt_root_id"])

NSC_py = pd.read_csv("./input/NSC_v"+v+".csv", sep=",")

NSC_id = NSC_py["NSC_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                              filter_in_dict={"pre_pt_root_id" : bridge_cand_id, 
                              "post_pt_root_id" : NSC_id},
                              materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                              filter_in_dict={"target_id": syn_id},
                              merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_04_K_NSC_input_from_clk_bridge_filtered_v"+v+".csv")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del datastack_name
del client
del clk_output
del bridge_cand_id
del NSC_py
del NSC_id
del syn_df
del syn_id
del valid_syn_df
del valid_id
del result
