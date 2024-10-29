#----Figure_01_E_F_preparation--------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure 01 Panel E and F
# for the connectivity analysis for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

clk_py = pd.read_csv(("./input/clk_neurons_v"+v+".csv"), sep=",")

clk_id = clk_py["clk_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                  filter_in_dict={"post_pt_root_id" : clk_id},
                                  materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                  filter_in_dict={"target_id": syn_id},
                                  merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/clk_neuron_input_filtered_v"+v+".csv")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del client
del clk_id
del clk_py
del datastack_name
del result
del syn_df
del valid_id
del valid_syn_df
del syn_id
