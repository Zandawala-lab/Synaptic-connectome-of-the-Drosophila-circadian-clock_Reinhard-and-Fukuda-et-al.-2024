#----Figure_03_E_preparation----------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure 03 Panel E
# for the connectivity analysis for Reinhard and Fukuda et al. 2024
# run first:
# 'preparation' section in "Figure_03_E.R"
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd
import navis
from fafbseg import flywire

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

pre_clk = pd.read_csv("./output/Figure_03/Figure_03_E_clk_input_80_v"+v+".csv")

OCG_id = pre_clk["pre_pt_root_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                filter_in_dict={"post_pt_root_id" : OCG_id},
                                materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                filter_in_dict={"target_id": syn_id},
                                merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_E_pre_clk_input_filtered_v"+v+".csv")




OCG_id = flywire.update_ids(id = [720575940609125641, 720575940622649633,
                              720575940632267112, 720575940633513645])["new_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                 filter_in_dict={"post_pt_root_id" : OCG_id},
                                 materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                  filter_in_dict={"target_id": syn_id},
                                  merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_E_pre_OCG02_input_filtered_v"+v+".csv")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del datastack_name
del client
del pre_clk
del OCG_id
del syn_df
del syn_id
del valid_syn_df
del valid_id
del result
