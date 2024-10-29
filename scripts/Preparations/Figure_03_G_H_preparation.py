#----Figure_03_G,H_preparation--------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure 03 Panel G and H
# for the connectivity analysis for Reinhard and Fukuda et al. 2024
# run first:
# 'preparation' section in "Figure_03_G_H.R"
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])
clk_input = pd.read_csv("./input/tmp/clk_neuron_input_filtered_v"+v+".csv")
bridge_cand_id = np.unique(clk_input["pre_pt_root_id"])

R1_6_id = pd.read_csv("./input/tmp/Figure_03_R1-6_ids_v"+v+".csv")

R1_6_id = R1_6_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1",
                                        filter_in_dict={"pre_pt_root_id": R1_6_id,
                                        "post_pt_root_id": bridge_cand_id },
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2",
                                        filter_in_dict={"target_id": syn_id},
                                        merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_R1_6_output_filtered_v"+v+".csv")


R7_id = pd.read_csv("./input/tmp/Figure_03_R7_ids_v"+v+".csv")

R7_id = R7_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                        filter_in_dict={"pre_pt_root_id": R7_id,
                                        "post_pt_root_id": bridge_cand_id},
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                        filter_in_dict={"target_id": syn_id},
                                        merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_R7_output_filtered_v"+v+".csv")

R8_id = pd.read_csv("./input/tmp/Figure_03_R8_ids_v"+v+".csv")

R8_id = R8_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                        filter_in_dict={"pre_pt_root_id": R8_id,
                                        "post_pt_root_id": bridge_cand_id},
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                        filter_in_dict={"target_id": syn_id},
                                        merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_R8_output_filtered_v"+v+".csv")

OC_id = pd.read_csv("./input/tmp/Figure_03_ocellar_retinula_cell_ids_v"+v+".csv")

OC_id = OC_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                       filter_in_dict={"pre_pt_root_id": OC_id},
                                       materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                       filter_in_dict={"target_id": syn_id},
                                       merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_OC_output_filtered_v"+v+".csv")



HB_id = pd.read_csv("./input/tmp/Figure_03_HB_ids_v"+v+".csv")

HB_id = HB_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                        filter_in_dict={"pre_pt_root_id": HB_id,
                                        "post_pt_root_id": bridge_cand_id},
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                        filter_in_dict={"target_id": syn_id},
                                        merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_HB_output_filtered_v"+v+".csv")


OCG_id = pd.read_csv("./input/tmp/Figure_03_OCG02c_ids_v"+v+".csv")

OCG_id = OCG_id["x"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                        filter_in_dict={"post_pt_root_id": OCG_id},
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                        filter_in_dict={"target_id": syn_id},
                                        merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_OCG_input_filtered_v"+v+".csv")

syn_df = client.materialize.query_table("synapses_nt_v1", 
                                        filter_in_dict={"pre_pt_root_id": OCG_id,
                                        "post_pt_root_id": bridge_cand_id},
                                        materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                                  filter_in_dict={"target_id": syn_id},
                                                  merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/Figure_03_GH_OCG_output_filtered_v"+v+".csv")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del datastack_name
del client
del clk_input
del bridge_cand_id
del R1_6_id
del OC_id
del R7_id
del R8_id
del HB_id
del OCG_id
del syn_df
del syn_id
del valid_syn_df
del valid_id
del result
