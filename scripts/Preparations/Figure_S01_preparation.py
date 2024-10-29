#----Figure_S01_preparation-----------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure S01 
# for the connectivity analysis for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

nuclei_df = client.materialize.query_table("nuclei_v1", materialization_version=v)
nuclei_df.to_csv("./input/tmp/nuclei_coords_v"+v+".csv")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del datastack_name
del client
del nuclei_df



