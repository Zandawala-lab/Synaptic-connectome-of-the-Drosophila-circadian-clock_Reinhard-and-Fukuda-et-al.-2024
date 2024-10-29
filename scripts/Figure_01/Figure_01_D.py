#----Figure_01_D ---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 01 Panel D for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# This is the code to obtain the clock neuron meshes used
# in Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import navis
import cloudvolume
import pandas as pd
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

clk_py = pd.read_csv(("./input/clk_neurons_v"+v+".csv"), sep=",")
clk_id = clk_py["clk_id"]

navis.patch_cloudvolume()

vol = cloudvolume.CloudVolume('graphene://https://prodv1.flywire-daf.com/segmentation/1.0/fly_v31',
use_https=True, progress=False)
 
m = vol.mesh.get(clk_id, as_navis=True)

navis.write_mesh(m, "./output/neurons_v"+v+"/clk",filetype="obj")

#clean up the Python environment, as otherwise variables may interfere with the R environment
del clk_py
del clk_id
del vol
del m

