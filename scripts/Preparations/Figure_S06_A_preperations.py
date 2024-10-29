#----Figure_S06A_preparation----------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure S06 Panel A to 
# obtain the clock neuron meshes used in Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import navis
import cloudvolume
import pandas as pd

clk_py = pd.read_csv(("./input/clk_neurons_hemibrain_v1-2.csv"), sep=",")
clk_id = clk_py["clk_id"]

navis.patch_cloudvolume()

vol = cloudvolume.CloudVolume('precomputed://gs://neuroglancer-janelia-flyem-hemibrain/v1.2/segmentation', use_https=True, progress=False)
 
m = vol.mesh.get(clk_id, as_navis=True, lod =2)

navis.write_mesh(m, "./output/hemibrain_mesh_raw/",filetype="obj")

#clean up the Python environment, as otherwise variables may interfere with the R environment

del clk_py
del clk_id
del vol
del m
