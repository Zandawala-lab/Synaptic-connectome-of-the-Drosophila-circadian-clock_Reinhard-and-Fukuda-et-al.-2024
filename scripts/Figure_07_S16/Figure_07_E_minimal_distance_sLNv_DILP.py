#----Figure_07_S16_preparation--------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code which has to be run in preparation for Figure 07 Panel E and 
# S18 for the connectivity analysis for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])
### define function-------------------------------------------------------------
def distance_filter_neuron(coords, root_ids, r_um=14*1000): #v783 skeleton coordinates are in nm while v680 are in µm
    neurons = root_ids
    print(neurons.shape, neurons.dtype)
    neuron_pt_kdtree = spatial.cKDTree(coords)
    clustered_neuron_pt = neuron_pt_kdtree.query_ball_point(coords, r=r_um) 
    close_pt = set()
    far_pt = []
    for i_cl in range(len(clustered_neuron_pt)):
        if len(clustered_neuron_pt[i_cl]) > 1: 
            local_neuron_pt = np.array(clustered_neuron_pt[i_cl]) 
            neuron_status = neurons[local_neuron_pt] != neurons[i_cl] 
            for id_ in local_neuron_pt[neuron_status]: 
                close_pt.add(id_)
            
            close_pt.add(i_cl) if sum(neuron_status) >0 else far_pt.append(i_cl)
        else: 
            far_pt.append(i_cl)
    far_pt = np.array(far_pt)
    close_pt = np.array(list(close_pt))
    assert len(far_pt) + len(close_pt) == len(clustered_neuron_pt)
    return close_pt, far_pt
#-------------------------------------------------------------------------------
df = pd.read_csv("./input/tmp/tmp_neuronlist_v"+v+".csv")
coords = np.array(df[["X", "Y", "Z"]], dtype="int64")
neurons = np.array(df["group"])
pt_ids = np.array(df["PointNo"], dtype=np.uint64)
unit = "nm"
r_um = 1000
close_pt = []
while len(close_pt)<1:
  close_pt, far_pt = distance_filter_neuron(coords, neurons, r_um=r_um)
  if len(close_pt)<1:
    r_um = r_um +1
print("minimal distance between s-LNv and m-NSC_DILP:", r_um/1000, "µm")
result = pd.DataFrame(data = [[r_um,unit]], index = [1], columns = ["distance","unit"])
if len(close_pt)>0:
  close_pt_ids = pt_ids[close_pt]
else:
  close_pt_ids = []
  
if len(far_pt)>0:
  far_pt_ids = pt_ids[far_pt]
else:
  far_pt_ids = []

result.to_csv("./input/distance_for_paracrine_signaling_um_v"+v+".csv")

df_close = df.query("PointNo in @close_pt_ids")

#df_close.to_csv("./output/close_pt_ids.csv")

df_far = df.query("PointNo in @far_pt_ids")

#df_far.to_csv("./output/far_pt_ids.csv")
