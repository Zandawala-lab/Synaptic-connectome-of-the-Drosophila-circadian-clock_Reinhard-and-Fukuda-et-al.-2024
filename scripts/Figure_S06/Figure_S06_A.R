#----Figure_S06_A---------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S06 Panel A for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
neurons = list.files(path = paste0(PATH_output,"hemibrain_mesh_raw/"),
                     full.names = FALSE, recursive = FALSE)
for (i in neurons) {
  neuron_tmp = readOBJ(paste0(PATH_output,"hemibrain_mesh_raw/",i))
  neuron_corrected = neuron_tmp/8 #scales coordinates right
  neuron_form = xform_brain(neuron_corrected*8/1000, reference= "FAFB14",
                      sample="JRCFIB2018F", via = JRC2018F)
  #neuron_form = mirror_fafb(neuron_form) #use if you want to mirror neurons
  open3d()
  plot3d(neuron_form,add = T)
  writeOBJ(paste0("./output/hemibrain_meshes/",i),
           withNormals = TRUE, withTextures = TRUE,
           separateObjects = TRUE,ids = NULL)
  close3d()
}


