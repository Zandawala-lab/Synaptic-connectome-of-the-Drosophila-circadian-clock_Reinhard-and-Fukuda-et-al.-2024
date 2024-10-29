#----Figure_S11_S12-------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure S11/ S12 for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("clk_neuron_input_filtered_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
}
synapses_input = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
                          col_types = cols(pre_pt_supervoxel_id = col_character(),
                                           pre_pt_root_id = col_character(),
                                           post_pt_supervoxel_id = col_character(),
                                           post_pt_root_id = col_character()))

if (!paste0("clk_neuron_output_filtered_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_04_preparations.py")
}
clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
synapses_output =  read_csv(paste0(PATH_input,"tmp/clk_neuron_output_filtered_v",v,".csv"),
                            col_types =  cols(pre_pt_supervoxel_id = col_character(),
                                              pre_pt_root_id = col_character(),
                                              post_pt_supervoxel_id = col_character(),
                                              post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
synapses_input$pre_x = xyzmatrix(synapses_input$pre_pt_position)[,1]
synapses_input$pre_y = xyzmatrix(synapses_input$pre_pt_position)[,2]
synapses_input$pre_z = xyzmatrix(synapses_input$pre_pt_position)[,3]
synapses_input$post_x = xyzmatrix(synapses_input$post_pt_position)[,1]
synapses_input$post_y = xyzmatrix(synapses_input$post_pt_position)[,2]
synapses_input$post_z = xyzmatrix(synapses_input$post_pt_position)[,3]

for (i in unique(clk$clk_name)) {
open3d()
plot3d(synapses_input[synapses_input$post_pt_root_id %in% clk[clk$clk_name %in% 
                              c(i) & clk$hemisphere=="right",]$clk_id,]$post_x,
       synapses_input[synapses_input$post_pt_root_id %in% clk[clk$clk_name %in% 
                              c(i) & clk$hemisphere=="right",]$clk_id,]$post_y,
       synapses_input[synapses_input$post_pt_root_id %in% clk[clk$clk_name %in% 
                              c(i) & clk$hemisphere=="right",]$clk_id,]$post_z,
       add=T,size = 0.5,type = "s")
view3d(zoom=1)
  writeOBJ(paste0(PATH_output,"Figure_S11_12/input_syn_",i,"_v",v,".obj"))
close3d()
}

synapses_output$pre_x = xyzmatrix(synapses_output$pre_pt_position)[,1]
synapses_output$pre_y = xyzmatrix(synapses_output$pre_pt_position)[,2]
synapses_output$pre_z = xyzmatrix(synapses_output$pre_pt_position)[,3]
synapses_output$post_x = xyzmatrix(synapses_output$post_pt_position)[,1]
synapses_output$post_y = xyzmatrix(synapses_output$post_pt_position)[,2]
synapses_output$post_z = xyzmatrix(synapses_output$post_pt_position)[,3]

for (i in unique(clk$clk_name)) {
  open3d()
  plot3d(synapses_output[synapses_output$pre_pt_root_id %in% clk[clk$clk_name %in% 
                                c(i) & clk$hemisphere=="right",]$clk_id,]$pre_x,
         synapses_output[synapses_output$pre_pt_root_id %in% clk[clk$clk_name %in% 
                                c(i) & clk$hemisphere=="right",]$clk_id,]$pre_y,
         synapses_output[synapses_output$pre_pt_root_id %in% clk[clk$clk_name %in% 
                                c(i) & clk$hemisphere=="right",]$clk_id,]$pre_z,
         add=T,size=0.5, type ="s")
  view3d(zoom=1)
  writeOBJ(paste0(PATH_output,"Figure_S11_12/output_syn_",i,"_v",v,".obj"))
  close3d()
}
