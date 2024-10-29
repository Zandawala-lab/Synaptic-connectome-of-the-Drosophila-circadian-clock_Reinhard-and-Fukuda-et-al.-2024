#----Figure_06_A----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 06 Panel A for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# load libraries needed for scRNAseq analysis ----------------------------------
library(Seurat)
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
set.seed(1313)
load(paste0(PATH_input,"scRNAseq_data_clk.rda"))
#-------------------------------------------------------------------------------
#clustering of clock neurons... cluster names were updated in final figure to 
# match the cluster names used in the rest of the manuscript
pdf(paste0(PATH_output,"Figure_06/Figure_06_A.pdf"),
    width = 10, height = 10)
DimPlot(sub.clusters, 
        reduction = 'tsne', 
        group.by = 'new_clusters', 
        pt.size = 0.5,
        label = T,
        label.box = TRUE,
        repel = TRUE,
        label.size = 4) +
  NoLegend() +
  labs(title = NULL)
dev.off()  