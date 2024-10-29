#----Figure_06_C----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 06 Panel C for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# load libraries needed for scRNAseq analysis ----------------------------------
library(Seurat)
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
set.seed(1313)
load(paste0(PATH_input,"scRNAseq_data_clk.rda"))
#-------------------------------------------------------------------------------
pdf(paste0(PATH_output,"Figure_06/Figure_06_C.pdf"),
    width = 30, height = 10)
FeaturePlot(sub.clusters, features = c("AstA","AstC","CCHa1","CNMa","Dh31","Dh44",
                                       "ITP","NPF","Pdf","Proc","sNPF","Trissin"),
            slot = "scale.data", reduction = "tsne", min.cutoff = "q10",
            max.cutoff = "q90", ncol = 6, pt.size = 0.1, 
            cols = c("#CCCCCC","#FF9900","#FF0000")) &
  theme(plot.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.size = unit (3,"mm"), 
        panel.background = element_rect(colour = "black", linewidth = 1, fill = NA),
        legend.box.spacing = unit(1,"mm")) & NoAxes()  
dev.off()  

