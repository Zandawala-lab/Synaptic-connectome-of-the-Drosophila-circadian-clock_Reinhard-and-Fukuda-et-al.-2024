#----Figure_06/07_preparations--------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code is needed as preparation for Figure 06/7 as preparation to run 
# the analysis for the connectivity analysis for Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# load libraries needed for scRNAseq analysis ----------------------------------
library(Seurat)
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
load(paste0(PATH_input,"scRNAseq_data_clk.rda"))
load(paste0(PATH_input,"scRNAseq_data_NSC.rda"))
#-------------------------------------------------------------------------------
#expression of neuropeptides and their receptors in clock neurons for dotplot---
genes_to_plot = c("tim","per","Clk", "cry", "Rh7","gl","Pdf","AstC","CNMa","ITP",
                   "Dh31","sNPF","AstA","Dh44","NPF","Trissin",
                   "CCHa1","Proc","Gpb5","Mip","Ms","Hug","Tk","AstA-R1","AstC-R1",
                   "AstC-R2","CCHa1-R","CNMaR","Dh31-R","Dh44-R1","Dh44-R2",
                   "NPFR","Pdfr","sNPF-R","TrissinR","Proc-R")
genes_fig_06 = c("tim","per","Clk", "cry", "Rh7","gl","Pdf","AstC","CNMa","ITP",
                 "Dh31","sNPF","AstA","Dh44","NPF","Trissin",
                 "CCHa1","Proc","Gpb5","Mip","Ms","Hug","Tk")
genes_fig_07 = c("AstA-R1","AstC-R1","AstC-R2","CCHa1-R","CNMaR","Dh31-R","Dh44-R1","Dh44-R2",
                 "NPFR","Pdfr","sNPF-R","TrissinR","Proc-R")

exp_mat = as.matrix(sub.clusters[["integrated"]]@scale.data[genes_to_plot,])
meta = sub.clusters@meta.data %>% 
  select(reorder_clusters)
meta = bind_cols(meta, as.data.frame(t(exp_mat)))
meta = pivot_longer(meta, -reorder_clusters, names_to="Gene",
                     values_to="Expression")
meta_summary_avrg_exp = meta %>%
  group_by(reorder_clusters, Gene) %>%
  summarise(avrg_expression = mean(Expression))

exp_mat = as.matrix(sub.clusters[["integrated"]]@data[genes_to_plot,])
meta = sub.clusters@meta.data %>% 
  select(reorder_clusters)
meta = bind_cols(meta, as.data.frame(t(exp_mat)))
meta = pivot_longer(meta, -reorder_clusters, names_to="Gene",
                    values_to="Expression")
meta_summary_perc_exp = meta %>%
  group_by(reorder_clusters, Gene) %>%
  summarise(pct_expression = sum(Expression > 0) / length(Expression) * 100)

meta_summary = left_join(meta_summary_avrg_exp, meta_summary_perc_exp, by=c("reorder_clusters","Gene"))

# updating cluster names to match the rest of the analysis:---------------------
meta_summary$reorder_clusters =  as.character(meta_summary$reorder_clusters)
meta_summary[meta_summary$reorder_clusters %in% c("1llnv"),]$reorder_clusters = "l-LNv"
meta_summary[meta_summary$reorder_clusters %in% c("2slnv"),]$reorder_clusters = "s-LNv"
meta_summary[meta_summary$reorder_clusters %in% c("3LN_ITP"),]$reorder_clusters = "LN_ITP"
meta_summary[meta_summary$reorder_clusters %in% c("4LNd_NPF"),]$reorder_clusters = "LNd_NPF"
meta_summary[meta_summary$reorder_clusters %in% c("5LNd_Trissin"),]$reorder_clusters = "LNd_Trissin"
meta_summary[meta_summary$reorder_clusters %in% c("6LPN"),]$reorder_clusters = "LPN"
meta_summary[meta_summary$reorder_clusters %in% c("7DN3"),]$reorder_clusters = "DN3"
meta_summary[meta_summary$reorder_clusters %in% c("8DN3_VGlut"),]$reorder_clusters = "DN3_VGlut"
meta_summary[meta_summary$reorder_clusters %in% c("9DN2"),]$reorder_clusters = "DN2"
meta_summary[meta_summary$reorder_clusters %in% c("10DN1p_a"),]$reorder_clusters = "DN1p_sNPF"
meta_summary[meta_summary$reorder_clusters %in% c("11DN1p_b"),]$reorder_clusters = "DN1p"
meta_summary[meta_summary$reorder_clusters %in% c("12DN1p_AstA"),]$reorder_clusters = "DN1p_AstA"
meta_summary[meta_summary$reorder_clusters %in% c("13DN1p_CNMa"),]$reorder_clusters = "DN1p_CNMa"
meta_summary[meta_summary$reorder_clusters %in% c("14DN1p_CNMaAstC"),]$reorder_clusters = "DN1p_CNMa_AstC"
meta_summary[meta_summary$reorder_clusters %in% c("15DN1p_Rh7"),]$reorder_clusters = "DN1p_Rh7"
meta_summary[meta_summary$reorder_clusters %in% c("16DN1a"),]$reorder_clusters = "DN1a"
colnames(meta_summary) = c("cluster","gene",colnames(meta_summary)[-c(1:2)])
#-------------------------------------------------------------------------------
write.csv(meta_summary[meta_summary$gene %in% genes_fig_06,], file = paste0(PATH_input,"tmp/",
            "neuropeptide_scale_data_summary_single_cell.csv"), row.names = F)
write.csv(meta_summary[meta_summary$gene %in% genes_fig_07,], file = paste0(PATH_input,"tmp/",
           "receptors_scaled_data_summary_single_cell.csv"),row.names = F)
#-------------------------------------------------------------------------------
#receptors for clock neuropeptides in NSC---------------------------------------
genes_to_plot = c("Crz","sNPF","Dh44","CG13248","CG13743","Ilp2","Ilp3","Ilp5",
                    "Tk","ITP","Dh31","ImpL2","AstA-R1","AstA-R2","AstC-R1",
                    "AstC-R2","CCHa1-R","CNMaR","Dh31-R","Dh44-R1","Dh44-R2",
                    "NPFR","Pdfr","Proc-R","sNPF-R","TrissinR")

exp_mat = as.matrix(NSCs[["RNA"]]@scale.data[genes_to_plot,])
meta = NSCs@meta.data %>% 
  select(new_cluster)
meta = bind_cols(meta, as.data.frame(t(exp_mat)))
meta = pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary_avrg_exp = meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(avrg_expression = mean(Expression))

exp_mat = as.matrix(NSCs[["RNA"]]@scale.data[genes_to_plot,])
meta = NSCs@meta.data %>% 
  select(new_cluster)
meta = bind_cols(meta, as.data.frame(t(exp_mat)))
meta = pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary_perc_exp = meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(pct_expression = sum(Expression > 0) / length(Expression) * 100)

meta_summary = left_join(meta_summary_avrg_exp, meta_summary_perc_exp, by=c("new_cluster","Gene"))
# updating cluster names to match the rest of the analysis:---------------------
meta_summary[meta_summary$new_cluster %in% c("CN"),]$new_cluster = "l-NSC_CRZ"
meta_summary[meta_summary$new_cluster %in% c("DH44"),]$new_cluster = "m-NSC_DH44"
meta_summary[meta_summary$new_cluster %in% c("IPCs"),]$new_cluster = "m-NSC_DILP"
meta_summary[meta_summary$new_cluster %in% c("l-NSC-DH31"),]$new_cluster = "l-NSC_DH31"
meta_summary[meta_summary$new_cluster %in% c("l-NSC-ITP"),]$new_cluster = "l-NSC_ITP"
colnames(meta_summary) = c("cluster","gene",colnames(meta_summary)[-c(1:2)])
write.csv(meta_summary, file = paste0(PATH_input,"tmp/","scaled_data_summary_NSC.csv"),row.names = F)
#-------------------------------------------------------------------------------