#----Figure_07_A----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 07 Panel A for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
clk_rec_sc_raw = read_delim(paste0(PATH_input,
                  "tmp/receptors_scaled_data_summary_single_cell.csv"), 
                   delim = ",", escape_double = FALSE, trim_ws = TRUE,
                   show_col_types = FALSE)
#-------------------------------------------------------------------------------
clk_rec_sc_raw$cluster = factor(clk_rec_sc_raw$cluster,
                                levels = c("s-LNv","l-LNv","LN_ITP","LNd_NPF",
                                           "LNd_Trissin","LPN","DN3","DN3_VGlut",
                                           "DN2","DN1p","DN1p_SNPF","DN1p_AstA",
                                           "DN1p_CNMa","DN1p_CNMa_AstC","DN1p_Rh7",
                                           "DN1a"))

clk_rec_sc_raw$gene = factor(clk_rec_sc_raw$gene,
                             levels = c("tim","per","Clk","cry","Rh7","gl","Pdfr",
                                        "AstC-R2","AstA-R1","Dh31-R","sNPF-R",
                                        "CNMaR","Dh44-R1","TrissinR","Proc-R",
                                        "CCHa1-R","Dh44-R2","NPFR","AstC-R1"))

ggplot(clk_rec_sc_raw[!clk_rec_sc_raw$gene %in% c("tim","per","Clk","cry","Rh7",
                                                  "gl"),],aes(y=cluster,x=gene))+
  geom_point(aes(col= avrg_expression, size= pct_expression),show.legend = T)+
  scale_colour_gradientn(colours = rev(c("grey","white","red","darkred", "black")),
                         values = c(1,0.7,0.4,0.19,0))+
  scale_size_continuous(limits = c(0,100),range = c(2,10),breaks = c(25,50,75,100))+
  theme(axis.text.x = element_text( angle = 45,hjust = 1,vjust = 1, colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "top",
        legend.box="vertical", legend.margin=margin(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_07/Figure_07A.pdf"),width=7,height =8)

