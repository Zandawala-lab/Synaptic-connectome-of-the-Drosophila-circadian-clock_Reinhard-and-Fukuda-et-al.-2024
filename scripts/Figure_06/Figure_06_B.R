#----Figure_06_B----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 06 Panel B for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("neuropeptide_scale_data_summary_single_cell.csv") %in% input_files) {
  source("./scripts/Preparations/Figure_06_07_preparations.R")
  source("./scripts/01_setup.R")
}

clk_pept_sc_raw = read_delim(paste0(PATH_input,
                    "tmp/neuropeptide_scale_data_summary_single_cell.csv"), 
                    delim = ",", escape_double = FALSE, trim_ws = TRUE,
                    show_col_types = FALSE)
#-------------------------------------------------------------------------------
clk_pept_sc_raw$cluster = factor(clk_pept_sc_raw$cluster,
                                 levels = c("s-LNv","l-LNv","LN_ITP","LNd_NPF",
                                            "LNd_Trissin","LPN","DN3","DN3_VGlut",
                                            "DN2","DN1p","DN1p_sNPF","DN1p_AstA",
                                            "DN1p_CNMa","DN1p_CNMa_AstC","DN1p_Rh7",
                                            "DN1a"))
clk_pept_sc_raw$gene = factor(clk_pept_sc_raw$gene,
                              levels = c("tim","per","Clk","cry","Rh7","gl","Pdf",
                                         "AstC","CNMa","ITP","Dh31","sNPF","AstA",
                                         "Dh44","NPF","Trissin","CCHa1","Proc",
                                         "Gpb5","Mip","Ms","Hug","Tk"))

ggplot(clk_pept_sc_raw,aes(y=cluster,x=gene))+
  geom_point(aes(col= avrg_expression, size= pct_expression),show.legend = T)+
  scale_colour_gradientn(colours = rev(c("grey","white","red","darkred", "black")),
                         values = c(1,0.7,0.4,0.13,0))+
  #scale_color_gradient2(low = "blue", high = "red",mid = "white")+ 
  scale_size_continuous(limits = c(0,100),range = c(2,10),breaks = c(25,50,75,100))+
  theme(axis.text.x = element_text( angle = 45,hjust = 1,vjust = 1, colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_06/Figure_06B.pdf"),width=15,height =7)

