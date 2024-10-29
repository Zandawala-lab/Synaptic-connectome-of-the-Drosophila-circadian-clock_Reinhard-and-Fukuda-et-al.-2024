#----Figure_06_G----------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 06 Panel G for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
if (!paste0("scaled_data_summary_NSC.csv") %in% input_files) {
  source("./scripts/Preparations/Figure_06_07_preparations.R")
  source("./scripts/01_setup.R")
}
NSC_rec_sc = read_delim(paste0(PATH_input,"tmp/scaled_data_summary_NSC.csv"), 
                        delim = ",", escape_double = FALSE,
                        trim_ws = TRUE, show_col_types = FALSE)
#-------------------------------------------------------------------------------
NSC_rec_sc$cluster = factor(NSC_rec_sc$cluster,
                      levels = rev(c("l-NSC_DH31","l-NSC_CRZ","l-NSC_ITP",
                                   "m-NSC_DILP","m-NSC_DH44")))
NSC_rec_sc$gene = factor(NSC_rec_sc$gene,
                    levels = c("Crz","sNPF","Dh44","CG13248","CG13743",
                                "Ilp2","Ilp3","Ilp5","Tk","ITP","Dh31","ImpL2",
                                "AstA-R1","AstA-R2","AstC-R1","AstC-R2",
                                "CCHa1-R","CNMaR","Dh31-R","Dh44-R1","Dh44-R2",
                                "NPFR","Pdfr","Proc-R", "sNPF-R","TrissinR"
                                        ))

ggplot(NSC_rec_sc,aes(y=cluster,x=gene))+
  geom_point(aes(col= avrg_expression, size= pct_expression),show.legend = T)+
  scale_colour_gradientn(colours = rev(c("grey","white","red","darkred", "black")),
                         values = c(1,0.7,0.4,0.2,0))+
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
ggsave(paste0(PATH_output,"Figure_06/Figure_06G.pdf"),width=14,height =5)


