#----Figure_07_E_S16------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the code for Figure 07 Panel E and Figure S16 for the connectivity 
# analysis for Reinhard and Fukuda et al. 2024
# swc_L2 files have a resolution of: 0.5 nodes/micron and were used for nBLAST scores
#-------------------------------------------------------------------------------
# set variables:----------------------------------------------------------------
l_LNv = "#444444"
s_LNv = "#000000"
LN_ITP = "#BD0023"
LNd = "#FE7E00"
LPN = "#838383"
DN3 = "#009817"
DN2 = "#00B4FF"
DN1p= "#003BBD"
DN1a = "#BD00B0"
ipsi_gap = 2
contra_gap = 10
grid.col = NULL
grid.col = c("left_s-LNv" = s_LNv,"left_l-LNv" = l_LNv,"left_5th-LNv" = LN_ITP,
             "left_LNd" = LNd,"left_LPN" = LPN,"left_DN3" = DN3,"left_DN2" = DN2,
             "left_DN1p" = DN1p,"left_DN1a" = DN1a,
             "right_DN1a" = DN1a,"right_DN1p" = DN1p_A,"right_DN2" = DN2,
             "right_DN3" = DN3,"right_LPN" = LPN,"right_LNd" = LNd,
             "right_5th-LNv" = LN_ITP,"right_l-LNv" = l_LNv,"right_s-LNv" = s_LNv)

clk = read_delim(paste0(PATH_input,"clk_neurons_v",v,".csv"),
                 col_types = cols(clk_id = col_character()),delim = ",")
if (v == "630") {
  clk = clk[!clk$clk_id == "720575940632004906",] #exclude the DN1a for which no skeleton was generated
}

clk[clk$clk_name %in%c("DN1pA","DN1pB","DN1pC","DN1pD","DN1pE"),]$clk_name = "DN1p"
clk[clk$clk_name %in% c("s-CPDN3A","s-CPDN3B","s-CPDN3C","s-CPDN3D","s-CPDN3E",
                        "l-CPDN3","APDN3"),] $clk_name = "DN3"
clk[clk$clk_name %in% c("LNd_CRYn", "LNd_CRYp"),] $clk_name = "LNd"
clk[clk$clk_id %in% flywire_latestid(rootid = c("720575940634984800",
            "720575940615233954"),version = v),] $clk_name = "LNd"
clk[clk$clk_name %in% c("LN_ITP"),] $clk_name = "5th-LNv"

NSC = read_delim(paste0(PATH_input,"NSC_v",v,".csv"),
                 col_types = cols(NSC_id = col_character()),delim = ",")
NSC_DILP = NSC[NSC$NSC_id %in% flywire_latestid(rootid = c("720575940625379859", 
     "720575940612923390", "720575940624064295",
     "720575940631884883", "720575940620932045",
     "720575940628199802", "720575940643539566", 
     "720575940618694827", "720575940611254681",
     "720575940615911380", "720575940650527222",
     "720575940603765280", "720575940620628957",
     "720575940622897639", "720575940614623455",
     "720575940628363820", "720575940623081400",
     "720575940632436307"), version = v),]

if (!paste0("clk_neuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Preparations/Figure_01_E_F_preparation.py")
  source("./scripts/01_setup.R")
}
if (!paste0("neuropeptide_scale_data_summary_single_cell.csv") %in% input_files) {
  source("./scripts/Preparations/Figure_06_07_preparations.R")
  source("./scripts/01_setup.R")
}
synapses_curated = read_csv(paste0(PATH_input,"tmp/clk_neuron_input_filtered_v",v,".csv"),
          col_types =  cols(pre_pt_supervoxel_id = col_character(),
          pre_pt_root_id = col_character(),
          post_pt_supervoxel_id = col_character(),
          post_pt_root_id = col_character()))
#-------------------------------------------------------------------------------
#layer one: synaptic connections------------------------------------------------
synapses_sum = synapses_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = length(post_pt_root_id))
synapses_sum_clk = synapses_sum[synapses_sum$n_synapses >= 5 & 
                synapses_sum$pre_pt_root_id %in% clk$clk_id,]

clk_join = clk
colnames(clk_join) = c("name_pre","pre_pt_root_id","hemisphere_pre")
synapses_sum_clk = left_join(synapses_sum_clk, clk_join, by = "pre_pt_root_id")

colnames(clk_join) = c("name_post","post_pt_root_id","hemisphere_post")
synapses_sum_clk = left_join(synapses_sum_clk, clk_join, by = "post_pt_root_id")

synapses_sum_clk_grouped = synapses_sum_clk %>%
  group_by(pre_pt_root_id,name_pre,
           hemisphere_pre,post_pt_root_id,
           name_post,hemisphere_post)%>%
  summarise(n_synapses_sum = sum(n_synapses,na.rm = T))
synapses_sum_clk_grouped[synapses_sum_clk_grouped$n_synapses_sum < 5,]$n_synapses_sum = 0
synapses_sum_clk_grouped[synapses_sum_clk_grouped$n_synapses_sum >= 5,]$n_synapses_sum = 1

synapses_sum_clk_grouped = synapses_sum_clk_grouped %>%
  group_by(name_pre,hemisphere_pre,name_post,hemisphere_post)%>%
  summarise(n_synapses_sum = sum(n_synapses_sum,na.rm = T))

synapses_sum_clk_grouped[synapses_sum_clk_grouped$n_synapses_sum >= 1,]$n_synapses_sum = 1
syn_con = synapses_sum_clk_grouped
syn_con$name_pre_hemisphere = paste(syn_con$hemisphere_pre,
syn_con$name_pre,sep = "_")
syn_con$name_post_hemisphere = paste(syn_con$hemisphere_post,
 syn_con$name_post,sep = "_")
syn_con$n_synapses_sum = as.numeric(syn_con$n_synapses_sum)

syn_con_empty = as.data.frame(matrix(nrow = 18, ncol = 8))
colnames(syn_con_empty) = colnames(syn_con)
syn_con_empty$name_pre = rep(c("s-LNv","l-LNv","5th-LNv","LNd","LPN", "DN3",
                               "DN2","DN1p","DN1a"),times = 2)
syn_con_empty$hemisphere_pre = "right"
syn_con_empty$name_post = rep(c("s-LNv","l-LNv","5th-LNv","LNd","LPN", "DN3",
                                "DN2","DN1p","DN1a"),times = 2)
syn_con_empty$hemisphere_post = rep(c("right", "left"),each = 9)
syn_con_empty$n_synapses_sum = NA
syn_con_empty$name_post_hemisphere = paste(syn_con_empty$hemisphere_post,syn_con_empty$name_post, sep = "_")
syn_con = syn_con[syn_con$hemisphere_pre == "right",]
syn_con = rbind(syn_con, syn_con_empty[!syn_con_empty$name_pre %in% syn_con$name_pre |
                                         !syn_con_empty$name_post_hemisphere %in% 
                                         syn_con$name_post_hemisphere,])
syn_con$name_pre = factor(syn_con$name_pre, levels = c("s-LNv","l-LNv","5th-LNv",
                                                       "LNd","LPN", "DN3","DN2",
                                                       "DN1p","DN1a"))
syn_con$name_post_hemisphere = factor(syn_con$name_post_hemisphere,
                                      levels = c("left_s-LNv","left_l-LNv",
                                                 "left_5th-LNv","left_LNd",
                                                 "left_LPN", "left_DN3","left_DN2",
                                                 "left_DN1p","left_DN1a",
                                                 "right_s-LNv","right_l-LNv",
                                                 "right_5th-LNv","right_LNd",
                                                 "right_LPN", "right_DN3",
                                                 "right_DN2","right_DN1p","right_DN1a"))

ggplot(syn_con,aes(x=name_pre,y=name_post_hemisphere,col=n_synapses_sum, size = 5))+
  geom_point()+
  scale_color_gradient(low = "black", high = "black", na.value = "white")+
  theme(axis.text.x = element_text( angle = 90,hjust = 0,vjust = 0.5, colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 14),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_07/Figure_07_E_Syn_con.pdf"),width = 10,
       height = 10,units = "cm")
#-------------------------------------------------------------------------------
# layer two: proximity of neurons ----------------------------------------------
# based on Nagy et al. 2019 (https://doi.org/10.1371/journal.pgen.1008158) s-LNv
# signal to m-NSC_DILP via PDF. The minimal distance between s-LNv and m-NSC_DILP 
# is used as threshold for proximity calculation. If neurons are close enough to 
# each other to be within the range of this threshold, they are considered close 
# enough for paracrine signaling. 

if (!paste0("prox_con_v",v,".csv") %in% input_files) {
  if (!paste0("distance_for_paracrine_signaling_um_v",v,".csv") %in% input_files) {
    if (!paste0("swc_v",v,"_L2") %in% input_files) {
      stop("swc files in input missing! Please download the swc_L2 files
            for the root ids from codex:\n
           https://codex.flywire.ai/api/download\n
           and save them in ./input/")
    }else{ 
      neuron_list_clk = as.data.frame(matrix(ncol = 9, nrow = 0))
      for (i in c("right","left")) {
        root_id_list_clk = clk$clk_id[clk$clk_name == "s-LNv" & clk$hemisphere == i]
        for (j in root_id_list_clk) {
          clk_neuron = read.neuron.swc(paste0(PATH_input,"swc_v",v,"_L2","/",j,".swc"))
          clk_neuron_d = clk_neuron$d
          clk_neuron_d$root_id = j
          clk_neuron_d$group = "clk"
          colnames(neuron_list_clk) = colnames(clk_neuron_d)
          neuron_list_clk = rbind(neuron_list_clk,clk_neuron_d)
        }
      } 
      root_id_list_NSC = NSC_DILP$NSC_id
      neuron_list_NSC = as.data.frame(matrix(ncol =9, nrow = 0))
      for (j in root_id_list_NSC) {
        NSC_neuron = read.neuron.swc(paste0(PATH_input,"swc_v",v,"_L2","/",j,".swc"))
        NSC_neuron_d = NSC_neuron$d
        NSC_neuron_d$root_id = j
        NSC_neuron_d$group = "NSC"
        colnames(neuron_list_NSC) = colnames(NSC_neuron_d)
        neuron_list_NSC = rbind(neuron_list_NSC,NSC_neuron_d)
      }
      neuron_list = rbind(neuron_list_clk, neuron_list_NSC)
      neuron_list$PointNo = c(1:length(neuron_list$PointNo))
      
      write.csv(neuron_list,paste0(PATH_input,"tmp/tmp_neuronlist_v",v,".csv"))
    }
    reticulate::source_python("./scripts/Figure_07_S16/Figure_07_E_minimal_distance_sLNv_DILP.py")
  source("./scripts/01_setup.R")
  }
  prox_con = as.data.frame(matrix(nrow = 0,ncol=5))
  colnames(prox_con) = c("from","hemisphere_from","to","hemisphere_to","dist")
    if (!paste0("swc_v",v,"_L2") %in% input_files) {
      stop("swc files in input missing! Please download the swc_L2 files
             for the root ids from codex:\n
             https://codex.flywire.ai/api/download\n
             and save them in ./input/")
    }else{
      
      for (i in unique(clk$clk_name)) {
        root_id_list_1 = clk$clk_id[clk$clk_name == i & clk$hemisphere == "right"]
        neuron_list_1 =as.data.frame(matrix(ncol = 9,nrow = 0))
        for (j in root_id_list_1) {
          neuron_1_tmp = read.neuron.swc(paste0(PATH_input,"swc_v",v,"_L2",
                                                "/",j,".swc"))
          neuron_1_tmp_d = neuron_1_tmp$d
          neuron_1_tmp_d$root_id = j
          neuron_1_tmp_d$group = "group_1"
          colnames(neuron_list_1) =colnames(neuron_1_tmp_d)
          neuron_list_1 = rbind(neuron_list_1,neuron_1_tmp_d)
        }
        for (j in unique(clk$clk_name)) {
          for (k in c("right","left")) {
            root_id_list_2 = clk$clk_id[clk$clk_name == j & clk$hemisphere == k]
            neuron_list_2 = as.data.frame(matrix(ncol = 9, nrow = 0))
            for (l in root_id_list_2) {
              neuron_2_tmp = read.neuron.swc(paste0(PATH_input,"swc_v",v,"_L2",
                                                    "/",l,".swc"))
              neuron_2_tmp_d = neuron_2_tmp$d
              neuron_2_tmp_d$root_id = l
              neuron_2_tmp_d$group = "group_2"
              colnames(neuron_list_2) =colnames(neuron_2_tmp_d)
              neuron_list_2 = rbind(neuron_list_2,neuron_2_tmp_d)
            }
            neuron_list = rbind(neuron_list_1,neuron_list_2)
            neuron_list$PointNo = c(1:length(neuron_list$PointNo))
            write.csv(neuron_list,paste0(PATH_input,"tmp/tmp_neuronlist_v",v,".csv"))
            reticulate::source_python("./scripts/Figure_07_S16/Figure_07_E_distance_filter_neurons.py")
            close_pt = read.csv(paste0(PATH_input,"tmp/tmp_close_pt_ids_v",v,".csv"))
            far_pt = read.csv(paste0(PATH_input,"tmp/tmp_far_pt_ids_v",v,".csv"))
            prox_con_tmp = as.data.frame(matrix(nrow = 1,ncol=5))
            colnames(prox_con_tmp) = c("from","hemisphere_from","to",
                                       "hemisphere_to","dist")
            # open3d()
            # plot3d(x = close_pt$X*1000, y = close_pt$Y*1000,
            # z = close_pt$Z*1000,col="magenta",add=T)
            # plot3d(x = far_pt$X*1000, y = far_pt$Y*1000,
            # z = far_pt$Z*1000,col="black",add=T)
            prox_con_tmp$from = i
            prox_con_tmp$hemisphere_from = "right"
            prox_con_tmp$to = j
            prox_con_tmp$hemisphere_to = k
            if (length(close_pt$PointNo)>0) {
              prox_con_tmp$dist = 1
            }else{
              prox_con_tmp$dist = 0
            }
            prox_con = rbind(prox_con,prox_con_tmp)
          }
        }
      }
    }
  prox_con$from = factor(prox_con$from, levels = c("s-LNv","l-LNv","5th-LNv",
                                                   "LNd","LPN", "DN3","DN2",
                                                   "DN1p","DN1a"))
  prox_con$hemi_to = paste(prox_con$hemisphere_to,prox_con$to,sep = "_")
  prox_con$hemi_to = factor(prox_con$hemi_to,
                            levels = c("left_s-LNv","left_l-LNv", "left_5th-LNv",
                                       "left_LNd", "left_LPN", "left_DN3", 
                                       "left_DN2","left_DN1p","left_DN1a",
                                       "right_s-LNv","right_l-LNv","right_5th-LNv",
                                       "right_LNd", "right_LPN", "right_DN3",
                                       "right_DN2","right_DN1p","right_DN1a"))
  write.csv(prox_con, paste0(PATH_input,"tmp/prox_con_v",v,".csv"),row.names = F)
  source("./scripts/01_setup.R")
}
prox_con = read.csv(paste0(PATH_input,"tmp/prox_con_v",v,".csv"))
ggplot(prox_con,aes(x=from,
  y=hemi_to,col=dist,size = 5))+
  geom_point()+
  scale_color_gradient(low = "white", high = "#ffe458ff", na.value = "grey")+
  theme(axis.text.x = element_text( angle = 90,hjust = 0,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 14),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_B_Prox_con.pdf"),width = 10,
       height = 10,units = "cm")
#-------------------------------------------------------------------------------
peptides = c("AstA","AstC", "CCHa1", "CNMa", "Dh31", "Dh44", "Gpb5", "Hug", "ITP",
             "Mip", "Ms", "NPF", "Pdf", "Proc", "sNPF", "Tk", "Trissin")
clk_pept_T2A = read_delim(paste0(PATH_input,"neuropeptide_T2A.csv"), 
        delim = ";", escape_double = FALSE,trim_ws = TRUE,
        show_col_types = FALSE)
clk_pept_T2A = clk_pept_T2A[,-c(3,4)]
clk_pept_ab = read_delim(paste0(PATH_input,"neuropeptide_antibody.csv"), 
       delim = ";", escape_double = FALSE, trim_ws = TRUE,
       show_col_types = FALSE)
clk_pept_ab = clk_pept_ab[,-c(3,4)]
clk_pept_sc_raw = read_delim(
  paste0(PATH_input,"tmp/neuropeptide_scale_data_summary_single_cell.csv"), 
           delim = ",", escape_double = FALSE, trim_ws = TRUE,
  show_col_types = FALSE)
clk_pept_sc = clk_pept_sc_raw[clk_pept_sc_raw$pct_expression >= 50 & 
                                clk_pept_sc_raw$avrg_expression >= 0.208 &
                                tolower(clk_pept_sc_raw$gene) %in%
                                tolower(peptides),]
clk_pept_sc_tmp = clk_pept_sc
clk_pept_sc_tmp$cluster = factor(clk_pept_sc_tmp$cluster,
                                 levels = c("s-LNv","l-LNv","LN_ITP","LNd_NPF",
                                            "LNd_Trissin","LPN","DN3","DN3_VGlut",
                                            "DN2","DN1p","DN1p_sNPF","DN1p_AstA",
                                            "DN1p_CNMa","DN1p_CNMa_AstC","DN1p_Rh7",
                                            "DN1a"))
ggplot(clk_pept_sc_tmp[!clk_pept_sc$gene %in% c("Gpb5","Mip","Tk","Hug"),],
       aes(y=cluster,x=gene))+
  geom_point(aes(col= avrg_expression, size= pct_expression),show.legend = T)+
  scale_colour_gradientn(colours = rev(c("white","red","darkred", "black")),
                         values = c(1,0.7,0.4,0))+
  scale_size_continuous(limits = c(0,100),range = c(4,7),breaks = c(50,75,100))+
  theme(axis.text.x = element_text( angle = 90,hjust = 1,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_A_Pept_sc_threshold.pdf"),
       width = 30, height = 12,units = "cm")

clk_pept_sc = droplevels.data.frame(clk_pept_sc)
clk_pept_sc[clk_pept_sc$cluster %in% c("DN1p_sNPF", "DN1p_AstA", "DN1p_CNMa",
                                       "DN1p_CNMa_AstC", "DN1p_Rh7"),]$cluster = "DN1p"
clk_pept_sc[clk_pept_sc$cluster %in% c("DN3_VGlut"),]$cluster = "DN3"
clk_pept_sc[clk_pept_sc$cluster %in% c("LNd_NPF", "LNd_Trissin"),]$cluster = "LNd"
clk_pept_sc_tmp = clk_pept_sc[clk_pept_sc$cluster %in% c("LN_ITP"),]
clk_pept_sc_tmp[clk_pept_sc_tmp$cluster %in% c("LN_ITP"),]$cluster = "LNd"
clk_pept_sc[clk_pept_sc$cluster %in% c("LN_ITP"),]$cluster = "5th-LNv"
clk_pept_sc = rbind(clk_pept_sc, clk_pept_sc_tmp)
clk_pept_sc = unique(clk_pept_sc[,-c(3,4)])
colnames(clk_pept_sc) = c("cluster", "peptide")

clk_pept =rbind(clk_pept_T2A,clk_pept_ab
     , clk_pept_sc)
clk_pept[clk_pept$peptide %in% c("Pdf"),]$peptide = "PDF"
clk_pept[clk_pept$peptide %in% c("Dh44"),]$peptide = "DH44"
clk_pept[clk_pept$peptide %in% c("Dh31"),]$peptide = "DH31"
clk_pept = clk_pept%>%
  group_by(cluster,peptide)%>%
  summarize(n_methods = length(peptide))

clk_pept$cluster = factor(clk_pept$cluster,
                          levels = c("s-LNv","l-LNv","5th-LNv",
                                     "LNd","LPN","DN3",
                                     "DN2","DN1p","DN1a"))

ggplot(clk_pept[!clk_pept$peptide %in% c("Hug","Mip","Tk","Gpb5"),],
       aes(x=peptide,y=cluster,col=n_methods,size= 5))+
  geom_point()+
  scale_color_gradient2(low = "#e4e4e4ff", high = "black",mid = "darkgrey",
                        midpoint = 2)+
  theme(axis.text.x = element_text( angle = 90,hjust = 1,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_A_Pept_overview.pdf"),
       width = 10, height = 10,units = "cm")
clk_pept = clk_pept[clk_pept$n_methods>=2,]

clk_rec_T2A = read_delim(paste0(PATH_input,"neuropeptide_receptor_T2A.csv"), 
       delim = ";", escape_double = FALSE, trim_ws = TRUE,
       show_col_types = FALSE)
clk_rec_T2A = clk_rec_T2A[,-c(4,5)]
clk_rec_imag = read_delim(
  paste0(PATH_input,"neuropeptide_receptor_functional.csv"), 
                          delim = ";", escape_double = FALSE,
                          trim_ws = TRUE, show_col_types = FALSE)
clk_rec_imag = clk_rec_imag[,-c(4)]
clk_rec_sc_raw = read_delim(
  paste0(PATH_input,"receptors_scaled_data_summary_single_cell_with_ligand.csv"), 
          delim = ";", escape_double = FALSE, trim_ws = TRUE,
          show_col_types = FALSE)
clk_rec_sc = clk_rec_sc_raw[clk_rec_sc_raw$pct_expression >= 50 &
                              clk_rec_sc_raw$avrg_expression >= 0.067 & 
                              tolower(clk_rec_sc_raw$ligand) %in%
                              tolower(peptides),]
clk_rec_sc_tmp = clk_rec_sc
clk_rec_sc_tmp$cluster = factor(clk_rec_sc_tmp$cluster,
                             levels = c("s-LNv","l-LNv","LN_ITP","LNd_NPF",
                                        "LNd_Trissin","LPN","DN3","DN3_VGlut",
                                        "DN2","DN1p","DN1p_sNPF","DN1p_AstA",
                                        "DN1p_CNMa","DN1p_CNMa_AstC","DN1p_Rh7",
                                        "DN1a"))
ggplot(clk_rec_sc_tmp,aes(y=cluster,x=gene))+
  geom_point(aes(col= avrg_expression, size= pct_expression),show.legend = T)+
  scale_colour_gradientn(colours = rev(c("white","red","darkred", "black")),
                         values = c(1,0.7,0.4,0))+
  scale_size_continuous(limits = c(0,100),range = c(4,7),breaks = c(50,75,100))+
  theme(axis.text.x = element_text( angle = 90,hjust = 1,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_A_Rec_sc_threshold.pdf"),
       width = 30, height = 13,units = "cm")
clk_rec_sc = as.data.frame(clk_rec_sc)
clk_rec_sc[clk_rec_sc$cluster %in% c("DN1p_sNPF", "DN1p_AstA", "DN1p_CNMa",
                                     "DN1p_CNMa_AstC", "DN1p_Rh7"),]$cluster = "DN1p"
clk_rec_sc[clk_rec_sc$cluster %in% c("DN3_VGlut"),]$cluster = "DN3"
clk_rec_sc[clk_rec_sc$cluster %in% c("LNd_NPF", "LNd_Trissin"),]$cluster = "LNd"
clk_rec_sc_tmp = clk_rec_sc[clk_rec_sc$cluster %in% c("LN_ITP"),]
clk_rec_sc_tmp[clk_rec_sc_tmp$cluster %in% c("LN_ITP"),]$cluster = "LNd"
clk_rec_sc[clk_rec_sc$cluster %in% c("LN_ITP"),]$cluster = "5th-LNv"
clk_rec_sc = unique(rbind(clk_rec_sc, clk_rec_sc_tmp))

clk_rec_sc = unique(clk_rec_sc[,-c(3,4)])
colnames(clk_rec_sc) = c("cluster", "receptor", "ligand")  

clk_rec =rbind(clk_rec_T2A
    , clk_rec_sc)
clk_rec[clk_rec$receptor %in% c("AstC-R2-RB","AstC-R2",
                                "AstC-R1"),]$receptor = "AstC-R1/2"
clk_rec[clk_rec$receptor %in% c("CCHa1R"),]$receptor = "CCHa1-R"
clk_rec[clk_rec$receptor %in% c("DH31-R-RC"),]$receptor = "Dh31-R"
clk_rec[clk_rec$receptor %in% c("NPFR-RA/C","NPFR-RB/D"),]$receptor = "NPFR"
clk_rec[clk_rec$receptor %in% c("sNPFR"),]$receptor = "sNPF-R"
clk_rec[clk_rec$receptor %in% c("Dh44-R1","Dh44-R2"),]$receptor = "Dh44-R1/2"

clk_rec_tmp =rbind(unique(clk_rec_T2A[-3])
               , unique(clk_rec_sc[-3]))
clk_rec_tmp[clk_rec_tmp$receptor %in% c("AstC-R2-RB",
                                        "AstC-R2"),]$receptor = "AstC-R2"
clk_rec_tmp[clk_rec_tmp$receptor %in% c("CCHa1R"),]$receptor = "CCHa1-R"
clk_rec_tmp[clk_rec_tmp$receptor %in% c("DH31-R-RC"),]$receptor = "Dh31-R"
clk_rec_tmp[clk_rec_tmp$receptor %in% c("NPFR-RA/C",
                                        "NPFR-RB/D"),]$receptor = "NPFR"
clk_rec_tmp[clk_rec_tmp$receptor %in% c("sNPFR"),]$receptor = "sNPF-R"

clk_rec_tmp = clk_rec_tmp%>%
  group_by(cluster,receptor)%>%
  summarize(n_methods = length(receptor))
clk_rec_tmp[clk_rec_tmp$receptor %in% c("Dh44-R1",
                                        "Dh44-R2"),]$receptor = "Dh44-R1/2"
clk_rec_tmp[clk_rec_tmp$receptor %in% c("AstC-R2-RB","AstC-R2",
                                        "AstC-R1"),]$receptor = "AstC-R1/2"
clk_rec_tmp = clk_rec_tmp%>%
  group_by(cluster,receptor)%>%
  summarize(n_methods = max(n_methods))
clk_rec_tmp$cl_rc = paste(clk_rec_tmp$cluster,clk_rec_tmp$receptor, sep = "_")
clk_rec_imag$ligand= NULL
clk_rec_imag$cl_rc = paste(clk_rec_imag$cluster,clk_rec_imag$receptor, sep = "_")

clk_rec_imag$n_methods = 3
clk_rec_tmp = rbind(clk_rec_tmp[!clk_rec_tmp$cl_rc %in% clk_rec_imag$cl_rc,],
                    clk_rec_imag)
clk_rec_tmp$cluster = factor(clk_rec_tmp$cluster,
                             levels = c("s-LNv","l-LNv","5th-LNv",
                                        "LNd","LPN","DN3",
                                        "DN2","DN1p","DN1a"))

ggplot(clk_rec_tmp,aes(x=receptor,y=cluster,col=n_methods,size= 5))+
  geom_point()+
  scale_color_gradient2(low = "#e4e4e4ff", high = "black",mid = "darkgrey",
                        midpoint = 2)+
  theme(axis.text.x = element_text( angle = 90,hjust = 1,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 20),
        axis.title = element_blank(),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_A_rec_overview.pdf"),width = 10,
       height = 10,units = "cm")
clk_rec$cl_rc = paste(clk_rec$cluster,clk_rec$receptor, sep = "_")
clk_rec_tmp$cl_rc = paste(clk_rec_tmp$cluster,clk_rec_tmp$receptor, sep = "_")
clk_rec$cluster=NULL
clk_rec$receptor=NULL
clk_rec = left_join(clk_rec,clk_rec_tmp)
clk_rec$cl_rc = NULL
clk_rec = clk_rec[clk_rec$n_methods>=2,] #here the method threshold is set

pept_con = as.data.frame(matrix(nrow = 0,ncol=3))
colnames(pept_con) = c("from", "to", "connection")

for (i in unique(clk$clk_name)) {
  pept_tmp = clk_pept[clk_pept$cluster == i,]$peptide
  partner_tmp = unique(clk_rec[tolower(clk_rec$ligand) %in%
                                 tolower(pept_tmp),]$cluster)
  nrow =length(partner_tmp)
  pept_con_tmp = as.data.frame(matrix(nrow = nrow,ncol=3))
  colnames(pept_con_tmp) = c("from", "to", "connection")
  pept_con_tmp$to = partner_tmp
  pept_con_tmp$from = i
  pept_con_tmp$connection = 1
  
  pept_con = rbind(pept_con, pept_con_tmp)
}

pept_con_left = pept_con
pept_con_left$hemisphere_from ="left"
pept_con_left$hemisphere_to ="left"
pept_con_left_II = pept_con_left
pept_con_left_II$hemisphere_to ="right"
pept_con_left = rbind(pept_con_left,pept_con_left_II)

pept_con_right = pept_con
pept_con_right$hemisphere_from ="right"
pept_con_right$hemisphere_to ="right"
pept_con_right_II = pept_con_right
pept_con_right_II$hemisphere_to ="left"
pept_con_right = rbind(pept_con_right,pept_con_right_II)

pept_con = rbind(pept_con_right,pept_con_left)

pept_con$from = factor(pept_con$from,
                       levels = c("s-LNv","l-LNv","5th-LNv","LNd", "LPN", "DN3",
                                  "DN2","DN1p","DN1a"))
pept_con$hemi_to = paste(pept_con$hemisphere_to,pept_con$to,sep = "_")
pept_con$hemi_to = factor(pept_con$hemi_to,
                          levels = c("left_s-LNv","left_l-LNv", "left_5th-LNv",
                                     "left_LNd", "left_LPN", "left_DN3", 
                                     "left_DN2","left_DN1p", "left_DN1a",
                                     "right_s-LNv","right_l-LNv", 
                                     "right_5th-LNv","right_LNd", "right_LPN",
                                     "right_DN3", "right_DN2","right_DN1p",
                                     "right_DN1a"))

plot = ggplot(pept_con[pept_con$hemisphere_from == "right",],
       aes(x=from,
           y=hemi_to,size=5,col = connection))+
  geom_point()+
  scale_color_gradient(low = "#74a8ffff", high = "#74a8ffff", na.value = "grey")+
  theme(axis.text.x = element_text( angle = 90,hjust = 0,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 14),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_A_Pept_con_2_method.pdf"),
       width = 10, height = 10,units = "cm")
#-------------------------------------------------------------------------------
#combining distance and peptide layer:------------------------------------------
colnames(prox_con) = c("from","hemisphere_from","to",
                       "hemisphere_to","connection" ,"hemi_to")
paracrine_con = rbind(prox_con,pept_con)
paracrine_con = paracrine_con %>%
  group_by(from, hemisphere_from,to,hemisphere_to,hemi_to)%>%
  summarize(connection = sum(connection))
paracrine_con$connection = ifelse(paracrine_con$connection == 2, yes= 1,no = 0)

ggplot(paracrine_con[paracrine_con$hemisphere_from == "right",],
       aes(x=from, size =5,
           y=hemi_to, col= connection))+
  geom_point()+
  scale_color_gradient(low = "white", high = "#00ac0eff", na.value = "grey")+
  theme(axis.text.x = element_text( angle = 90,hjust = 0,vjust = 0.5,
                                    colour = "black"),
        axis.text.y = element_text( angle = 0,hjust = 0,colour = "black"),
        text = element_text(size = 14),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black")
  )
ggsave(paste0(PATH_output,"Figure_S16/Figure_S16_C_Paracrine_con.pdf"),width = 10,
       height = 10,units = "cm")
paracrine_con_right = paracrine_con[paracrine_con$hemisphere_from == "right",]
paracrine_con_left = paracrine_con_right
paracrine_con_left$hemisphere_from = "left"
paracrine_con_left$hemisphere_to = ifelse(paracrine_con_left$hemisphere_to == "right",
                                          yes = "left", "right")
paracrine_con = rbind( paracrine_con_right,paracrine_con_left)
paracrine_con$hemi_to = paste(paracrine_con$hemisphere_to,
                              paracrine_con$to, sep ="_")
paracrine_con$hemi_from = paste(paracrine_con$hemisphere_from,
                                paracrine_con$from, sep ="_")
col = as.data.frame(grid.col)
col$names = row.names(col)
colnames(col) = c("color","hemi_from")
paracrine_con = left_join(paracrine_con,col,by = "hemi_from")
col_tmp = paracrine_con$color

for (i in unique(clk$clk_name)) {
  
  paracrine_con$color = col_tmp
  paracrine_con[paracrine_con$from != i,]$color = adjustcolor(paracrine_con[paracrine_con$from != i,]$color,alpha.f = 0)
  
  circos.par(start.degree = -96.5,
             gap.after = c(
               "left_s-LNv" = ipsi_gap,"left_l-LNv" = ipsi_gap,
               "left_5th-LNv" = ipsi_gap,"left_LNd" = ipsi_gap,
               "left_LPN" = ipsi_gap,"left_DN3" = ipsi_gap,
               "left_DN2" = ipsi_gap,"left_DN1p" = ipsi_gap,
               "left_DN1a" = contra_gap,
               
               "right_DN1a" = ipsi_gap,"right_DN1p" = ipsi_gap,
               "right_DN2" = ipsi_gap,"right_DN3" = ipsi_gap,
               "right_LPN" = ipsi_gap,"right_LNd" = ipsi_gap,
               "right_5th-LNv" = ipsi_gap,
               "right_l-LNv" = ipsi_gap,
               "right_s-LNv" = contra_gap
             )
  )
  pdf(paste0(PATH_output,"Figure_07/","Figure_07_E_",i,
             "_paracrine_connectivity_diagram_v",v,".pdf"),width = 5, height = 5)
  chordDiagram(paracrine_con[,c(7,5,6)], scale = F, link.lwd = 5, 
               link.border  = paracrine_con$color,
               order = c("left_s-LNv","left_l-LNv",
                         "left_5th-LNv","left_LNd","left_LPN",
                         "left_DN3","left_DN2","left_DN1p","left_DN1a",
                         "right_DN1a","right_DN1p","right_DN2",
                         "right_DN3","right_LPN","right_LNd",
                         "right_5th-LNv","right_l-LNv","right_s-LNv"),
               col = paracrine_con$color,
               grid.col = grid.col,
               link.visible = paracrine_con$hemisphere_from == "right",
               directional = 1,
               direction.type = c( "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = c("grid"),
               annotationTrackHeight = mm_h(7),
               abline(v = 0, lty = 2, col = "#000000",lwd = 3),
               preAllocateTracks = list()
  )
  circos.clear()
  dev.off()  
}
