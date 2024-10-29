#----01_setup-------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This is the setup file for the connectivity analysis for 
# Reinhard and Fukuda et al. 2024
# Here all necessary libraries are loaded and all used functions are initialized. 
# 
# To run the analyzes execute the scripts according to the Figure panel. 
#-------------------------------------------------------------------------------
#libraries: --------------------------------------------------------------------
library(natverse)
library(data.table)
library(tidyverse)
library(coconat)
library(circlize)
library(ggpattern)
library(fafbseg)
library(rglplus)
library(dendroextras)
#general variables:-------------------------------------------------------------
PATH_input = "./input/"
PATH_output = "./output/"
input_files = list.files(path = PATH_input, full.names = FALSE, recursive = FALSE)
input_files_tmp = list.files(path = paste0(PATH_input,"tmp/"),
                            full.names = FALSE, recursive = FALSE)
input_files = c(input_files,input_files_tmp)
v = read_delim(paste0(PATH_input,"version.csv"),
                     col_types = cols(version = col_character()),delim = ";")
v = v$version[1]

grey_scale = colorRampPalette(c("white","black"))
red_scale = colorRampPalette(c("white","#b20000ff"))
l_LNv = "#444444"
s_LNv = "#000000"
LN_ITP = "#BD0023"
LNd_CRYp = "#ffe200ff"
LNd_CRYn = "#FE7E00"
LPN = "#838383"
l_CPDN3 = "#8EF444"
APDN3 = "#8EFFC1"
s_CPDN3A = "#009817"
s_CPDN3B = "#bbff48ff"
s_CPDN3C = "#59ff73ff"
s_CPDN3D = "#4da55bff"
s_CPDN3E = "#518100ff"
DN2 = "#00B4FF"
DN1p_A = "#003BBD"
DN1p_B = "#8100FF"
DN1p_C = "#003BBD"
DN1p_D = "#003BBD"
DN1p_E = "#003BBD"
DN1a = "#BD00B0"

#functions:---------------------------------------------------------------------

seqlast <- function (from, to, by) 
{
 vec <- do.call(what = seq, args = list(from, to, by))
 if ( tail(vec, 1) != to ) {
 return(c(vec, to))
 } else {
 return(vec)
 }
}
