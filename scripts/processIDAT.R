#### Process IDAT files
library(tidyr)
library(dplyr)
library(illuminaio)
library(DNAcopy)
library(GWASTools)
library(GenomicRanges)
library(magrittr)
library(bracer)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <path to red idat> <path to green idat> <output table>', call.=FALSE)
}


idat_file_path_red <- args[1]
idat_file_path_green <- args[2]
red_channel <- readIDAT(idat_file_path_red)
green_channel <- readIDAT(idat_file_path_green)
red_intensity <- red_channel$Quants
green_intensity <- green_channel$Quants
total_intensity <- red_intensity + green_intensity
write.table(total_intensity, file = args[3], quote = F, sep = "\t")
