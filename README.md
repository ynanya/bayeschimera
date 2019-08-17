# bayeschimera
NGS-based estimation of mixed chimera status after allogeneic stem cell transplantation.
This is beta-version. 


#  Basic Concepts
This program estimates the mixed chimerism for the bam file from blood or marrow cells after allogeneic stem cell transplatation. 


#  

## BAM files
This program takes two arguments: PATH_TO_RECIPIENT_BAM_FILE and PATH_TO_CHIMERA_BAM_FILE.

PATH_TO_RECIPIENT_BAM_FILE: path to bam file of the DNA where 100% recipient chimerism are assumed. The tumor samples before transplant or normal control (oral swab, skin fibroblast, non-tumor-containing T cells,...) 
PATH_TO_CHIMERA_BAM_FILE: path to bam file of blood or marrow cells after allogeneic transplantation. 


Both of these bam files should include data for SNPs  


## Lsit of SNPs
The list of SNPs to be utilized for calculation should be provided in data folder.
The default name of this list is 'snps_list.txt', but you can specifically indicate the path to SNPs list file as the command option. 
An example list is provided in the data folder, and modify the data following this format. 


#  Software Requirements
R version 3.6 or later
R packages: rstan, ggmcmc, tidyverse

Execute the following command once to install these packages. 

if(!require("rstan")){install.packages("rstan")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("ggmcmc")){install.packages("ggmcmc")}

#  Usage of the command
chimerism.sh [OPTIONS] PATH_TO_RECIPIENT_BAM_FILE  PATH_TO_CHIMERA_BAM_FILE

OPTIONS:
[-d dir] : path to SNPs list file.
[-h]: show this message.


# Environment
This program was developed under CentOS 7, R version 3.6, PHP 7.2.13.
