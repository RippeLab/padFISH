# This script is used to run Cellpose on selected images or lists of image files (batch mode)
# RUN ON RSTUDIO SERVER from CURRY CLUSTER

## Step 1 ##
# Extract a folder list containing the files that should be processed by Cellpose
# Run Cellpose on specific files located in these folders (e.g. only DAPI channel) 
# Use the --imgfilter parameter to only run cellpose on file names that contain a user-defined string

## Step 2 ##
# A "cellpose R function" is defined, to which a single folder path can be handed over.
# The function created a bash script and submits it to the cluster via the qsub system
# The bash script executes Cellpose on the specified folder and/or files (see --imgfilter)
# The bash script creates a copy of itself and stores it in the inputfolder = outputfolder

# Define input folders for Cellpose (without "/" at the end)
folder <- "/input/folder/path"
setwd(folder)
folderlist <- list.dirs(folder)
folderlist

# Optional:
#folderfiles <- list.files(folder, recursive = TRUE)
#folderfiles

# Subset folderlist further (Optional), such as _Fused folder obtained from pre-processing (with stitched images)
folderlist <- folderlist[grep("_Fused",folderlist)]
folderlist

## R cellpose function (create bash script for cluster submission) ##
# Adjust --img_filter to select specific input files (_Fused_405nm_corr.tif = DAPI images for segmentation)
# --diameter depends on your image size/size of nuclei (HUVECs: 60x objective = 150, 100x objective = 200)
# The pretrained_model cyto works best for most nuclei/images

cellpose<-function(folder){
  sh_script<-c('#!/bin/sh -x',
               '#PBS -N Rcellpose',
               '#PBS -l walltime=1:00:00',
               '#PBS -l nodes=1:ppn=8',
               '#PBS -l mem=32gb',
               '#PBS -V',
               '#PBS -q gpu',
               'module load conda ',
               'module load cuda ',
               'conda activate cellpose-cuda',
               paste0('cellpose --use_gpu --dir ',folder,' --img_filter "_Fused_405nm_corr" --pretrained_model cyto --diameter 150 --save_png --no_npy --verbose'),
               'exit',
               'exit')
  write.table(sh_script,paste0(folder,'cellposeScript.sh'),row.names = F,quote = F,col.names = F)
  
  system(paste0('qsub ',folder,'cellposeScript.sh' ),wait = F)
}

##### Run Cellpose on selected folders/files (check folderlist and select) ####
# Note: Cellpose requires a lot of computational resources
# Do not submit too many jobs in parallel to leave space on the cluster.

# Run cellpose on N folders in parallel (e.g. 1-2-3)
cellpose(folder = paste0(folderlist[1],"/"))
cellpose(folder = paste0(folderlist[2],"/"))
cellpose(folder = paste0(folderlist[3],"/"))
