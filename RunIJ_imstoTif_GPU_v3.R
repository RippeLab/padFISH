# This script is used to run ImageJ macro to convert .ims files to .tif on a selected folder
# All the output of this script will be stored in the input folder
# Output: converted Tif images, bash script, ImageJ log file (IJ_imstoTif.ocurryjobnumber) and ImageJ error file (IJ_imstoTif.ecurryjobnumber)

## Step 1 ##
# Copy the macro to your LSDF2 (database) folder and change the path of the macro

## Step 2 ##
# Extract a folder list containing the files that should be processed by ImageJ
# Choose if you want to save z-stack (unprojected) .tif images and maximum projected .tif images

## Step 3 ##
# A "imstoTif R function" is defined, to which a single folder path can be handed over
# The function creates a bash script and submits it to the cluster via the qsub system
# The bash script executes ImageJ and run the macro on the specified folder 
# The bash script creates a copy of itself and stores it in the folder

# Define where the macro is located on LSDF2 (database)
macro <- "/path/to/FIJImacro/imstoTif_headless_v2.ijm"

# Define input folders for images (without "/" at the end)
folder <- "/output/folder"
setwd(folder)
folderlist <- list.dirs(folder)
# inspect folderlist to later select subfolders
folderlist

# "true" if you want to save z-stack (unprojected) tif files and/or maximum projected tif files
# use "false" if you do not want to save 
zstack <- "false"
maxproj <- "true"
# this is specifically for dSTORM to lower the image size
# for resolution one can choose either "1024x1024" or "2048x2048", this depends on the settings you used during imaging
crop <- "false"
resolution <- "1024x1024"

## R imstoTif function (create bash script for cluster submission) ##
# Adjust the walltime and memory required for the job, this depends on the size and number of your images
# If you request too much resources it will likely take longer for Currycluster to schedule your job

imstoTif <- function(folder){
  sh_script<-c('#!/bin/bash',
               '#PBS -N IJ_imstoTif',
               '#PBS -l walltime=2:00:00',
               '#PBS -l nodes=1:ppn=8',
               '#PBS -l mem=50gb',
               '#PBS -V',
               '#PBS -q gpu',
               'module load fiji',
               paste0('ImageJ-linux64 --headless --run ',macro, ' \'dir=\"',folder,  '\", saveZstack=',zstack, ', saveMaxproj=',maxproj, ', saveCrop=',crop, ', imres=\"',resolution,'\"\''),
               'exit')
  write.table(sh_script,paste0(folder,'/imstoTifScript.sh'),row.names = F,quote = F,col.names = F)
  
  system(paste0('qsub ',folder,'/imstoTifScript.sh' ), wait = F)
}

# Run imstoTif function on the folder directly
#imstoTif(folder)

#OR:

# Run on a subfolder (select from folderlist above)
imstoTif(folder = paste0(folderlist[3]))
imstoTif(folder = paste0(folderlist[4]))
imstoTif(folder = paste0(folderlist[5]))
