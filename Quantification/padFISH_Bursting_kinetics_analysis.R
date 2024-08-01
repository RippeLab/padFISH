### Analysis for Bursting kinetics - calculate Integrated Intensity in cell nuclei from padFISH  #####
# This script quantifies the intensities of target genes (BIRC2, NFKBIA, SELE) to assess their transcription bursting kinetics ###
# It uses the "quantNuclei.R" custom function to compute the sum of intensity in the pre-segmented cellpose nuclear mask
# Thresholding of each channel is applied to define active or inactive transcriptional burst is set from the intensity distribution  and visual inspection

# Input data (stitched images):
#1 cellpose nuclear masks segmented from DAPI channel (_corr_cp_masks.png)
#2 (optional) ROIs masks to exclude from analysis (e.g. laser overexposure), pre-selected and binarized in FIJI (_ROIs_masks_inverted.png)
#3 channel of interest in TIF format for target gene (Channel.tif)


# Load required packages and functions
library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
source("/filepath/quantNuclei_v01.R")

# Assign working directory ("/"in the end)
input_dir <- "/filepath/R_input/"
# Retrieve files in the directory 
all_files <- list.files(input_dir, recursive = TRUE)
# Select files for nuclear masks and masks for ROIs (areas to exclude)
nucmasks_image <- all_files[grep("_corr_cp_masks.png$",all_files)]
nucmasks_ROIs <- all_files[grep("_ROIs_masks_inverted.png$",all_files)]
# Load images
nucmasks_ori <- readImage(paste0(input_dir, nucmasks_image))
nucmasks_exclude <- readImage(paste0(input_dir, nucmasks_ROIs))

# Reformat nuclear masks (only for cellpose mask)
nucmasks_ori <- Image(as.numeric(as.factor(nucmasks_ori))-1,dim = dim(nucmasks_ori))
# Look at nuclei masks image 
display(nucmasks_ori)
display(nucmasks_exclude)

##### Filtering steps for the nuclear masks ####

# Filtering Step 1:  Exclude regions of interest of the nuclear mask for the analysis (overexposure, overlap of cells...)

# OPTIONAL: Invert nucmasks_exclude (if not already done in FIJI)
#invert_nucmasks_excl <- 1 - nucmasks_exclude
#display(invert_nucmasks_excl)
#nucmasks_exclude <- invert_nucmasks_excl

# If mask is already inverted (ROIs = 0), continue from here:
# Multiply the nuclear mask with the inverted exclude mask
nucmasks_filt1 <- nucmasks_ori * nucmasks_exclude
display(nucmasks_filt1)
# Save the final mask and re-assign it to nucmasks fo further steps
writeImage(nucmasks_filt1, "/filepath/R_output/Xmin_nucmasks_filt1.png")
nucmasks <- nucmasks_filt1

# Filtering Step 2:  Remove nuclei that are cut at the edges of the image

#set the border size on which to exclude nuclei (10-20)
dims <- dim(nucmasks)
border_size <- 10
top_row <- nucmasks[1:border_size,]
left_column <- nucmasks[,1:border_size]
bottom_row <- nucmasks[(dims[1]-border_size):dims[1],]
right_column <- nucmasks[,(dims[2]-border_size):dims[2]]
#assign border
border <- c(top_row,left_column,bottom_row,right_column)
ids <- unique(border[which(border != 0)])
#apply rmObjects (works on integer)
storage.mode(nucmasks) <- 'integer'
nucmasks_filt2 <- rmObjects((nucmasks), ids)
storage.mode(nucmasks_filt2) <- 'double'
#display image to check result 
display(nucmasks_filt2)
# If nuclei at borders are still included, repeat from: nucmasks <- nucmasks_filt1 and change border_size

# Save filtered image (include "/2^16 max pixel value so that the nuclei IDs will be displayed in FIJI), include +1 to start cell count from ID=1
writeImage((nucmasks_filt2+1)/2^16, "/filepath/R_output/Xmin_nucmasks_final.tiff",type="tiff",bits.per.sample = 16)


###### Quantify Sum of Intensity across cell nuclei #######

#print initial number of cells used in the analysis (save for statistics)
Ncell_initial <- length(unique(as.vector(nucmasks_filt2)))
print(Ncell_initial)
# Assign filtered image to nuclei_masks before proceeding 
nuclei_masks <- nucmasks_filt2

####### ASSIGN TARGET GENE for QUANTIFICATION STEP (Change file names and output accordingly) #######
# Assign gene/channel image for analysis
gene_image <- all_files[grep("_Channel_Xnm_corr.tif$",all_files)]
gene_channel <- readImage(paste0(input_dir, gene_image))

#Define function quantNuclei to calculate integrated intensity (sum) within nuclei
gene_IntegrInt <- quantNuclei(nuclei_masks, gene_channel, sum) %>% data.frame() %>% setNames('IntegratedIntensity')
#View(gene_IntegrInt)

# create data frame with x and y coordinates of nuclear masks (and nuclear area) 
nucfeat <- cbind(computeFeatures.moment(nuclei_masks),computeFeatures.shape(nuclei_masks))
nucfeat <- data.frame("x"=round(nucfeat[,"m.cx"]),"y"=round(nucfeat[,"m.cy"]),"nucleus_area"=nucfeat[,"s.area"])
nucfeat$maskID <- as.integer(rownames(nucfeat))
gene_IntegrInt$maskID <- as.integer(rownames(gene_IntegrInt))
# bind nucfeat and gene_IntegrInt by the maskID
fulldf <- full_join(nucfeat, gene_IntegrInt, by="maskID")
# Save resulting data frames (change gene name accordingly)
write.csv(gene_IntegrInt, "/filepath/R_output/gene_Xmin_Sum_IntegrInt_table.csv", row.names = TRUE)
write.csv(fulldf, "/filepath/R_output/gene_Xmin_Sum_IntegrInt_XY_full_table.csv", row.names = TRUE)


#### Plot Integrated Intensity distribution ####

# Read the integrated intensity data from the saved CSV file or assign 'geneIntegrInt' in the environment to int_table)
int_table <- gene_IntegrInt
#int_table <- read.csv("/filepath/R_output/gene_Xmin_Sum_IntegrInt_table.csv")
# Inspect intensity table and check x-limit and bin width
hist(int_table$'IntegratedIntensity')
# Save intensity plot
out_folder <- "/filepath/R_output/gene_Xmin_Sum_IntegrInt_plot.pdf"
pdf(out_folder)
hist(int_table$'IntegratedIntensity')
dev.off()

# Plot histogram (select bin-width and x-lim values accordingly for better visualization)
ggplot(int_table, aes(x = `IntegratedIntensity`)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  xlim(0, 100) +  # Set x-axis limits
  labs(title = "Integrated Intensity Distribution",
       x = "Sum intensity",
       y = "N counts (cells)") +
  theme_minimal()

# Define final x axis limits and binwidth, then save plot in the output folder
SumInt_Plot <- ggplot(int_table, aes(x = `IntegratedIntensity`)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  xlim(0, 100) +  # Set x-axis limits
  labs(title = "Integrated Intensity Distribution",
       x = "Sum intensity",
       y = "N counts (cells)") +
  theme_minimal()
ggsave(filename = "/filepath/R_output/gene_Xmin_Sum_IntegrInt_plot_bin.pdf", plot = SumInt_Plot)



###### Inspect cells to define intensity threshold for the transcriptional burst ######
# (After files have been saved in output folder) #
library(EBImage)

# retrieve and display the final mask used for the analysis fro the R output folder
test_input <- "/filepath/R_output/"
test_files <- list.files(test_input, recursive = TRUE)
test_masks <- test_files[grep("_nucmasks_final.tiff$",test_files)]
final_masks <- readImage(paste0(test_input, test_masks))
final_nucmasks <- Image(as.numeric(as.factor(final_masks))-1,dim = dim(final_masks))
display(final_nucmasks)

# check singular cell ID position --> cell ID = row number, check in the csv table which ID corresponds to the FIJI mask.
cell_ID <- 3    # insert cell ID of interest
cell_mask <- final_nucmasks == cell_ID
display(cell_mask)

#OR: if mask is already in the environment (nucmasks_filt2)
cell_ID <- 3    # insert cell ID of interest
cell_mask <- nucmasks_filt2 == cell_ID
display(cell_mask)

