### Integrated Intensity IF #####

###############################################################################################
# Load required packages and functions
library(EBImage)
library(tidyr)
library(dplyr)
library(ggplot2)
source("/Volumes/sd17B002/Irene/00_scripts/01_Quantification/intensity/quantNuclei_v01.R")
###############################################################################################

# Assign working directory ("/"in the end!)
input_dir <- "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_input/"
# Retrieve files in the directory 
all_files <- list.files(input_dir, recursive = TRUE)
# Select files for nuclear masks (ROIs = areas to exclude from the original nuclear mask)
nucmasks_image <- all_files[grep("_corr_cp_masks.png$",all_files)]
#optional
nucmasks_ROIs <- all_files[grep("_inverted.png$",all_files)]
# Load images
nucmasks_ori <- readImage(paste0(input_dir, nucmasks_image))
nucmasks_exclude <- readImage(paste0(input_dir, nucmasks_ROIs))

# Reformat nuclear masks (necessary for cellpose mask)
nucmasks_ori <- Image(as.numeric(as.factor(nucmasks_ori))-1,dim = dim(nucmasks_ori))
#nucmasks_exclude <- Image(as.numeric(as.factor(nucmasks_exclude))-1,dim = dim(nucmasks_exclude))
# Look at nuclei masks image 
display(nucmasks_ori)
nucmasks <- nucmasks_ori

# OPTIONAL Pre-Filtering Step:  Exclude regions of interest of the nuclear mask for the analysis (overexposure, stitching errors...)
## Invert nucmasks_exclude (if not already done in FIJI)
#invert_nucmasks_excl <- 1 - nucmasks_exclude
#display(invert_nucmasks_excl)
#nucmasks_exclude <- invert_nucmasks_excl

# If mask is already inverted (ROIs = 0), continue from here:
# Multiply the nuclear mask with the inverted exclude mask
nucmasks_filt1 <- nucmasks_ori * nucmasks_exclude
display(nucmasks_filt1)
# Save the final mask and re-assign it to nucmasks fo further steps
writeImage(nucmasks_filt1, "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_output/00_240611_HUVEC_p65-Ser276_240min_nucmasks_filt1.png")
nucmasks <- nucmasks_filt1


#check if nucmasks is binarized (number of cells and pixel area will be displayed in Console)
#table(nucmasks_ori)

######################## Filtering steps for the nuclear masks ######################
# Remove nuclei that are cut at the edges of the image

dims<-dim(nucmasks)
# Set the border size on which to exclude nuclei (border of 1 pixel is NOT enough) - use 10 or more
border_size <- 20
top_row <- nucmasks[1:border_size,]
left_column <- nucmasks[,1:border_size]
bottom_row <- nucmasks[(dims[1]-border_size):dims[1],]
right_column <- nucmasks[,(dims[2]-border_size):dims[2]]
# Assign border
border <- c(top_row,left_column,bottom_row,right_column)
ids <- unique(border[which(border != 0)])
# Transform image as integer before removing objects (somehow rmObjects does not work on double anymore), then change back to double
storage.mode(nucmasks) <- 'integer' 
nucmasks_final <- rmObjects((nucmasks), ids)
storage.mode(nucmasks_final) <- 'double'

# Look at image after nuclei are removed
display(nucmasks_final)
# If nuclei at borders are still included, repeat from: nucmasks <- nucmasks_filt1 or nucmasks <- nucmasks_ori and change border_size (ex- 15)
# Save filtered image (include "/2^16 max pixel value so that the nuclei IDs will be displayed in FIJI), include +1 to start cell count from ID=1
writeImage((nucmasks_final+1)/2^16, "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_output/01_240611_HUVEC_p65-Ser276_240min_nucmasks_final.tiff",type="tiff",bits.per.sample = 16)

############################### Quantify Intensity across cell nuclei ###############################

#OPTIONAL: print Number of cells that will be used from now on in the analysis (see result in Console)
Ncell_initial <- length(unique(as.vector(nucmasks_final)))
print(Ncell_initial)

# Assign filtered image to nuclei_masks before proceeding 
nuclei_masks <- nucmasks_final

########## ASSIGN CHANNEL for QUANTIFICATION STEP (Change file names and output accordingly) #######
# Assign gene image
IF_image <- all_files[grep("_Fused_637nm_corr.tif$",all_files)]
IF_channel <- readImage(paste0(input_dir, IF_image))

####### STEP 1: Define function to calculate sum intensity #####

# Use quantNuclei to calculate integrated intensity (sum intensity) within nuclei
p65_IntegrInt <- quantNuclei(nuclei_masks, IF_channel, sum) %>% data.frame() %>% setNames('IntegratedIntensity')

# create data frame with x and y coordinates of nuclear masks (and area) 
nucfeat <- cbind(computeFeatures.moment(nuclei_masks),computeFeatures.shape(nuclei_masks))
nucfeat <- data.frame("x"=round(nucfeat[,"m.cx"]),"y"=round(nucfeat[,"m.cy"]),"nucleus_area"=nucfeat[,"s.area"])
nucfeat$maskID <- as.integer(rownames(nucfeat))
p65_IntegrInt$maskID <- as.integer(rownames(p65_IntegrInt))
# bind nucfeat and gene_IntegrInt by the maskID
fulldf <- full_join(nucfeat, p65_IntegrInt, by="maskID")

# Save resulting data frames (change gene name accordingly)
write.csv(p65_IntegrInt, "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_output/02_240611_HUVEC_p65-Ser276_240min_intensity_table.csv", row.names = TRUE)
write.csv(fulldf, "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_output/03_240611_HUVEC_p65-Ser276_240min_intensity_table_XY_full.csv", row.names = TRUE)
## Plot Integrated Intensity distribution 

# Read the integrated intensity data from the saved CSV file (or assign 'geneIntegrInt' in the environment to int_table)
int_table <- p65_IntegrInt
#int_table <- read.csv("/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_stitched_adjusted/_analysis/240min/R_output/02_240611_HUVEC_p65-Ser276_240min_intensity_table.csv")

# Inspect intensity table and check x-limit and bin width for later
hist(int_table$'IntegratedIntensity')

# Save intensity plot
out_folder <- "/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/240min/R_output/04_240611_HUVEC_p65-Ser276_240min_intensity_plot.pdf"
pdf(out_folder)
hist(int_table$'IntegratedIntensity')
dev.off()

#plot conditions together
p65_0h_inttable <- read.csv("/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65_pan/_analysis/0h/R_output/02_240611_HUVEC_p65_0h_intensity_table.csv")
p65_30min_inttable <- read.csv("/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65_pan/_analysis/30min/R_output/02_240611_HUVEC_p65_30min_intensity_table.csv")
p65_240min_inttable <- read.csv("/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65_pan/_analysis/240min/R_output/02_240611_HUVEC_p65_240min_intensity_table.csv")

p65_0h_df <- p65_0h_inttable %>% select(IntegratedIntensity)
colnames(p65_0h_df)<- "IntegratedIntensity"
p65_0h_df$condition <- "1.WT" 
p65_0h_df

p65_30min_df <- p65_30min_inttable %>% select(IntegratedIntensity)
colnames(p65_30min_df)<- "IntegratedIntensity"
p65_30min_df$condition <- "2.TNFa+30min" 
p65_30min_df

p65_240min_df <- p65_240min_inttable %>% select(IntegratedIntensity)
colnames(p65_240min_df)<- "IntegratedIntensity"
p65_240min_df$condition <- "3.TNFa+240min" 
p65_240min_df

df_conditions <- rbind(p65_0h_df, p65_30min_df, p65_240min_df)

ggplot(df_conditions, aes(x=condition, y=log(IntegratedIntensity, 10)))+geom_boxplot()+labs(y= "Integrated Intensity")
ggplot(df_conditions, aes(x=condition, y=IntegratedIntensity))+geom_boxplot()+labs(y= "Integrated Intensity")
final_plot <- ggplot(df_conditions, aes(x=condition, y=IntegratedIntensity))+geom_boxplot()+labs(y= "Integrated Intensity")
ggsave(paste0("/Volumes/sd17B002/Irene/05_IF/HUVEC/2024-06-11/p65-Ser276/_analysis/", "final_240611_HUVEC_p65-Ser276_conditions_plot.png"), final_plot, dpi=300)

ggplot(df_conditions, aes(x = condition, y = log(IntegratedIntensity, 10))) +
  geom_boxplot() +
  labs(y = "Integrated Intensity") +
  ylim(2.0, NA)

###############

ggplot(df_conditions, aes(x=condition, y=log(IntegratedIntensity, 10)))+geom_boxplot(outlier.shape = NA)+labs(y= "Nuclear Intensity (log10)")+scale_y_continuous(limits = quantile(df_conditions$"Integrated Intensity", c(0.1, 0.9)))
ggplot(df_conditions, aes(x=condition, y=log(IntegratedIntensity, 10)))+geom_boxplot(outlier.shape = NA)+labs(y= "Nuclear Intensity (log10)")+scale_y_continuous(limits = quantile(df_conditions$"Integrated Intensity", c(0.2, 0.8)))
ggplot(df_conditions, aes(x=condition, y=log(IntegratedIntensity, 10)))+geom_boxplot(outlier.shape = NA)+labs(x= "Condition", y= "Nuclear Intensity (log10)")+scale_y_continuous(limits = c(2.4, 3.6))

geom_boxplot(outlier.shape = NA)

scale_y_continuous(limits = quantile(dfr$y, c(0.1, 0.9)))








