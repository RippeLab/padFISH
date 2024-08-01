### This script quantifies the expression of target genes in the CXCL cluster to assess their co-expression ###
# It uses the "quantNuclei.R" custom function to quantify the signal in a pre-selected subnuclear masks coresponding to the co-expression locus across target genes
# Thresholding of each channel is applied to define transcription as active or inactive based on intensity distribution inspection and visual inspection

# Input data (stitched images):
#1 cellpose nuclear masks segmented from DAPI channel (_corr_cp_masks.png)
#2 pre-segmented co-expression masks across target genes (_locus_masks.tif)
#3 multimeasure result table from FIJI across all 5 channels used in the padFISH experiment (Multimeasure_result_table.csv)

# load required packages
library(EBImage)
library(tidyr)
library(dplyr)
source("/filepath/quantNuclei_v01.R")

#assign directory where files are located
nucmask_dir <- "/select/filepath/"        #select segmented nuclear mask from cellpose - (filename: _corr_cp_masks.png)
locusmask_dir <- "/select/filepath/"      #select pre-segmented co-expression locus mask from ilastik+FIJI - (filename: _locus_masks.tif)

#retrieve all files in directory 
allfilesmask <- list.files(nucmask_dir, recursive = TRUE)
allfileslocus <- list.files(locusmask_dir, recursive = TRUE)
#select nucmask file 
nucmask_sel <- allfilesmask[grep("_corr_cp_masks.png$", allfilesmask)]
#select pre-segmented locus mask
locusmask_sel <- allfileslocus[grep("_locus_masks.tif$", allfileslocus)]

#load images
nucmasks <- readImage(paste0(nucmask_dir, nucmask_sel))
locusmasks_ori <- readImage(paste0(locusmask_dir, locusmask_sel))
#reformat nuclear masks and display in Viewer
nucmasks <-Image(as.numeric(as.factor(nucmasks))-1,dim = dim(nucmasks))
display(nucmasks)


#remove nuclei on edge of image
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
nucmasks <- rmObjects((nucmasks), ids)
storage.mode(nucmasks) <- 'double'
#display image to check result (If needed, reassign nucmasks and repeat with different border_size)
display(nucmasks)

#print initial number of cells used in the analysis (save for statistics)
Ncell_initial <- length(unique(as.vector(nucmasks)))
print(Ncell_initial)

#Assign co-expression locus masks and assign labels
locusmasks_bin <- locusmasks_ori
locus_labels <- bwlabel(locusmasks_bin)

#assign locusID (coexpression masks) within their corresponding nuclear maskID
locusID_nuc <- quantNuclei(nucmasks, locus_labels, unique) %>% data.frame() %>% setNames('locusID')
locusID_nuc$maskID <- as.integer(rownames(locusID_nuc))
locusID_nuc <- locusID_nuc %>% unnest(cols='locusID') %>% filter(locusID != 0)
#check locusID_nuc 
View(locusID_nuc)

#Intensity quantification 
#load intensity table obtained from FIJI Multimeasure (previously applied on the 5-channel image)
intensity_table <- read.csv(paste0(locusmask_dir,"Multimeasure_result_table.csv"))
#first create on the initial file a new column with the row numbers
inttable <- intensity_table %>% group_by(Ch) %>% 
  mutate(row_number=row_number())
#then subset the table to only have the row number, mean and channel number
inttable2 <- subset(inttable, select = c(row_number,Mean,Ch))
#then create a wide table, where each channel gets its own column
all_ch_inttable <- inttable2 %>%
  pivot_wider(names_from = Ch, values_from = Mean, names_prefix = "mean_ch")

#merge the results tables and apply names to columns (order of channels is defined from the multi-channel image used)
final_table <- merge(locusID_nuc,all_ch_inttable,by.x="locusID",by.y = "row_number")
colnames(final_table) <- c("locus", "cell", "CXCL3", "DAPI", "CXCL1", "CXCL8", "CXCL2")
#exclude DAPI from the following analysis
final_table[, "DAPI"] <- NULL
#add size info per co-expression locus 
colnames(intensity_table) <- c("locus", "area", colnames(intensity_table)[3:ncol(intensity_table)])
final_table <- merge(final_table, intensity_table[, 1:2], by="locus")
#add overall intensity per locus
final_table$intensity <- final_table[, grep("CXCL", colnames(final_table))] %>% rowSums()  
#save intensity table in work directory 
write.csv(final_table, "/output/filepath/_final_table_intensity.csv", row.names = FALSE)

#Plot intensity values and area from alleles to define filtering parameters
library(ggplot2)
ggplot(final_table, aes(x = area, y = intensity)) + geom_point()
ggplot(final_table, aes(x = area)) + geom_histogram()
ggplot(final_table, aes(x = intensity)) + geom_histogram()
#inspect also with log10
ggplot(final_table, aes(x = area %>% log10(.))) + geom_histogram()
ggplot(final_table, aes(x = intensity %>% log10(.))) + geom_histogram()

#Save plots as interim data in output folder
#save intensity/area scatter plot
plot_intensity_area <- ggplot(final_table, aes(x = area, y = intensity)) + geom_point()
ggsave(filename = "/output/filepath/plot_intensity_area.pdf", plot = plot_intensity_area)
#save area histogram
histogram_area <- ggplot(final_table, aes(x = area)) + geom_histogram()
ggsave(filename = "/output/filepath/histogram_area.pdf", plot = histogram_area)
#save intensity histogram
histogram_intensity <- ggplot(final_table, aes(x = intensity)) + geom_histogram()
ggsave(filename = "/output/filepath/histogram_intensity.pdf", plot = histogram_intensity)

#Filtering step:
#check how many cells have more than 2 loci in the unfiltered table (final_table)
sum(table(final_table[, "cell"]) > 2)      
#define exclusion parameters based on locus area and sum of intensity. Needed to help filter out the cells with more than 2 loci)
#define the values from inspection of the area and intensity plots to exclude outliers
final_table$allele[final_table$area < 10 | final_table$area > 200 | final_table$intensity > 15000 ]  #insert values
#filter the final_table based on both area and intensity thresholds
final_table_filt1 <- final_table[final_table$area > 10 & final_table$area < 200 , ]
final_table_filt2 <- final_table_filt1[final_table_filt1$intensity < 15000 , ]

### check how many cells have more than 2 loci in the filtered table 
sum(table(final_table_filt2[, "cell"]) > 2)  
### Check which cells have more than 2 loci, then exclude them from the table
cellsMore3 <- which(table(final_table_filt2[, "cell"]) > 2) %>% names(.) %>% as.numeric(.)
cellsMore3
final_table_filt3 <- final_table_filt2[!(final_table_filt2$cell %in% cellsMore3),]
#write Final_Table after cell filtering steps
Final_Table <- final_table_filt3[final_table_filt3[, "cell"] != 0, ]
write.csv(Final_Table, "/output/filepath/_Final_Table.csv", row.names = FALSE)
#check number of cells considered in the analysis (save for statistics)
sum(table(Final_Table[, "cell"]))

#Define intensity threshold for each genes to define if transcription is active or inactive
#Adjust threshold by both inspecting the intensity distribution and from visual inspection of images
CXCL1_cutoff <- 500    #insert value
hist(Final_Table[, "CXCL1"], breaks = 100)
abline(v = CXCL1_cutoff, col = "red")  #check position of threshold on the distribution
text(x = CXCL1_cutoff*2, y = 10, 
     labels = paste(sum(final_table[, "CXCL1"] >= CXCL1_cutoff) / nrow (final_table) * 100, "%"), col = "red")

CXCL2_cutoff <- 800   #insert value
hist(Final_Table[, "CXCL2"], breaks = 100)
abline(v = CXCL2_cutoff, col = "red")   #check position of threshold on the distribution
text(x = CXCL2_cutoff*2, y = 10, 
     labels = paste(sum(final_table[, "CXCL2"] >= CXCL2_cutoff) / nrow (final_table) * 100, "%"), col = "red")

CXCL3_cutoff <- 300    #insert value
hist(Final_Table[, "CXCL3"], breaks = 100)
abline(v = CXCL3_cutoff, col = "red")    #check position of threshold on the distribution
text(x = CXCL3_cutoff*1.5, y = 30, 
     labels = paste(sum(final_table[, "CXCL3"] >= CXCL3_cutoff) / nrow (final_table) * 100, "%"), col = "red")

CXCL8_cutoff <- 1000    #insert value
hist(Final_Table[, "CXCL8"], breaks = 100)
abline(v = CXCL8_cutoff, col = "red")   #check position of threshold on the distribution
text(x = CXCL8_cutoff*5, y = 30, 
     labels = paste(sum(final_table[, "CXCL8"] >= CXCL8_cutoff) / nrow (final_table) * 100, "%"), col = "red")

#Binarize the gene activation for 1 = active and 0 = inactive to obtain the co-expression pattern of the 4 genes of interest
### assign Final_Table to binary_table to plot the binary gene combinations
binary_table <- Final_Table
binary_table[, "CXCL2"] <- ifelse(binary_table[, "CXCL2"] >= CXCL2_cutoff, 1, 0)
binary_table[, "CXCL1"] <- ifelse(binary_table[, "CXCL1"] >= CXCL1_cutoff, 1, 0)
binary_table[, "CXCL3"] <- ifelse(binary_table[, "CXCL3"] >= CXCL3_cutoff, 1, 0)
binary_table[, "CXCL8"] <- ifelse(binary_table[, "CXCL8"] >= CXCL8_cutoff, 1, 0)

binary_table[, "category"] <- ""
for (i in 1:nrow(binary_table)){
  binary_table[i, "category"] <- paste(binary_table[i, "CXCL1"], binary_table[i, "CXCL2"], collapse = ",")
  binary_table[i, "category"] <- paste(binary_table[i, "category"], binary_table[i, "CXCL3"], collapse = ",")
  binary_table[i, "category"] <- paste(binary_table[i, "category"], binary_table[i, "CXCL8"], collapse = ",")
}

#save binary table in output directory, then plot binary_table and save it in the same directory
write.csv(binary_table, "/output/filepath/_binary_table.csv", row.names = FALSE)
ggplot(binary_table) + geom_bar(aes(x = category)) + theme_minimal()
binary_plot <- ggplot(binary_table) + geom_bar(aes(x = category)) + theme_minimal()
ggsave(filename = "/output/filepath/_binary_table_plot.pdf", plot = binary_plot)

#average number of observations (co-expressions) on the number of loci, get the percentage, write stats table to combine data sets
stats <- binary_table %>% count(category) %>% mutate(perc = n/nrow(binary_table)*100)
stats
binary_plot_perc <- ggplot(stats) + geom_col(aes(x = category, y = perc)) + theme_minimal()
binary_plot_perc
ggsave(filename = "/output/filepath/Repx_binary_plot_perc.pdf", plot = binary_plot_perc)
write.csv(stats, "/output/filepath/Repx_averaged_stats.csv", row.names = FALSE)


##### After  stats table are generated for all three conditions (0-30-240 min) - plot data for visualization #####

library(tidyr)
library(dplyr)
library(ggplot2)
library(matrixStats)


#load stats tables from replicate 1 - 2 - 3 
stats1 <- read.csv(paste0("/output/filepath/Rep1_averaged_stats.csv"))
stats2 <- read.csv(paste0("/output/filepath/Rep2_averaged_stats.csv"))
stats3 <- read.csv(paste0("/output/filepath/Rep3_averaged_stats.csv"))

# merge stats1, stats2, stats3 which come from each data set
stats <- merge(stats1, stats2, by = "category", suffixes = c("_run1", "_run2"))
colnames(stats3)[2:3] <- paste0(colnames(stats3)[2:3], "_run3")
stats <- merge(stats, stats3, by = "category")
stats


#OPTIONAL: ## create average between the three runs from the percentage value. Obtain then standard deviation (for bar plot)
## NOTE: rowMeans and rowSds require matrix or array ###
#stats$avg_perc <- rowMeans(stats[, c("perc_run1", "perc_run2", "perc_run3")])
#stats$sd_perc <- rowSds(stats[, c("perc_run1", "perc_run2", "perc_run3")])
#stats$avg_perc <- as.matrix(stats[, c("perc_run1", "perc_run2", "perc_run3")]) %>% rowMeans(.)
#stats$sd_perc <- as.matrix(stats[, c("perc_run1", "perc_run2", "perc_run3")]) %>% rowSds(.)
#write.csv(stats, "/output/filepath/_stats_avg_sd_perc.csv", row.names = FALSE)



