# padFISH
Collection of scripts used for Preprocessing and Quantification for padFISH data set in the study "Two different chromatin modules regulate proinflammatory gene expression" from Seufert et al.

> padFISH - multiplex smFISH of nascent RNAs 
padFISH is a multiplexed single-molecule fluorescence in situ hybridization (smFISH) method that combines elements from hybridization-based in situ sequencing (HybISS) [1] and single-cell resolution in situ hybridization on tissues (SCRINSHOT) [2]. The padFISH technique enables high-resolution tracing of nascent RNAs using intronic padlock probes (PLPs) that target cDNA, followed by rolling circle amplification (RCA) for signal enhancement. 

> Preprocessing
1.	For the preprocessing of raw images (Imaris format) and metadata, image stacks were first transformed into maximum projected TIF files by running the RunIJ_imstoTif_GPU_v3.R  script and imstoTif_headless_v2.ijm FIJI macro.
2.	The 02_Preprocessing_v1.6.ijm FIJI macro performed flatfield correction, chromatic aberration correction and stitching (via Grid/Collection Stitching FIJI plugin). Stitched images were used as input in all subsequent analysis.
3.	Nuclei segmentation on DAPI images was performed with Cellpose 2 [3] using the pretrained cyto model with diameter 150 and 200 for 60x and 100x objectives, respectively. To run Cellpose in batch the RunCellpose_batch_GPU_v3.R was used.
4.	Cell nuclei at the borders of the image or that displayed overexposure in individual channels were filtered out in R prior to further analysis (see Quantification section). 

> Quantification
1.	Individual channel images and Cellpose nuclear masks were used as input in R to quantify image features in regions corresponding to nuclear masks using the custom function quantNuclei_v01.R [4]. 
2.	For padFISH transcriptional bursting analysis, we computed the sum of fluorescence intensities in each nucleus using the padFISH_Bursting_kinetics_analysis.R script. 
3.	For the padFISH co-expression analysis at the CXCL cluster, the padFISH_CXCL_coexpression_analysis.R script was used.
4.	To quantify Immunofluorescence (IF) data, we used IF_Integrated_Intensity_v2.R.


References 
1.	Gyllborg, D. et al. Hybridization-based in situ sequencing (HybISS) for spatially resolved transcriptomics in human and mouse brain tissue. Nucleic Acids Res 48, e112 (2020)
2.	Sountoulidis, A. et al. SCRINSHOT enables spatial mapping of cell states in tissue sections with single-cell resolution. PLoS Biol 18, e3000675 (2020)
3.	Pachitariu, M. & Stringer, C. Cellpose 2.0: how to train your own model. Nat Methods 19, 1634-1641 (2022)
4.	Frank, L. Perturbing and imaging nuclear compartments to reveal mechanisms of transcription regulation and telomere maintenance. PhD thesis, Heidelberg University (2023)
