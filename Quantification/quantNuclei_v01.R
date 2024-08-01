#' This function applies any function/operation to quantify image features in regions corresponding nuclear masks
#' @param nuc_masks Image of nuclear masks (bg = 0, nuclei = 1 to n), 
#' @param im Image used for quantification
#' @param operation Function/operation run on the pixels regions of all nuclear masks
#' @returns Quantification results

quantNuclei <- function(nuc_masks, im, operation) {
  if(sum(as.numeric(dim(nuc_masks)==dim(im)))==2){
    print(paste0(Sys.time()," Quantifying intensities in masks..."))
    # Run operation on regions corresponding to each nuclear mask pixels
    quantification_results <- tapply(im,nuc_masks,operation)
    # Remove quantification results of pixel region with value 0 (= inverse of nuclear masks)
    quantification_results <- quantification_results[-which(as.numeric(names(quantification_results))==0)]
    # Return results
    print(paste0(Sys.time()," Done."))
    return(quantification_results)
  }else{
    warning("Aborted. Images do not have the same dimensions!")
  }
}
