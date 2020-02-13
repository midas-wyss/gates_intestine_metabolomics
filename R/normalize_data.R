#' normalize metabolomics data
#' 
#' Performs normalization of metabolomics data for each given metabolite across all samples. In this normalization, the media of each metabolite is set to 1 and the data is then corrected to this.
#'
#' @param data_matrix A matrix of metabolomics intensities
#' @param protein_data Logical flag on the existence of protein abundance data for a given sample, to perform normalization to protein concentration. Defaults to `FALSE`.
#' @param protein_quant A vector with the protein concentration values. Length of vector needs to be the same size as number of samples.
#' @return A normalized intensity matrix
normalize_data <- function(data_matrix, protein_data = FALSE, protein_quant = NULL) {
  
  if(!protein_data) {
    # count NAs
    x1 <- length(which(is.na(data_matrix))) 
    
    if (x1 != 0) {
      stop("Found NA values. Perform data imputation first.")
    } 
    
    # calculate median
    met_medians <- apply(data_matrix, 1, median, na.rm = TRUE)
    
    # perform normalization
    met_norm <- matrix(0, ncol = ncol(data_matrix), nrow = nrow(data_matrix))
    for (i in seq(1, nrow(data_matrix))) {
      met_norm[i, ] <- data_matrix[i, ] / met_medians[i]
    }    
    
  } else if (protein_data & is.null(protein_quant)) {
    
    stop("Need protein quantification.")
    
  } else {
    
    # count NAs
    x1 <- length(which(is.na(data_matrix))) 
    
    if (x1 != 0) {
      stop("Found NA values. Perform data imputation first.")
    } 
    
    if (length(protein_quant) != ncol(data_matrix)) {
      stop("Sizes of protein data and samples don't match. Please revise.")
    } else {
      # normalize to protein abundance
      for (i in seq(1, ncol(data_matrix))) {
        data_matrix[, i] <- data_matrix[, i] / protein_quant[i]
      }
      
      # calculate median
      met_medians <- apply(data_matrix, 1, median, na.rm = TRUE)
      
      # perform normalization
      met_norm <- matrix(0, ncol = ncol(data_matrix), nrow = nrow(data_matrix))
      for (i in seq(1, nrow(data_matrix))) {
        met_norm[i, ] <- data_matrix[i, ] / met_medians[i]
      }
    }
  }
  
  return(met_norm)
  
}