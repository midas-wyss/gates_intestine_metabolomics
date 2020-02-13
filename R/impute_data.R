#' impute_data
#' Data imputation for metabolomics. Given a matrix of intensities for a metabolomics experiment, this function will replace empty or `NA` values.
#' @param intensity_matrix A matrix of intensities for a metabolomics experiment
#' @return An imputed matrix with no `NA` values
impute_data <- function(intensity_matrix) {
  
  # error rate (defaults to 1%)
  err <- 0.01 
  
  # find minimal intensity (assumes this is the lowest value that can be measured above limit of detection)
  min_int <- min(intensity_matrix, na.rm = TRUE)
  
  # define error interval, given error rate
  err_int <- c(min_int - (min_int * err), min_int + (min_int * err))
  
  
  # get minimum of columns (samples)
  # min_samps <- apply(intensity_matrix, 2, min, na.rm = TRUE)
  
  # get average and sd of mins
  # avg_min <- mean(min_samps)
  # sd_min <- sd(min_samps)
  
  # find all NA values
  all_nas <- which(is.na(intensity_matrix))
  
  # generate values for NAs
  imp_values <- runif(n = length(all_nas), min = err_int[1], max = err_int[2])
  
  # get imputed values
  # imp_values <- rnorm(n = length(all_nas), mean = avg_min, sd = sd_min)
  
  # put them in matrix
  intensity_matrix[all_nas] <- imp_values
  
  # return imputed matrix
  return(intensity_matrix)
}