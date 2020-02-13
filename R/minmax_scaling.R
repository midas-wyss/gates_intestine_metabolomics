#' metabolomics scale
#' Scale metabolomics data based on a min-max normalization. The min-max normalization is given by:
#' 
#' \deqn{M = \frac{x - min(x)}{max(x) - min(x)}
#' 
#' @param data_matrix A matrix to be scaled
#' @return Min-max scaled matrix
minmax_scaling <- function(data_matrix) {
  a1 <- apply(data_matrix, 2, min)
  a2 <- apply(data_matrix, 2, max)
  
  mm <- matrix(0, ncol = ncol(data_matrix), nrow = nrow(data_matrix))
  for (i in seq(1, ncol(data_matrix))) {
    mm[, i] <- (data_matrix[, i] - a1[i]) / (a2[i] - a1[i])
  }
  
  return(mm)
}