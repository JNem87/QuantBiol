# data = data to process - data frame, matrix (possibly data table?)
# rem.comp. = component to remove

# Based on Lindsay I Smith 2002
# http://www.cs.otago.ac.nz/cosc453/student_tutorials/principal_components.pdf

# It's usually preferable to use the basic svd function, since the function below is fairly slow.
# However, svd() can in rare cases fail to converge and quit with an error message.
# The below solution should avoid this.

remove_principal_component <- function(data, rem.comp.) {
  # preparation: extract means, center data
  data_colmeans <- apply(data, 2, mean)
  data_centered <- apply(data, 2, function(c) {c - mean(c)})
  
  # extraction of eigen values and vectors from data
  data_eigen <- eigen(cov(data_centered))$vectors
  
  # construct the new data set with the selected components removed
  newdata <- t(data_eigen[,-rem.comp.]) %*% t(data_centered)
  
  # reconstruct old data
  reconstructed <- t((solve(t(data_eigen[,-rem.comp.])) %*% newdata) + data_colmeans)
  
  dimnames(reconstructed) <- dimnames(data)
  
  if (class(data) == 'data.frame') {
    return(as.data.frame(reconstructed))
  } else {
    return(reconstructed)
  }
}