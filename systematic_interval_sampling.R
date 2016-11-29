# The function below performs systematic interval sampling. Each column of a data frame x is sorted by value (0s are omitted) and divided into n bins.
# One value out of each bin is then randomly drawn.

systematic_sampling <- function(x,n){
  for(currcol in 1:ncol(x)) {
    
    # Prepare a temporary data frame with which to later sample x.
    to_process <- data.frame(values = x[,currcol], ID = rownames(x))
    to_process <- to_process[order(to_process$values),]
    to_process <- to_process[which(to_process$values != 0),]
    
    # Randomly generate cutoffs for bins. This assumes that due to n not being a divisor of the number of observations
    # some bins will have to be larger than others. Which bins have which size is randomly assigned below.
    if (nrow(to_process) < n) {x[,currcol] <- 0} else if (nrow(to_process) == n) {} else {
      binbreaks = rep(floor(nrow(to_process)/n),n)
      binbreaks[sample(1:length(binbreaks),nrow(to_process)%%n, replace = F)] = ceiling(nrow(to_process)/n)
      binbreaks = c(0,cumsum(binbreaks))
      
      # Randomly sample observation to retain within the bins.
      to_keep = sapply(1:(length(binbreaks)-1),function(x){
        if (binbreaks[x+1]-(binbreaks[x]+1) == 0) {binbreaks[x+1]} else {
          sample((binbreaks[x]+1):binbreaks[x+1],1)
        }
      })
      
      # Subset the temporary data frame.
      to_process <- to_process[to_keep,]
    }
    
    # Subset x.
    x[-which(rownames(x) %in% to_process$ID),currcol] <- 0
  }
  ## Drop columns with no values above 0.
  x <- x[,colSums(x) != 0]
}