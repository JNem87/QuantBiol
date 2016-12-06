reduceToUniverse <- function(cpx_matrix, geneUniverse){
        hUniv <- hash(geneUniverse, 1:length(geneUniverse))
        newmat <<- lapply(cpx_matrix, function(oldSet){
                newSet <- oldSet[has.key(oldSet, hUniv)]
        })
}

