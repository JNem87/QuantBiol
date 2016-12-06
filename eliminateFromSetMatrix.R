eliminateFromSetMatrix <- function(cpx_matrix, genematrix, set_ids){
        hSet <- hash(names(cpx_matrix), 1:length(cpx_matrix))
        hGene <- hash(names(genematrix), 1:length(genematrix))
        newmat <- NULL
        tmp_list <- lapply(set_ids, function(setA){
                removeList <<- cpx_matrix[hSet[[setA]]][[1]]
                newGenes <- lapply(removeList, function(geneA){
                        removeSets <<- unlist(genematrix[hGene[[geneA]]])
                        newmat <<- lapply(removeSets, function(remSet){
                                cpx_matrix[[remSet]][!cpx_matrix[[remSet]] %in% geneA]
                        })
                        newmat
                })
                # newmat <<- lapply(cpx_matrix, function(oldSet){
                #         newSet <- oldSet[!oldSet %in% removeList]
                # })
        })

        return(newmat)
}
# remSet <- "25"
# setA <- "9"
#
# set_ids <- "9"
# geneA <- "ENSMUSG00000019256"
# cpx_matrix <- myTerms
# genematrix <- myGenes
# cpx_matrix[1]
# cpx_matrix[[1]][!cpx_matrix[[1]] %in% "ENSMUSG00000028982"]
