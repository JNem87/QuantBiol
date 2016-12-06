eliminateFrom2by2Matrix <- function(cpx_matrix, genematrix, gene_ids, hitList){
        hGene <- hash(names(genematrix), 1:length(genematrix))
        hHits <- hash(hitList, 1:length(hitList))
        tmp_list <- lapply(gene_ids, function(geneA){
                newmat <- lapply(genematrix[hGene[[geneA]]][[1]], function(setA){
                        mat2 <- cpx_matrix[[setA]]
                        if (has.key(hHits, key = geneA)){
                                mat2[1,1] <- mat2[1,1] - 1
                                mat2[2,1] <- mat2[2,1] + 1
                        } else {
                                mat2[1,2] <- mat2[1,2] - 1
                                mat2[2,2] <- mat2[2,2] + 1
                        }
                        cpx_matrix[[setA]] <<- mat2
                })
        })
        return(cpx_matrix)
}
