transformSettoList <- function(input_set){
        #Create a list object with sets as names and a vector of elements in
        #each list element
        geneList <- lapply(input_set[,2], function(x){
                unlist(strsplit(x, ","))
        })
        names(geneList) <- input_set[,1]
        return(geneList)
}
