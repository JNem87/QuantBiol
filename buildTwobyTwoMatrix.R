buildTwobyTwoMatrix <- function(inputMatrix, minSetSize = 3, myGeneUniverse  = NULL){

        #Create a list object with sets as names and a vector of elements in
        #each list element
        if (!is.null(myGeneUniverse)) {
                hHitList <- hash(myGeneUniverse, 1:length(myGeneUniverse))
        }

        genesInUniverse <- NULL

        geneList <- lapply(inputMatrix[,2], function(x){
                unlist(strsplit(x, ","))
        })
        names(geneList) <- inputMatrix[,1]


        #keep only genes in universe
        if (!is.null(myGeneUniverse)) {
                genesInUniverse <- lapply(geneList, function(x) {
                        x[has.key(x,hHitList)]
                })
                genesInUniverse <- Filter(length, genesInUniverse)
                geneList <- genesInUniverse
        }



        #count genes
        if (!is.null(myGeneUniverse)) {
                geneCounts <- lapply(geneList, function(x) {
                        sum(has.key(x,hHitList))
                })
        } else {
                geneCounts <- lapply(geneList, length)
        }


        geneList <- geneList[geneCounts > (minSetSize-1)]
        allGenes <- sort(unique(unlist(strsplit(inputMatrix[,2], ","))))
        return(list(setlist=geneList, numberOfGenesinSets=allGenes))
}
