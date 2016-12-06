buildDataMatrix <- function(inputMatrix, minSetSize = 3){

        setList <- inputMatrix[,1]
        geneList <- sort(unique(unlist(strsplit(inputMatrix[,2], ","))))

        set_matrix <- matrix(data = FALSE, ncol = length(setList), nrow = length(geneList), dimnames = list(geneList, setList))

        hGeneList <- hash(geneList, 1:length(geneList))
        hSetList <- hash(setList, 1:length(setList))

        apply(inputMatrix, 1, function(x) {
                iter_genes <- unlist(strsplit(as.character(x[2]), ","))
                for (i in 1:length(iter_genes)){
                        set_matrix[hGeneList[[iter_genes[i]]], hSetList[[x[1]]]] <<- TRUE
                }
        })
        (set_matrix <- set_matrix[,!(colSums(set_matrix) < minSetSize)])
}
