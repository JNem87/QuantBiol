buildTwobyTwoMatrix2 <- function(cpx_matrix, hitList, myGeneUniverse  = NULL, minSetSize = 4, maxSetSize = 10000, minHitNum=2){

        if (!is.null(myGeneUniverse)) {
                hUniverse <- hash(myGeneUniverse, 1:length(myGeneUniverse))

                hHitList <- hash(hitList, 1:length(hitList))

                hHitList <- hHitList[names(has.key(names(hHitList), hUniverse))[has.key(names(hHitList), hUniverse)]]
        } else{
                hHitList <- hash(hitList, 1:length(hitList))
        }



        length_hits <- length(hitList)
        length_universe <- length(myGeneUniverse)
        #sum(has.key(cpx_matrix[[1]],hHitList))
        numhits <- lapply(cpx_matrix, function(x){
                hits <- sum(has.key(unique(x),hHitList))
                if (hits >= minHitNum & length(x) >= minSetSize & length(x) <= maxSetSize) {
                        f11 <- hits
                        f12 <- length(x)-f11
                        f21 <- length_hits-f11
                        f22 <- length_universe - f11 - f12 -f21
                        matrix(c(f11,f12,f21,f22),nrow = 2, byrow = T)
                } else {
                        NULL
                }
        })
        numhits <- Filter(length, numhits)

        numhits
}
