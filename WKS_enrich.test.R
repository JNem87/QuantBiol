WKS_enrich.test <- function(myTerms, hitList, threshold=1){

        geneUniverse <- nrow(hitList)

        hHitlist <- hash(rownames(hitList), 1:geneUniverse)

        myTerms2 <- lapply(myTerms, function(x){
                if (length(x) > 4){
                        return(x)
                } else {
                        return(NULL)
                }
        })

        myTerms2 <- Filter(length, myTerms2)

        hitList2 <- hitList[,1]
        names(hitList2) <- row.names(hitList)

        profvis({
                system.time(ks.res <- WKS.test(w=hitList2, db=myTerms2, alternative = "less"))
        })

        system.time(ks.res <- WKS.test(w=hitList2, db=myTerms2, alternative = "less"))

        system.time(ks.res <- WKS.test(w=hitList2, db=myTerms2, alternative = "less",ld = ld))

        ld <- ks.res[[2]]
        ks.res <- ks.res[[1]]

        return(as.list(ks.res))
}
