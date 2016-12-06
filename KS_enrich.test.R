KS_enrich.test <- function(myTerms, hitList, threshold=0.01){

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

        # geneIDs <- myTerms[["GO:2001275"]]
        ks.res <- lapply(myTerms2, function(geneIDs){
                #print(geneIDs)
                hasIndex <- has.key(geneIDs,hHitlist)
                if (sum(hasIndex) > 0){
                        quantiles <- values(hHitlist[geneIDs[hasIndex]])/geneUniverse
                        ks.test(x=quantiles, y=punif, alternative = "greater")
                } else {
                        return(NULL)
                }
        })

        ks.pvals <- lapply(ks.res, function(x){
                x$p.value
        })

        ks.pvals

}


