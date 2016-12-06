KS_enrich.elim.test <- function(myTerms, hitList, threshold=0.01){

        geneUniverse <- nrow(hitList)

        hHitlist <- hash(rownames(hitList), 1:geneUniverse)

        small_p <- 1e-100

        best_sets <- NULL

        myTerms2 <- lapply(myTerms, function(x){
                if (length(x) > 4){
                        return(x)
                } else {
                        return(NULL)
                }
        })
        myTerms2 <- Filter(length, myTerms2)

        i <- 1

        while (length(myTerms2)>0 & small_p <= threshold){
                cat("Elimination round:",i, "remaining test sets:", length(myTerms2),"\n")
                i <- i+1
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

                ks.pvals <- unlist(ks.pvals)

                best_sets <- c(best_sets, ks.pvals[ks.pvals == min(ks.pvals)])
                small_p <- min(ks.pvals)

                myTerms2 <- eliminateFromSetMatrix(myTerms2, names(ks.pvals[ks.pvals == min(ks.pvals)]))

                myTerms2 <- lapply(myTerms2, function(x){
                        if (length(x) > 4){
                                return(x)
                        } else {
                                return(NULL)
                        }
                })

                myTerms2 <- Filter(length, myTerms2)
        }


        best_sets

}


