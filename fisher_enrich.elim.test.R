#perform fisher enrichment test on the provided gene sets
fisher_enrich.elim.test <- function(cpx_matrix, hitList, myGeneUniverse, minSetSize = 4, threshold=.05, genematrix, maxSetSize=10000, minHitNum=2){

        numhits <- buildTwobyTwoMatrix2(cpx_matrix, hitList, myGeneUniverse, minSetSize, maxSetSize, minHitNum)
        numhits_orig <- numhits

        # numhits3 <- lapply(numhits, function(nowzero){
        #         if(nowzero[1,1] == 0){
        #                 return(NULL)
        #         }
        #         else {
        #                 return(nowzero)
        #         }
        # })

        #numhits <- Filter(length, numhits3)

        best_ps <- NULL
        i = 1
        small_p <- 1e-100

        while (length(numhits)>0 & small_p <= threshold)
        {
                cat("Elimination round:",i, "remaining test sets:", length(numhits),"\n")
                i <- i+1
                pvalues <- lapply(numhits, function(mat){
                        fisher.test(mat, alternative = "greater")$p.value
                })

                small_p <- min(unlist(pvalues))

                best_ps <- c(best_ps, pvalues[pvalues==min(unlist(pvalues))])

                #genesToRemove <- hitList[hitList %in% unlist(cpx_matrix[names(pvalues[pvalues==min(unlist(pvalues))])])]
                genesToRemove <- unlist(cpx_matrix[names(pvalues[pvalues==min(unlist(pvalues))])])

                #hitList <- hitList[!hitList %in% genesToRemove]

                numhits2 <- eliminateFrom2by2Matrix(numhits, genematrix, genesToRemove, hitList)

                numhits3 <- lapply(numhits2, function(nowzero){
                        if(nowzero[1,1] < minHitNum | nowzero[1,2] < 1){
                                return(NULL)
                        }
                        else {
                                return(nowzero)
                        }
                })

                numhits <- Filter(length, numhits3)

                if (i%%50 == 0L) {
                        print(best_ps)
                }

        }

        #best_ps
        #best_ps[best_ps<=threshold]
        fdr_vals <- p.adjust(unlist(best_ps), method = "BH")
        list(best_ps[best_ps<=threshold],numhits_orig[names(numhits_orig) %in% names(best_ps[best_ps<=threshold])], fdr_vals[best_ps<=threshold])
}
