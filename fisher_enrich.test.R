#perform fisher enrichment test on the provided gene sets
fisher_enrich.test <- function(cpx_matrix, hitList, myGeneUniverse, minSetSize = 4, threshold=.05, maxSetSize=10000, minHitNum=2){

        numhits <- buildTwobyTwoMatrix2(cpx_matrix, hitList, myGeneUniverse, minSetSize, maxSetSize, minHitNum)

        pvalues <- lapply(numhits, function(mat){
                fisher.test(mat, alternative = "greater")$p.value
        })

        fdr_vals <- p.adjust(unlist(pvalues), method = "BY")

        list(pvalues[pvalues<=threshold],numhits[names(numhits) %in% names(pvalues[pvalues<=threshold])], fdr_vals[pvalues<=threshold])
}
