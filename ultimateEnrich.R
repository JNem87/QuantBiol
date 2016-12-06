ultimateEnrich <- function(hitList, minSetSize = 4, threshold=.05, type="GO", test = "fisher", maxSetSize=10000, GeneUniverse= NULL, minHitNum = 2, ...){

        cat("Start UltimateEnrich\n")
        cat("Settings:\n")
        cat("minimum set size =",minSetSize, "\n")
        cat("threshold =",threshold, "\n")
        cat("Enrich for =", type, "\n")
        cat("Tests =", test, "\n")

        myTerms <- NULL

        if ("GO" %in% type){
                if (!file.exists("data/goLookup.RData")){
                        cat("Retrieve GO sets from database:\n")
                        myTerms <- getGOMatrix(goType = "BP")
                        myGenes <- transformSettoList(myTerms[[2]])
                        myTerms <- myTerms[[1]]
                        cat("Retrieve GO sets done.\n\n")
                } else {
                        load("data/goLookup.RData")
                        myGenes <- transformSettoList(goLookup[[2]])
                        myTerms <- goLookup[[1]]
                }
        }

        if ("PATHWAY" %in% type){
                if (!file.exists("data/pathwayLookup.RData")){
                        cat("Retrieve PATHWAY sets from database:\n")
                        myTerms <- getPathwayMatrix()
                        myGenes <- transformSettoList(myTerms[[2]])
                        myTerms <- myTerms[[1]]
                        cat("Retrieve PATHWAY sets done.\n\n")
                } else {
                        load("data/pathwayLookup.RData")
                        myGenes <- transformSettoList(pathwayLookup[[2]])
                        myTerms <- pathwayLookup[[1]]
                }
        }

        if ("COMPLEX" %in% type){
                if (!file.exists("data/complexLookup.RData")){
                        cat("Retrieve COMPLEX sets from database:\n")
                        myTerms <- getComplexMatrix()
                        myGenes <- transformSettoList(myTerms[[2]])
                        myTerms <- myTerms[[1]]
                        cat("Retrieve COMPLEX sets done.\n\n")
                } else {
                        load("data/complexLookup.RData")
                        myGenes <- transformSettoList(complexLookup[[2]])
                        myTerms <- complexLookup[[1]]
                }

        }



        if ("fisher" %in% test){
                cat("Generate 2x2 matrices for fisher test:\n")
                myMatrix <- buildTwobyTwoMatrix(myTerms, minSetSize = minSetSize, myGeneUniverse = GeneUniverse)
                cat("Generate 2x2 matrices done.\n\n")

                if (is.null(GeneUniverse)){
                        myGeneUniverse <- myMatrix[[2]]
                } else {
                        myGeneUniverse <- GeneUniverse
                }

                myMatrix <- myMatrix[[1]]

                fisher_hitlist <- hitList[rownames(hitList) %in% myGeneUniverse,,drop=F]
                fisher_hitlist <- rownames(fisher_hitlist[,,drop=F])
                #fisher_hitlist <- rownames(fisher_hitlist[fisher_hitlist < threshold,,drop=F],)

                cat("Perform fisher test:\n")
                fisher.res <- fisher_enrich.test(myMatrix, fisher_hitlist, myGeneUniverse, threshold = threshold, maxSetSize=maxSetSize, minHitNum = minHitNum)
                cat("Perform fisher test done.\n\n")
        }

        if ("fisher_elim" %in% test){
                cat("Generate 2x2 matrices for fisher test:\n")
                myMatrix <- buildTwobyTwoMatrix(myTerms, minSetSize = minSetSize, myGeneUniverse = GeneUniverse)
                cat("Generate 2x2 matrices done.\n\n")

                if (is.null(GeneUniverse)){
                        myGeneUniverse <- myMatrix[[2]]
                } else {
                        myGeneUniverse <- GeneUniverse
                }

                myMatrix <- myMatrix[[1]]

                fisher_hitlist <- hitList[rownames(hitList) %in% myGeneUniverse,,drop=F]
                fisher_hitlist <- rownames(fisher_hitlist[,,drop=F],)
                #fisher_hitlist <- rownames(fisher_hitlist[fisher_hitlist < threshold,,drop=F],)

                cat("Perform fisher test:\n")
                fisher.res <- fisher_enrich.elim.test(myMatrix, fisher_hitlist, myGeneUniverse, threshold = threshold, genematrix = myGenes, maxSetSize=maxSetSize, minHitNum=minHitNum)
                cat("Perform fisher test done.\n\n")
        }


        if ("ks" %in% test){
                myTerms <- transformSettoList(myTerms)

                myGeneUniverse <- unique(unlist(myTerms))

                myGeneUniverse <- myGeneUniverse[myGeneUniverse %in% rownames(hitList)]
                hitList <- hitList[rownames(hitList) %in% myGeneUniverse,,drop=F]
                #print(nrow(hitList)

                myTerms <- reduceToUniverse(myTerms, myGeneUniverse)


                cat("Perform ks test:\n")
                fisher.res <- KS_enrich.test(myTerms, hitList, threshold = threshold)
                cat("Perform ks test done.\n\n")
        }

        if ("ks_elim" %in% test){
                myTerms <- transformSettoList(myTerms)

                myGeneUniverse <- unique(unlist(myTerms))

                myGeneUniverse <- myGeneUniverse[myGeneUniverse %in% rownames(hitList)]

                myGenes <- myGenes[myGeneUniverse]

                myTerms <- reduceToUniverse(myTerms, myGeneUniverse)

                cat("Perform ks test:\n")
                fisher.res <- KS_enrich.elim.test(myTerms, myGenes, hitList, threshold = threshold)
                cat("Perform ks test done.\n\n")
        }


        if ("wks" %in% test){
                myTerms <- transformSettoList(myTerms)

                myGeneUniverse <- unique(unlist(myTerms))

                myGeneUniverse <- myGeneUniverse[myGeneUniverse %in% rownames(hitList)]
                hitList <- hitList[rownames(hitList) %in% myGeneUniverse,,drop=F]
                #print(nrow(hitList)

                myTerms <- reduceToUniverse(myTerms, myGeneUniverse)

                hitList2 <- as.vector(hitList[,1])

                names(hitList2) <- rownames(hitList)

                cat("Perform ks test:\n")
                fisher.res <- WKS.test(hitList2, myTerms, alternative = "greater")
                cat("Perform ks test done.\n\n")
        }


        tmpRes <- c(type, fisher.res)

        FDRpass <- !unlist(tmpRes[[4]]) > 1
        if(sum(FDRpass) == 0) {
                return(NULL)
        }
        tmpRes[[2]] <- tmpRes[[2]][as.vector(which(FDRpass))]
        tmpRes[[3]] <- tmpRes[[3]][as.vector(which(FDRpass))]
        tmpRes[[4]] <- tmpRes[[4]][as.vector(which(FDRpass))]

        tmpResTable <- unlist(tmpRes[[2]])
        tmpResTable2 <- sapply(tmpRes[[3]], function(x){
                c(x[1,1], sum(x[1,1], x[1,2]), x[1,1] / sum(x[1,1], x[1,2]))
        })

        tmpResTab <- cbind(names(tmpResTable), tmpResTable, t(tmpResTable2)[names(tmpResTable),])
        if(length(names(tmpResTable)) == 1) {
                tmpResTab <- as.data.frame(t(c(names(tmpResTable), tmpResTable, t(tmpResTable2)[names(tmpResTable),])))
        }

        colnames(tmpResTab) <- c("PATHID", "pval", "hits", "total", "cov")

        tmpResTab <- as.data.frame(tmpResTab)

        if ("GO" %in% type){
                f <- data.frame(set=names(tmpRes[[2]]), pvalue=unlist(tmpRes[[2]]))

                #f[order(f$pvalue, decreasing = T),]

                cls <- columns(Mus.musculus)
                cls <- cls[cls %in% c("GOID", "TERM")]


                #Define GOID to be the keytype to use
                kts <- keytypes(Mus.musculus)
                kt <- kts[kts=="GOID"]

                #get all keys
                ks <- as.character(f[,1])
                #ks <- ks[1:10]

                #query all ensembl_ids for all GOIDs
                res <- select(Mus.musculus, keys=ks, columns=cls, keytype=kt)

                res <- merge(res, f, by.x="GOID", by.y="set")

                colnames(res) <- c("PATHID", "PATHNAME")
                tmpResDFmerge <- merge(res, tmpResTab)

                tmpResDFmerge$pval <- -log10(as.numeric(as.character(tmpResDFmerge$pval)))
                tmpResDFmerge$hits <- as.numeric(as.character(tmpResDFmerge$hits))
                tmpResDFmerge$cov <- as.numeric(as.character(tmpResDFmerge$cov))


                tmpResDFmerge <- tmpResDFmerge[order(tmpResDFmerge$pval),]
        }

        if ("PATHWAY" %in% type){
                cls <- columns(reactome.db)
                cls <- cls[cls %in% c("PATHID", "PATHNAME")]

                #Define GOID to be the keytype to use
                kts <- keytypes(reactome.db)
                kt <- kts[kts=="PATHID"]

                #get all keys
                #ks <- AnnotationDbi::keys(reactome.db, keytype=kt)
                # ks <- ks[1:10]
                ks=names(sort(unlist(tmpRes[[2]])))

                #query all ensembl_ids for all GOIDs
                res <- select(reactome.db, keys=ks, columns=cls, keytype=kt)

                #tmpResTab
                #tmpResDF <- data.frame(PATHID=names(tmpRes), pvalue=unlist(tmpRes[[1]]))

                colnames(res) <- c("PATHID", "PATHNAME")
                tmpResDFmerge <- merge(res, tmpResTab)

                tmpResDFmerge$pval <- -log10(as.numeric(as.character(tmpResDFmerge$pval)))
                tmpResDFmerge$hits <- as.numeric(as.character(tmpResDFmerge$hits))
                tmpResDFmerge$cov <- as.numeric(as.character(tmpResDFmerge$cov))


                tmpResDFmerge <- tmpResDFmerge[order(tmpResDFmerge$pval),]
                tmpResDFmerge$PATHNAME <- substr(tmpResDFmerge$PATHNAME, 15, nchar(tmpResDFmerge$PATHNAME))
        }

        if ("COMPLEX" %in% type){
                allComplexes <- read.csv2(file = "data2/allComplexes.csv")

                compleatComplexes <- read.delim("data2/compleatComplexes.txt", header=F, stringsAsFactors = F)

                annoComplex <- allComplexes[allComplexes$Complex.id %in% names(unlist(tmpRes[[2]])),c(1,2)]
                if (nrow(annoComplex) != 0){
                        annoComplex$Complex.name <- substr(annoComplex$Complex.name,1,100)
                        annoComplex$Complex.name <- paste0(annoComplex$Complex.name, "(CORUM)")
                }

                annoComplex2 <- compleatComplexes[compleatComplexes$V1 %in% names(sort(unlist(tmpRes))),c("V1", "V6")]
                colnames(annoComplex2) <- colnames(annoComplex)
                annoComplex2$Complex.name <- substr(annoComplex2$Complex.name,1,100)
                annoComplex2$Complex.name <- paste0(annoComplex2$Complex.name, "(COMPLEAT)")

                annoComplex <- rbind(annoComplex,annoComplex2)

                colnames(annoComplex) <- c("PATHID", "PATHNAME")

                tmpResDFmerge <- merge(tmpResTab, annoComplex)

                tmpResDFmerge <- tmpResDFmerge[!tmpResDFmerge$PATHNAME == "Unknown(COMPLEAT)",]

                tmpResDFmerge$pval <- -log10(as.numeric(as.character(tmpResDFmerge$pval)))
                tmpResDFmerge$hits <- as.numeric(as.character(tmpResDFmerge$hits))
                tmpResDFmerge$cov <- as.numeric(as.character(tmpResDFmerge$cov))

                tmpResDFmerge <- tmpResDFmerge[order(tmpResDFmerge$pval),]

        }

        colnames(myTerms)[1] <- "Complex.id"

        hitGenesinSets <- lapply(as.character(tmpResDFmerge$PATHID), function(x){
                tmpSet <- unlist(strsplit(myTerms[myTerms$Complex.id == x,2],","))
                hitSet <- tmpSet[tmpSet %in% rownames(hitList)]
                list(hitSet, tmpSet)
        })

        names(hitGenesinSets) <- tmpResDFmerge$PATHID

        hitGenesinSetStrings <- lapply(hitGenesinSets, function(x){
                c(paste(x[[1]], collapse = ","), paste(x[[2]], collapse = ","))
        })

        hitGenesinSetTab <- do.call("rbind", hitGenesinSetStrings)

        colnames(hitGenesinSetTab) <- c("hitGenes", "SetGenes")

        hitGenesinSetMerge <- merge(tmpResDFmerge, hitGenesinSetTab, by.x="PATHID", by.y="row.names")



        # if (!is.null(ks.res)){
        #         return(ks.res)
        # }
        #
        # if (!is.null(fisher.res)){
        #         return(fisher.res)
        # }
        #
        # return(c(type, fisher.res))
        return(list(type, tmpResDFmerge, hitGenesinSetMerge))
}




# cpx.elim.pval.list <- data.frame()
# cpx_id_list <- colnames(cpx_matrix)

# fisher.test(as.factor(hitList), as.logical(cpx_matrix[,1]), alternative = "greater")
# first_fisher <- lapply(as.data.frame(cpx_matrix),function(x){
#         fisher.test(as.factor(hitList), as.logical(x), alternative = "greater")
# })
#system.time(fisher.test(as.factor(hitList), as.logical(cpx_matrix[,1]), alternative = "greater"))

#get standard fisher-values first
# first_fisher <- apply(as.data.frame(cpx_matrix),2, function(x){
#         fisher.test(hitList, x, alternative = "greater")
# })
#
# first_fisher <- lapply(first_fisher, function(x){
#         return(x$p.value)
# })
# first_fisher <- do.call(c, first_fisher)
# first_fisher <- as.data.frame(cbind(p_value=first_fisher, names(first_fisher)))
# colnames(first_fisher) <- c("pvalue_fisher", "Complex.id")
#
# first_fisher <- transform(first_fisher, pvalue_fisher =as.numeric(as.character(first_fisher$pvalue_fisher)))

# while (ncol(cpx_matrix)>0){
#         cpx.pvals <- apply(as.data.frame(cpx_matrix),2, function(x){
#                 fisher.test(as.factor(hitList), as.logical(x), alternative = "greater")
#         })
#
#         tmp.cpx.pval.res <- lapply(cpx.pvals, function(x){
#                 return(x$p.value)
#         })
#
#         tmp.cpx.pval.res <- do.call(c, tmp.cpx.pval.res)
#
#         tmp.cpx.pval.res <- as.data.frame(cbind(p_value=tmp.cpx.pval.res, names(tmp.cpx.pval.res)))
#
#         colnames(tmp.cpx.pval.res)[2] <- "Complex.id"
#
#         tmp.cpx.pval.res <- transform(tmp.cpx.pval.res, p_value=as.numeric(tmp.cpx.pval.res$p_value))
#
#         tmp_cpx_id <- tmp.cpx.pval.res$Complex.id[tmp.cpx.pval.res$p_value == min(tmp.cpx.pval.res$p_value)]
#
#         #   cpx.elim.pval.list <<- rbind(cpx.elim.pval.list,  c(complex_id=paste(tmp_cpx_id,collapse = ","),  p_value=unique(tmp.cpx.pval.res$p_value[tmp.cpx.pval.res$Complex.id %in% tmp_cpx_id])))
#         cpx.elim.pval.list <<- rbind(cpx.elim.pval.list,  cbind(tmp_cpx_id, tmp.cpx.pval.res$p_value[tmp.cpx.pval.res$Complex.id %in% tmp_cpx_id]))
#
#         if (length(tmp_cpx_id) > 1){
#                 cpx_matrix[rowSums(cpx_matrix[,which(colnames(cpx_matrix) %in% tmp_cpx_id)]) > 0,] <- 0
#         }
#         else {
#                 cpx_matrix[cpx_matrix[,which(colnames(cpx_matrix) == tmp_cpx_id)] > 0,] <- 0
#         }
#
#
#         cpx_id_list <- cpx_id_list[-which(colSums(cpx_matrix) < 4)]
#         cpx_matrix <- cpx_matrix[,-which(colSums(cpx_matrix) < 4)]
#
#         if (class(cpx_matrix) != "matrix") {
#                 cpx_matrix <- as.matrix(cpx_matrix)
#                 colnames(cpx_matrix) <- cpx_id_list
#         }
# }
# fisher2 <-  merge(first_fisher, cpx.elim.pval.list, by="Complex.id", all.x = T)
# return(fisher2[order(fisher2$p_value),])
















# head(first_fisher[order(first_fisher$pvalue_fisher),],10)
#
# ##for loop for the permutation based distribution
# for (i in 1:10) {
#         print(i)
#         cpx_matrix <- cpx_matrix_tmp
#         sig.genes.for.fisher <- sample(sig.genes.for.fisher)
#
#         cpx.elim.pval.list <- data.frame()
#         cpx_id_list <- colnames(cpx_matrix_tmp)
#
#         while (ncol(cpx_matrix)>0){
#                 cpx.pvals <- apply(as.data.frame(cpx_matrix),2, function(x){
#                         fisher.test(sig.genes.for.fisher, as.logical(x), alternative = "greater")
#                 })
#
#                 tmp.cpx.pval.res <- lapply(cpx.pvals, function(x){
#                         return(x$p.value)
#                 })
#
#                 tmp.cpx.pval.res <- do.call(c, tmp.cpx.pval.res)
#
#                 tmp.cpx.pval.res <- as.data.frame(cbind(p_value=tmp.cpx.pval.res, names(tmp.cpx.pval.res)))
#
#                 colnames(tmp.cpx.pval.res)[2] <- "Complex.id"
#
#                 tmp.cpx.pval.res <- transform(tmp.cpx.pval.res, p_value=as.numeric(tmp.cpx.pval.res$p_value))
#
#                 tmp_cpx_id <- tmp.cpx.pval.res$Complex.id[tmp.cpx.pval.res$p_value == min(tmp.cpx.pval.res$p_value)]
#
#                 #   cpx.elim.pval.list <<- rbind(cpx.elim.pval.list,  c(complex_id=paste(tmp_cpx_id,collapse = ","),  p_value=unique(tmp.cpx.pval.res$p_value[tmp.cpx.pval.res$Complex.id %in% tmp_cpx_id])))
#                 cpx.elim.pval.list <<- rbind(cpx.elim.pval.list,  cbind(tmp_cpx_id, tmp.cpx.pval.res$p_value[tmp.cpx.pval.res$Complex.id %in% tmp_cpx_id]))
#
#                 if (length(tmp_cpx_id) > 1){
#                         cpx_matrix[rowSums(cpx_matrix[,which(colnames(cpx_matrix) %in% tmp_cpx_id)]) > 0,] <- 0
#                 }
#                 else {
#                         cpx_matrix[cpx_matrix[,which(colnames(cpx_matrix) == tmp_cpx_id)] > 0,] <- 0
#                 }
#
#
#                 cpx_id_list <- cpx_id_list[-which(colSums(cpx_matrix) < 4)]
#                 cpx_matrix <- cpx_matrix[,-which(colSums(cpx_matrix) < 4)]
#
#                 if (class(cpx_matrix) != "matrix") {
#                         cpx_matrix <- as.matrix(cpx_matrix)
#                         colnames(cpx_matrix) <- cpx_id_list
#                 }
#         }
#
#         colnames(cpx.elim.pval.list) <- c("Complex.id", "p_value")
#
#         random.cpx.list <- c(random.cpx.list, as.numeric(cpx.elim.pval.list$p_value))
# }


