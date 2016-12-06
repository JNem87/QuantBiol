getGOMatrix <- function(goType = "BP"){
        #Look for columns GOID and ENSEMBL in database
        cls <- columns(Mus.musculus)
        cls <- cls[cls %in% c("GOID","ENSEMBL")]

        #Define GOID to be the keytype to use
        kts <- keytypes(Mus.musculus)
        kt <- kts[kts=="ENSEMBL"]

        #get all keys
        ks <- AnnotationDbi::keys(Mus.musculus, keytype=kt)
        # ks <- ks[1:10]

        #query all ensembl_ids for all GOIDs
        res <- select(Mus.musculus, keys=ks, columns=cls, keytype=kt)

        if (goType == "ALL"){
                res2 <- res
        } else {
                res2 <- res[res$ONTOLOGY==goType & !is.na(res$ENSEMBL),]
        }

        #Aggregate by GOIDs
        goAggregate <- aggregate(res2$ENSEMBL, by = list(res2$GOID), FUN=paste,collapse=",")
        colnames(goAggregate) <- c("GOID", "EnsemblID")
        geneAggregate <- aggregate(res2$GOID, by = list(res2$ENSEMBL), FUN=paste,collapse=",")
        colnames(geneAggregate) <- c("EnsemblID", "GOID")

        goLookup <- list(goAggregate, geneAggregate)
        save(goLookup, file = "data/goLookup.RData")

        return(list(goAggregate, geneAggregate))
}









