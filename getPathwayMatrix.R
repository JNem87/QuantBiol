getPathwayMatrix <- function(){
        cls <- columns(reactome.db)
        cls <- cls[cls %in% c("ENTREZID","PATHID")]

        #Define GOID to be the keytype to use
        kts <- keytypes(reactome.db)
        kt <- kts[kts=="ENTREZID"]

        #get all keys
        ks <- AnnotationDbi::keys(reactome.db, keytype=kt)
        # ks <- ks[1:10]

        #query all ensembl_ids for all GOIDs
        res <- select(reactome.db, keys=ks, columns=cls, keytype=kt)

        #translate ids to ensembl
        cls <- columns(Mus.musculus)
        cls <- cls[cls %in% c("ENTREZID","ENSEMBL")]

        #Define GOID to be the keytype to use
        kts <- keytypes(Mus.musculus)
        kt <- kts[kts=="ENTREZID"]

        #get all keys
        ks <- AnnotationDbi::keys(Mus.musculus, keytype=kt)
        # ks <- ks[1:10]

        #query all ensembl_ids for all GOIDs
        res2 <- select(Mus.musculus, keys=ks, columns=cls, keytype=kt)

        res3 <- merge(res, res2, by.x="ENTREZID", by.y="ENTREZID")

        res3 <- res3[,c(3,2)]
        #Remove NAs
        res3 <- res3[!is.na(res3$ENSEMBL),]
        res3 <- unique(res3)

        pathAggregate <- aggregate(res3$ENSEMBL, by = list(res3$PATHID), FUN=paste,collapse=",")
        colnames(pathAggregate) <- c("PATHID", "ENTREZID")
        pathAggregate2 <- aggregate(res3$PATHID, by = list(res3$ENSEMBL), FUN=paste,collapse=",")
        colnames(pathAggregate2) <- c("ENTREZID", "PATHID")

        pathwayLookup <- list(pathAggregate, pathAggregate2)
        save(pathwayLookup, file="data/pathwayLookup.RData")

        return(list(pathAggregate, pathAggregate2))
}
