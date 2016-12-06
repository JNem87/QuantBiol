getComplexMatrix <- function(){
        #load complexes
        allComplexes <- read.csv2(file = "data2/allComplexes.csv")

        head(allComplexes)


        #query homologs for proteins
        human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

        #query ensembl entrez lookup
        #mouse.lookup <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), mart=mouse)
        mouse.lookup <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "mgi_symbol", "uniprot_swissprot"), mart=mouse)

        #generate full gene list
        #full.genes.for.fisher <- count_TTAA_datamerge[,c(1,4)]
        # head(gene_ttaa_site_merge_plus_p_val)
        # full.genes.for.fisher <- gene_ttaa_site_merge_plus_p_val[,c(1,5)]

        #retrieve homolog IDs for all entrez_gene_ids for human and rat
        human.homologs <- getLDS(attributes = "entrezgene", mart = human, attributesL = "ensembl_gene_id", martL = mouse)
        rat.homologs <- getLDS(attributes = "entrezgene", mart = rat, attributesL = "ensembl_gene_id", martL = mouse)


        #subselect specific organisms for homolog mapping
        ratComplexes <- allComplexes[allComplexes$organism == "Rat",]
        humanComplexes <- allComplexes[allComplexes$organism == "Human",]
        mouseComplexes <- allComplexes[allComplexes$organism == "Mouse",]

        head(mouseComplexes$subunits.Entrez.IDs)

        #map mouse complexes
        mouse.trans <- lapply(as.character(mouseComplexes$subunits.Entrez.IDs), function(x) {
                mouse.list <- unlist(strsplit(x,","))
                unique(mouse.lookup$ensembl_gene_id[mouse.lookup$entrezgene %in% mouse.list])
        })

        mouseComplexes$subunits.Ensembl.IDs <- unlist(lapply(mouse.trans, function(x) {
                if (length(x) > 1){
                        return(paste(x,collapse=","))
                }
                return("")
        }))

        #get rid of all empty or 1 gene complexes
        mouseComplexes <- mouseComplexes[nchar(mouseComplexes$subunits.Ensembl.IDs) != 0,]

        #map rat complexes
        rat.homs <- lapply(as.character(ratComplexes$subunits.Entrez.IDs), function(x) {
                rat.list <- unlist(strsplit(x,","))
                unique(rat.homologs$Ensembl.Gene.ID[rat.homologs$EntrezGene.ID %in% rat.list])
        })

        ratComplexes$subunits.Ensembl.IDs <- unlist(lapply(rat.homs, function(x) {
                if (length(x) > 1){
                        return(paste(x,collapse=","))
                }
                return("")
        }))

        #get rid of all empty or 1 gene complexes
        ratComplexes <- ratComplexes[nchar(ratComplexes$subunits.Ensembl.IDs) != 0,]

        #map human complexes
        human.homs <- lapply(as.character(humanComplexes$subunits.Entrez.IDs), function(x) {
                human.list <- unlist(strsplit(x,","))
                unique(human.homologs$Ensembl.Gene.ID[human.homologs$EntrezGene.ID %in% human.list])
        })

        humanComplexes$subunits.Ensembl.IDs <- unlist(lapply(human.homs, function(x) {
                if (length(x) > 1){
                        return(paste(x,collapse=","))
                }
                return("")
        }))

        #get rid of all empty or 1 gene complexes
        humanComplexes <- humanComplexes[nchar(humanComplexes$subunits.Ensembl.IDs) != 0,]

        #merge all comples mouse plus homologs from human and rat
        mouseComplexes2 <-  rbind(mouseComplexes, ratComplexes, humanComplexes)

        #Add a source Column
        mouseComplexes2$Source <- rep("CORUM", nrow(mouseComplexes2))

        save(mouseComplexes2, file = "data/corum_complexes.RData")

        head(mouseComplexes2)
        colnames(mouseComplexes2)


        compleatComplexes <- read.delim("data2/compleatComplexes.txt", header=F, stringsAsFactors = F)


        compleatComplexes <- compleatComplexes[,c(1,5,6)]
        colnames(compleatComplexes) <- c("Complex.id", "subunits.Entrez.IDs", "Complex.name")

        head(compleatComplexes$subunits.Entrez.IDs)
        compleatComplexes$subunits.Entrez.IDs <- gsub(pattern = ";", replacement = ",", compleatComplexes$subunits.Entrez.IDs)


        #map human complexes
        human.homs2 <- lapply(as.character(compleatComplexes$subunits.Entrez.IDs), function(x) {
                human.list <- unlist(strsplit(x,","))
                unique(human.homologs$Ensembl.Gene.ID[human.homologs$EntrezGene.ID %in% human.list])
        })

        compleatComplexes$subunits.Ensembl.IDs <- unlist(lapply(human.homs2, function(x) {
                if (length(x) > 1){
                        return(paste(x,collapse=","))
                }
                return("")
        }))

        #get rid of all empty or 1 gene complexes
        compleatComplexes <- compleatComplexes[nchar(compleatComplexes$subunits.Ensembl.IDs) != 0,]

        compleatComplexes$Source <- rep("COMPLEAT", nrow(compleatComplexes))

        # save(compleatComplexes, file="data/compleatComplexes.RData")

        complexResults <- rbind(mouseComplexes2[,c(1,14)], compleatComplexes[,c(1,4)])


        complexResults2 <- do.call("rbind",apply(complexResults, 1, function(comp){
                cbind(comp[1], unlist(strsplit(comp[2],",")))
        }))

        head(complexResults2)
        colnames(complexResults2) <- c("ComplexID", "EnsemblID")
        complexResults2 <- as.data.frame(complexResults2)

        geneAggregate <- aggregate(complexResults2$ComplexID, by = list(complexResults2$EnsemblID), FUN=paste,collapse=",")
        colnames(geneAggregate) <- c("EnsemblID", "ComplexID")

        complexLookup <- list(complexResults, geneAggregate)
        save(complexLookup, file = "data/complexLookup.RData")
        return(list(complexResults, geneAggregate))
}
