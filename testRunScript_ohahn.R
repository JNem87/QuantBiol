source("R/getGOMatrix.R")
source("R/getGeneMatrix.R")
source("R/buildDataMatrix.R")
source("R/ultimateEnrich.R")
source("R/buildTwobyTwoMatrix.R")
source("R/buildTwobyTwoMatrix2.R")
source("R/fisher_enrich.test.R")
source("R/fisher_enrich.elim.test.R")
source("R/transformtosetlist.R")
source("R/eliminatefrom2by2matrix.R")
source("R/eliminateFromSetMatrix.R")
source("R/getPathwayMatrix.R")
source("R/getComplexMatrix.R")
source("R/reduceToUniverse.R")
source("R/KS_enrich.test.R")
source("R/KS_enrich.elim.test.R")
source("R/wks.r")
require(Mus.musculus)
require(hash)
require(reactome.db)
require(hom.Hs.inp.db)
require(org.Hs.eg.db)
require(biomaRt)
require(RJSONIO)
require(reshape2)
require(ggplot2)
require(DOSE)
#require(profvis)

#test dataset
myTestData <- read.delim("~/Desktop/Oliver/Commonlist_genes.txt", header = F)
myTestData1 <- myTestData[,1, drop=FALSE]
head(myTestData)
#myHitList <- rownames(myTestData1[myTestData1 < 1e-20,,drop=F])
#hitList <- rownames(myTestData1[myTestData1 < 1e-20,,drop=F])

#myTestData1 <- myTestData[order(myTestData[,1]),,drop=F]

geneUniv <- read.delim("~/Desktop/Oliver/expressed_list.txt", header = F)
geneUniv <- as.character(geneUniv[,1])

#myTestData1 <- myTestData1[,1]

myTestData1 <- data.frame(rep(0.05, nrow(myTestData)), row.names = myTestData[,1])
rownames(myTestData1) <- myTestData[,1]
nrow(myTestData1)
head(myTestData1)
rownames(myTestData1)
myTestData1 <- -log10(myTestData1) +1
colnames(pathwayLookup[[1]])
# system.time(tmpRes <- ultimateEnrich(myTestData1,test = "ks", threshold = 0.05, type = "PATHWAY", maxSetSize=200))
# system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.01, type = "GO", maxSetSize=200, minHitNum = 3))
# system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.05, type = "PATHWAY", maxSetSize=200))
# system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher", threshold = 0.01, type = "GO", maxSetSize=200, minHitNum = 3))
# system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.01, type = "COMPLEX", maxSetSize=200))
#

system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher", threshold = 0.01, type = "GO", maxSetSize=200, minHitNum = 3))
system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher", threshold = 0.05, type = "PATHWAY", maxSetSize=150, minHitNum = 3))
system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher", threshold = 0.05, type = "COMPLEX", maxSetSize=200, minHitNum = 3))

system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.01, type = "GO", maxSetSize=200, minHitNum = 3, GeneUniverse = geneUniv))
system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.05, type = "PATHWAY", maxSetSize=200, minHitNum = 3, GeneUniverse = geneUniv))
system.time(tmpRes <- ultimateEnrich(myTestData1,test = "fisher_elim", threshold = 0.05, type = "COMPLEX", maxSetSize=200, minHitNum = 2, GeneUniverse = geneUniv))



system.time(tmpRes <- ultimateEnrich(myTestData1,test = "wks", threshold = 0.05, type = "PATHWAY", maxSetSize=200))


hist(unlist(tmpRes), breaks=200)
hist(unlist(tmpRes), breaks=200)


#load("data/pathwayLookup.RData")

FDRpass <- !unlist(tmpRes[[4]]) > 0.1
tmpRes[[2]] <- tmpRes[[2]][as.vector(which(FDRpass))]
tmpRes[[3]] <- tmpRes[[3]][as.vector(which(FDRpass))]
tmpRes[[4]] <- tmpRes[[4]][as.vector(which(FDRpass))]

tmpResTable <- unlist(tmpRes[[2]])
tmpResTable2 <- sapply(tmpRes[[3]], function(x){
        c(x[1,1], sum(x[1,1], x[1,2]), x[1,1] / sum(x[1,1], x[1,2]))
})

tmpResTab <- cbind(names(tmpResTable), tmpResTable, t(tmpResTable2)[names(tmpResTable),])

colnames(tmpResTab) <- c("PATHID", "pval", "hits", "total", "cov")

tmpResTab <- as.data.frame(tmpResTab)


#saveRDS(unlist(tmpRes), "gogoKonstantina_KSResults.rds")

#head(pathwayLookup[[1]])

head(sort(unlist(tmpRes)))

allComplexes <- read.csv2(file = "data2/allComplexes.csv")

compleatComplexes <- read.delim("data2/compleatComplexes.txt", header=F, stringsAsFactors = F)

annoComplex <- allComplexes[allComplexes$Complex.id %in% names(unlist(tmpRes[[2]])),c(1,2)]
annoComplex$Complex.name <- paste0(annoComplex$Complex.name, "(CORUM)")

annoComplex2 <- compleatComplexes[compleatComplexes$V1 %in% names(sort(unlist(tmpRes))),c("V1", "V6")]
colnames(annoComplex2) <- colnames(annoComplex)
annoComplex2$Complex.name <- substr(annoComplex2$Complex.name,1,100)
annoComplex2$Complex.name <- paste0(annoComplex2$Complex.name, "(COMPLEAT)")

annoComplex <- rbind(annoComplex,annoComplex2)

colnames(annoComplex) <- c("PATHID", "PATHNAME")

tmpResDFmerge <- merge(tmpResTab, annoComplex)

tmpResDFmerge <- tmpResDFmerge[!tmpResDFmerge$PATHNAME == "Unknown(COMPLEAT)",]

# compleatComplexes[compleatComplexes$V1 %in% names(sort(unlist(tmpRes))),c("V1", "V6")]
# substr(compleatComplexes[compleatComplexes$V1 %in% names(sort(unlist(tmpRes))),c("V1", "V6")]$V6,1,100)

compleatComplexes[compleatComplexes$V1 %in% "HC9681",]
substr(compleatComplexes[compleatComplexes$V1 %in% names(sort(unlist(tmpRes))),c("V1", "V6")]$V6,1,100)


complexLookup[[1]][complexLookup[[1]]$Complex.id ==  "HC9681",]
complexLookup[[1]]


#GOTEST
f <- data.frame(set=names(tmpRes[[2]]), pvalue=unlist(tmpRes[[2]]))

f[order(f$pvalue, decreasing = T),]

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

head(res[order(res$pvalue, decreasing = F),],20)

#################################################################
#PATHWAY results
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

# ggplot(tmpResDFmerge, aes(x=cov, y=factor(PATHNAME, levels = tmpResTab$PATHNAME), size=hits, color=pval)) +
#         geom_point() + scale_fill_gradient(low="red", high="blue") +
#         theme_bw()
#

load("data/goLookup.RData")
myTerms <- goLookup[[1]]
#
load("data/pathwayLookup.RData")
myTerms <- pathwayLookup[[1]]
#
load("data/complexLookup.RData")
myTerms <- complexLookup[[1]]

head(myTerms)
colnames(myTerms)[1] <- "Complex.id"

hitGenesinSets <- lapply(tmpResDFmerge$PATHID, function(x){
        tmpSet <- unlist(strsplit(myTerms[myTerms$Complex.id == x,2],","))
        hitSet <- tmpSet[tmpSet %in% rownames(myTestData1)]
        list(hitSet, tmpSet)
})

names(hitGenesinSets) <- tmpResDFmerge$PATHID

hitGenesinSetStrings <- lapply(hitGenesinSets, function(x){
        c(paste(x[[1]], collapse = ","), paste(x[[2]], collapse = ","))
})

hitGenesinSetTab <- do.call("rbind", hitGenesinSetStrings)

colnames(hitGenesinSetTab) <- c("hitGenes", "SetGenes")

hitGenesinSetMerge <- merge(tmpResDFmerge, hitGenesinSetTab, by.x="PATHID", by.y="row.names")

write.table(hitGenesinSetMerge[order(hitGenesinSetMerge$pval, decreasing = T),], "~/Desktop/Oliver/GOElim.tsv", sep = "\t", row.names = F, quote = F)
write.table(hitGenesinSetMerge[order(hitGenesinSetMerge$pval, decreasing = T),], "~/Desktop/Oliver/PathwayElim.tsv", sep = "\t", row.names = F, quote = F)
write.table(hitGenesinSetMerge[order(hitGenesinSetMerge$pval, decreasing = T),], "~/Desktop/Oliver/ComplexElim.tsv", sep = "\t", row.names = F, quote = F)


# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/GO_Fisher.png", width = 20, height = 15, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/GO_Fisher.pdf", width = 10, height = 15)
#pdf("~/Desktop/Oliver/GOElim.pdf", width = 10, height = 15)
ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHNAME, levels = tmpResDFmerge$PATHNAME), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Go Terms - Fisher Test") +
        ggtitle("Significant Go Terms")
dev.off()

# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/GO_Fisher_results.csv", quote=F, sep ="\t", row.names = F)
#
# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/GO_Fisher_Elim.png", width = 15, height = 10, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/GO_Fisher_Elim.pdf", width = 10, height = 7)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/GO_Fisher_Elim_crispr_24h_cluster_01.pdf", width = 10, height = 7)
pdf("~/Desktop/Oliver/GOElim_BG.pdf", width = 10, height = 15)
ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHNAME, levels = tmpResDFmerge$PATHNAME), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Go Terms - Fisher Elimination Test") +
        ggtitle("Significant Go Terms")
dev.off()
#
# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/GO_Fisher_Elim_results.csv", quote=F, sep ="\t", row.names = F)
#
# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/Pathway_Fisher.png", width = 15, height = 10, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Pathway_Fisher.pdf", width = 10, height = 15)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Pathway_Fisher_crispr_24h_cluster_01.pdf", width = 10, height = 15)
ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHNAME, levels = tmpResDFmerge$PATHNAME), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Pathways - Fisher Test") +
        ggtitle("Significant Pathways")
# dev.off()
#
# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/Pathway_Fisher_results.csv", quote=F, sep ="\t", row.names = F)
#
# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/Pathway_Fisher_Elim.png", width = 15, height = 10, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Pathway_Fisher_Elim.pdf", width = 10, height = 5)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Pathway_Fisher_Elim_crispr_24h_cluster_01.pdf", width = 10, height = 5)
pdf("~/Desktop/Oliver/PathwayElim_BG.pdf", width = 10, height = 10)
ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHNAME, levels = tmpResDFmerge$PATHNAME), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Pathways - Fisher Elimination Test") +
        ggtitle("Significant Pathways")
dev.off()
#
# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/Pathway_Fisher_Elim_results.csv", quote=F, sep ="\t", row.names = F)
#
# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/Complex_Fisher.png", width = 20, height = 15, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Complex_Fisher.pdf", width = 15, height = 15)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Complex_Fisher_crispr_24h_cluster_01.pdf", width = 7, height = 15)

ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHID, levels = tmpResDFmerge$PATHID), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Complexes - Fisher Test") +
        ggtitle("Significant Complexes")
# dev.off()
#
# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/Complex_Fisher_results.csv", quote=F, sep ="\t", row.names = F)
#
# png("/data/user/mgarmha1/Projects/Transposon/graphs/results/EnrichmentPlots/Complex_Fisher_Elim.png", width = 15, height = 10, units = "in", res = 130)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Complex_Fisher_Elim.pdf", width = 15, height = 8)
# pdf("/data/user/mgarmha1/Projects/Thesis/figures/Complex_Fisher_Elim_crispr_24h_cluster_01.pdf", width = 5, height = 5)
pdf("~/Desktop/Oliver/ComplexElim_BG.pdf", width = 15, height = 12)
ggplot(tmpResDFmerge, aes(x=pval, y=factor(PATHNAME, levels = tmpResDFmerge$PATHNAME), size=as.factor(hits), color=cov)) +
        geom_point() +
        guides(fill=guide_legend(title = "test")) +
        scale_colour_gradient(low="blue", high="red", name="Coverage", guide = "colourbar") +
        scale_size_discrete(name="# of hits") +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Significant Complexes - Fisher Elimination Test") +
        ggtitle("Significant Complexes")
dev.off()
#
# write.table(tmpResDFmerge, file = "../Transposon/results/ElEnR results/Complex_Fisher_Elim_results.csv", quote=F, sep="\t", row.names = F)


tmpResTab$pval <- -log2(as.numeric(as.character(tmpResTab$pval)))
tmpResTab$hits <- as.numeric(as.character(tmpResTab$hits))
tmpResTab$cov <- as.numeric(as.character(tmpResTab$cov))

tmpResTab <- tmpResTab[order(tmpResTab$pval),]

ggplot(tmpResTab, aes(x=pval, y=factor(PATHID, levels = tmpResTab$PATHID), size=hits, color=cov)) +
        geom_point() +
        scale_fill_gradient(low="blue", high="red") +
        theme_bw() +
        xlab("-log2(P-Value)") +
        ylab("Significant Pathways - Fisher Test")



colnames(tmpResDFmerge)
tmpResDFmergeMelt <- melt(tmpResDFmerge, id_vars=c("PATHID", "PATHNAME"))



##Learn hits

load("data/complexLookup.RData")
myGenes <- transformSettoList(complexLookup[[2]])
myTerms <- complexLookup[[1]]

hitGenesinSets <- lapply(tmpResDFmerge$PATHID, function(x){
        tmpSet <- unlist(strsplit(myTerms[myTerms$Complex.id == x,2],","))
        hitSet <- tmpSet[tmpSet %in% rownames(myTestData1)]
        list(hitSet, tmpSet)
})

names(hitGenesinSets) <- tmpResDFmerge$PATHID

hitGenesinSetStrings <- lapply(hitGenesinSets, function(x){
        c(paste(x[[1]], collapse = ","), paste(x[[2]], collapse = ","))
})

hitGenesinSetTab <- do.call("rbind", hitGenesinSetStrings)

colnames(hitGenesinSetTab) <- c("hitGenes", "SetGenes")

hitGenesinSetMerge <- merge(tmpResDFmerge, hitGenesinSetTab, by.x="PATHID", by.y="row.names")






####################




f <- data.frame(set=names(tmpRes), pvalue=unlist(tmpRes))

head(f)

f[order(f$pvalue, decreasing = T),]

cls <- columns(Mus.musculus)
cls <- cls[cls %in% c("GOID", "TERM")]


#Define GOID to be the keytype to use
kts <- keytypes(Mus.musculus)
kt <- kts[kts=="GOID"]

#get all keys
ks <- as.character(f[,1])
ks <- ks[1:10]

#query all ensembl_ids for all GOIDs
res <- select(Mus.musculus, keys=ks, columns=cls, keytype=kt)

res <- merge(res, f, by.x="GOID", by.y="set")

res[order(res$pvalue, decreasing = F),]





numhits <- buildTwobyTwoMatrix2(cpx_matrix = myMatrix, hitList = fisher_hitlist, myGeneUniverse = myGeneUniverse, minSetSize = 2, maxSetSize = 200)
