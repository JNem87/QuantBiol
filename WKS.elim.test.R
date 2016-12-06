#----------------------------------------------------------------
#               Weighted Kolmogorov-Smirnov test
#                    K. Charmpi, B. Ycart,
#               Version 1.0, February 19, 2015
#----------------------------------------------------------------

#---------------------- Variables -------------------------------
#   nsim          integer, number of simulations
#   ndiscr        integer, number of discretization points
#   w             named vector of numeric. It contains the
#                 weights, with the corresponding genes being
#                 given as names
#   weights       named vector of numeric, containing the weights
#                 ordered in decreasing order or a continuous function of it
#   gene.list     vector of character, containing the genes
#                 (names of weights)
#   db            list of vector of character, database of gene sets
#   gene.set      vector of character, gene set
#   gene.sets     vector of character, gene set names
#   ranked        logical, indicates whether the true values in w should
#                 be replaced by their ranks or not
#   alternative   character in "two.sided", "greater", "less". Specifies
#                 the alternative hypothesis of the WKS test
#   s             vector of numeric, a sample
#   pv1,pv2       named vectors of numeric, p-values
#   kwd           character. A keyword, usually a tissue
#                 (e.g. "liver","breast","lung", etc.)
#   nid           integer, number of points to be identified on a graph
#   v1,v2,pv      vectors of character
#   lens          list of numeric, containing the gene set sizes
#   title         character, title of a graphic
#   xlab          character, label for the x-axis of a graphic
#   ylab          character, label for the y-axis of a graphic

#---------------------- Functions -------------------------------

#------------------- Main functions------------------------------
# Statistical tests
# WKS.test(w,db,nsim=10000,ndiscr,ranked=FALSE,alternative="two.sided")
# GSEA.test(w,db,nsim=1000)
# Graphical comparison of two p-value vectors
# paired.pvalues(pv1,pv2,kwd,nid=0,title="Paired p-values",
#                  xlab="-log10 p-values",
#                  ylab="-log10 p-values")
# Graphical plot of gene sets
# cumulated.weights(w,db,gene.sets,title="Cumulated weights")
# enrichment.plot(w,db,gene.sets,xlab="genes",ylab="ES")

#----------------- Auxiliary functions--------------------------
# Calculation of the test statistics
# WKS.EnrichScore(gene.list, gene.set, weights,alternative)
# GSEA.EnrichScore(gene.list, gene.set, weights)
# Simulation for the limiting stochastic process
# lim.distr(weights,nsim,ndiscr,alternative)
# rbrmotint(weights,ndiscr)
# comparison of gene sets
# order.matching(v1,v2)
#--------------- global graphic parameters ---------------------
colpoint <- "blue3"                  # color for points
coltri <- "red3"                     # color for triangles
colline <- "black"                   # color for lines
colstep <- "blue3"                   # color for step function
colid <- "green4"                    # color for identifying
cexp <- 0.2                          # point size
lth <- 1.5                           # line thickness
lty <- 2                             # line type
cexa <- 1.1                          # axes size
cexl <- 1.1                          # labels size
cexm <- 1.3                          # title size
cext <- 0.7                          # text size
sctick <- 0.1                        # scale for ticks

#---------------------------------------------------------------
#                    Statistical tests
#---------------------------------------------------------------


#---------------------------------------------------------------
#                        WKS test
#---------------------------------------------------------------

WKS.elim.test <-
        function(w,
                 db,
                 nsim = 10000,
                 ndiscr,
                 ranked = FALSE,
                 alternative = "two.sided") {
                #   Takes a named vector w and a database db. Performs
                #   the WKS test with alternative "two.sided", "less", or "greater".
                #   The limiting distribution is estimated through nsim Monte-Carlo
                #   simulations, the number of discretization points being ndiscr.
                #   If ndiscr is missing, then the number of discretization points
                #   is equal to the size of w.
                #   Calculates the p-values (upper bound of the one-sided confidence
                #   interval with level 0.95), sorts them and returns them as a vector.
                #   If ranked=TRUE, the weights are  replaced by their ranks, and the
                #   pre-calculated F1.rds and F1os.rds values for the two sided and
                #   one sided WKS test respectively, are used.
                #
                #   Usage: l <- readRDS("liver.rds")
                #          C2 <- readRDS("C2.rds")
                #          pvwks <- WKS.test(l,C2);
                #          pvwks <- WKS.test(l,C2,alternative="greater");
                #          pvwks <- WKS.test(l,C2,ranked=TRUE)
                #          pvwks <- WKS.test(l,C2,ranked=TRUE,alternative="less")
                #
                correl.vector <-
                        sort(w, decreasing = TRUE)# sort the vector in decreasing order
                lg <- length(correl.vector)
                gene.list <- names(correl.vector)       # all genes
                #restructure set matrix

                db2 <- lapply(db, function(x){
                        rownames(db)[x]
                })

                db <- db2

                if (ranked) {
                        if (alternative == "two.sided")
                                ld <- readRDS("F1.rds")
                        else
                                ld <- readRDS("F1os.rds")
                        correl.vector <-
                                rank(correl.vector) / lg # replace true values by ranks
                        ots <- {
                                unlist(lapply(db, function(x) {
                                        # observed test statistics
                                        WKS.EnrichScore(gene.list,
                                                        x,
                                                        correl.vector,
                                                        alternative)
                                }))
                        }
                        ordots <-
                                ots[order(ots, decreasing = TRUE)]  # ordered in decreasing order
                        ordots <- na.omit(ordots)               # eliminate NA values
                        # (corresponding to gene sets with
                        # an empty intersection with the list)
                        pvalsm <- {
                                unlist(lapply(ordots, function(x) {
                                        length(which(ld >= x))
                                }))
                        }
                        # upper bound of the one-sided confidence
                        pvalsci <-
                                unlist(lapply(pvalsm, function(x)
                                        # interval with level 0.95
                                {
                                        prop.test(x, 1000000, alternative = "less")$conf.int[2]
                                }))
                } else{
                        ots <- {
                                unlist(lapply(db, function(x) {
                                        # observed test statistics
                                        WKS.EnrichScore(gene.list,
                                                        x,
                                                        abs(correl.vector),
                                                        alternative)
                                }))
                        }
                        ordots <-
                                ots[order(ots, decreasing = TRUE)]  # ordered in decreasing order
                        ordots <- na.omit(ordots)               # eliminate NA values
                        # (corresponding to gene sets with
                        # an empty intersection with the list)
                        if (missing(ndiscr)) {
                                ld <- {
                                        lim.distr(
                                                weights = abs(correl.vector),
                                                nsim = nsim,
                                                ndiscr = lg,
                                                alternative
                                        )
                                } # estimation of limiting distribution
                        } else{
                                ld <- {
                                        lim.distr(
                                                weights = abs(correl.vector),
                                                nsim = nsim,
                                                ndiscr = ndiscr,
                                                alternative
                                        )
                                }
                        }
                        pvalsm <- {
                                unlist(lapply(ordots, function(x) {
                                        length(which(ld >= x))
                                }))
                        }
                        # upper bound one-sided confidence
                        pvalsci <-
                                unlist(lapply(pvalsm, function(x)
                                        # interval with level 0.95
                                {
                                        prop.test(x, nsim, alternative = "less")$conf.int[2]
                                }))
                }

                return(pvalsci)
        }                                       # end function WKS.test


#---------------------------------------------------------------
#                        GSEA test
#---------------------------------------------------------------


GSEA.test <- function(w, db, nsim = 1000) {
        #   Takes a named vector w and a database db. Performs the GSEA test.
        #   The number of Monte-Carlo simulations for every gene set
        #   is nsim. Calculates the p-values (upper bound of the one-sided confidence
        #   interval with level 0.95), sorts them and returns them as a vector.
        #
        #   Usage: dcol <- readRDS("liver.rds")
        #          C2 <- readRDS("C2.rds")
        #          pvgsea <- GSEA.test(dcol,C2)
        #
        weights <-
                sort(w, decreasing = TRUE) # sort weights in decreasing order
        gene.list <- names(weights)        # gene names
        db <-
                reduce.symbols(gene.list, db) # reduce database to symbols in gene.list
        Ng <- length(db)                   # number of gene sets
        Obs.ES <- vector(length = Ng, mode = "numeric")
        lg <- vector(length = Ng, mode = "numeric")
        lens <- lapply(db, length)
        for (i in 1:Ng) {
                Obs.ES[i] <- {
                        GSEA.EnrichScore(gene.list,
                                         db[[i]], weights)
                }
        }
        phi <- matrix(nrow = Ng, ncol = nsim)
        for (r in 1:nsim) {
                rdb <- random.database(lens, gene.list)# take random gene sets
                for (i in 1:Ng) {
                        # calculate test statistic
                        phi[i, r] <- {
                                GSEA.EnrichScore(gene.list,
                                                 rdb[[i]], weights)
                        }
                }
        }

        p.vals <- matrix(0, nrow = Ng, ncol = 1)

        for (i in 1:Ng) {
                ES.value <- Obs.ES[i]
                if (ES.value >= 0) {
                        # estimate significance according to the
                        # sign of the observed test statistic
                        temp <- phi[i, which(phi[i, ] >= 0)]
                        p.vals[i, 1] <- {
                                signif(length(which(
                                        temp >= ES.value
                                )) /
                                        length(temp),
                                digits = 5)
                        }
                } else {
                        temp <- phi[i, which(phi[i, ] < 0)]
                        p.vals[i, 1] <- {
                                signif(length(which(
                                        temp <= ES.value
                                )) /
                                        length(temp),
                                digits = 5)
                        }
                }
        }
        indna <- which(is.na(p.vals))             # replace the na p-values
        p.vals[indna] <- 0                        # with zero

        # upper bound of the one-sided confidence
        p.vals <-
                unlist(lapply(p.vals * nsim, function(x)
                        # interval with level 0.95
                {
                        prop.test(x, nsim, alternative = "less")$conf.int[2]
                }))
        names(p.vals) <- names(db)
        p.valss <- sort(p.vals)                   # sort p-values
        return(p.valss)
}                                         # end function GSEA.test

random.database <- function(lens, gene.list) {
        #   Returns random vectors of sizes given in list
        #   lens out of gene.list.
        #
        return(lapply(lens, function(n) {
                sample(gene.list, n)
        }))
}                                       # end function random.database



reduce.symbols <- function(pv, db) {
        #   Reduces database db to symbols present in vector pv.
        #   Returns the reduced database.
        #
        rdb <- {
                lapply(db, function(v) {
                        # apply to all gene sets
                        return(find.matching(v, pv))
                })
        }     # reduce gene set to symbols in pv
        l <- lapply(rdb, length)                 # new gene set lengths
        rdb <- rdb[which(l > 0)]                  # remove empty gene sets
        return(rdb)
}                                       # end function reduce.symbols

find.matching <- function(v1, v2) {
        #    returns the character chains common to
        #    v1 and v2, any two vectors of character chains
        #
        return(v2[which(v2 %in% v1)])
}                                      # end function find.matching

order.matching <- function(v1, v2) {
        #    returns the indices of those entries of v1 found in v2,
        #    v1 and v2 being any two vectors of charater chains
        #
        return(which(v1 %in% v2))
}                                      # end function order.matching

#---------------------------------------------------------------
#             Calculation of the test statistics
#---------------------------------------------------------------

#------------------ auxiliary functions --------------------------

normal.function <- function(weights) {
        #   Takes a vector weights which contains the discretized
        #   values of a function g. Approximates the integral and
        #   normalizes it so that it sums up to 1.
        #
        distr1 <-
                cumsum(weights)               # calculate (approximate) the integral and
        distr <-
                distr1 / distr1[length(weights)] # normalize it so that it sums up to 1
        return(distr)
}                                       # end function normal.function


calc_max_diff <- function(weights, s, alternative) {
        #   Takes a sample s and a vector weights.
        #   Calculates and returns the value of the test
        #   statistic used for the test with alternative
        #   hypothesis alternative.
        #
        ss <- sort(s)                           # sorted sample
        lg <- length(s)                         # sample size
        wss <- weights[ss]                      # corresponding weights
        rts <-
                cumsum(wss)                      # random value of the test statistic
        rtsn <- rts / rts[lg]                     # same normalized
        distr <- normal.function(weights)       # normalized integral
        trval <-
                distr[ss]                      # expected weighted distribution
        switch(
                alternative,
                "two.sided" = {
                        diff1 <- rtsn - trval      # difference at the discontinuity points
                        if (lg >= 2) {
                                diff2 <- {
                                        c(trval[1], trval[2:lg] - rtsn[1:(lg - 1)])
                                }
                                # and at one point before
                        } else{
                                diff2 <- trval[1]
                        }
                        diff <- abs(c(diff1, diff2))
                },
                "greater" = {
                        diff <- rtsn - trval
                },
                "less" = {
                        if (lg >= 2) {
                                diff <- {
                                        c(trval[1], trval[2:lg] - rtsn[1:(lg - 1)])
                                }
                                # at one point before
                        } else{
                                diff <- trval[1]
                        }
                }
        )
        m <- max(diff)                          # take the maximum
        return(m)
}                                       # end function calc_max_diff


WKS.EnrichScore <-
        function(gene.list, gene.set, weights, alternative) {
                #   Calculates the value of the WKS test statistic
                #   when the alternative hypothesis is alternative for
                #   given numeric data weights, indexed by the
                #   genes in gene.list and a given set of genes gene.set.
                #   Returns the calculated observed test statistic.
                #   Usage: dcol <- readRDS("kidney.rds")
                #          C2 <- readRDS("C2.rds")
                #          weights <- sort(dcol,decreasing=TRUE)
                #          gene.list <- names(weights)
                #          gene.set <- C2[[16]]
                #          WKS.EnrichScore(gene.list, gene.set, abs(weights),
                #                         alternative="two.sided")
                #          WKS.EnrichScore(gene.list, gene.set, abs(weights),
                #                         alternative="greater")
                #
                cgi <- order.matching(gene.list, gene.set)# common gene indices
                if (length(cgi) > 0)
                        m <- {
                                calc_max_diff(weights = weights,
                                              s = cgi,
                                              alternative)
                        }                      # calculate test statistic
                else
                        m <- NA
                return(m * sqrt(length(cgi)))             # scale it
        }                                       # end function WKS.EnrichScore


GSEA.EnrichScore <- function(gene.list, gene.set, weights) {
        #   Calculates the value of the GSEA test statistic for
        #   given numeric data weights, ranked in decreasing order, indexed by the
        #   genes in gene.list and a given set of genes gene.set.
        #   Returns the calculated observed test statistic.
        #   Part of the functions encoded by the Broad Institute, have been used.
        #   Usage: dcol <- readRDS("kidney.rds")
        #          C2 <- readRDS("C2.rds")
        #          weights <- sort(dcol,decreasing=TRUE)
        #          gene.list <- names(weights)
        #          gene.set <- C2[[16]]
        #          GSEA.EnrichScore(gene.list, gene.set, weights)
        #
        # get signs
        tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
        no.tag.indicator <- 1 - tag.indicator
        N <- length(gene.list)                # length of gene list
        Nh <- length(gene.set)                # gene set length
        Nm <-
                N - Nh                         # difference between the two
        correl.vector <- abs(weights)
        sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
        norm.tag    <- 1.0 / sum.correl.tag     # normalization factor
        norm.no.tag <- 1.0 / Nm                 # 1/(N-Nh)
        # calculation of the test statistic
        RES <- {
                cumsum(
                        tag.indicator * correl.vector * norm.tag
                        - no.tag.indicator * norm.no.tag
                )
        }
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
                ES <- signif(max.ES, digits = 5)
        } else {
                ES <- signif(min.ES, digits = 5)
        }
        return(ES)
}                                        # end function GSEA.EnrichScore


#----------------------------------------------------------------
#        Simulation for the limiting stochastic process
#----------------------------------------------------------------

#------------------ auxiliary function --------------------------

rbrmotint <- function(weights, ndiscr) {
        #   Takes a vector weights.
        #   Simulates a trajectory of the stochastic integral
        #   of weights with respect to BM in ndiscr
        #   discretization points. Returns the result.
        #
        incr <- rnorm(n = ndiscr, sd = 1 / sqrt(ndiscr))# BM increments
        lg <- length(weights)
        subseq <- trunc(seq(1, lg, lg / ndiscr))
        mult <- incr * weights[subseq]            # multiply increments
        res <- cumsum(mult)                     # sum of increments
        return(res)
}                                       # end function rbrmotint

#--------------------- main function ----------------------------

lim.distr <- function(weights, nsim, ndiscr, alternative) {
        #   Takes a vector weights. Simulates nsim trajectories
        #   of the limiting stochastic process (centered BM integral).
        #   The number of discretization points for each is ndiscr.
        #   Returns a vector with the realizations of the limiting
        #   random variable (maximum of the absolute simulated trajectory
        #   when alternative="two.sided" or maximum of the simulated
        #   trajectory else).
        #   Usage: dcol <- readRDS("kidney.rds")
        #          weights <- sort(dcol,decreasing=TRUE)
        #          lg <- length(weights)
        #          ld <- lim.distr(abs(weights),1000,lg,alternative="two.sided")
        #          plot(ecdf(ld))
        #          ld2 <- lim.distr(abs(weights),1000,lg,alternative="greater")
        #          lines(ecdf(ld2),col=2)
        #
        lg <- length(weights)
        subseq <- trunc(seq(1, lg, lg / ndiscr))
        distr1 <-
                cumsum(weights)               # approximate the integral and
        distr <-
                distr1 / distr1[lg]              # normalize it so that it sums up to 1
        norm_fact <- distr1[lg] / lg             # normalization factor

        fres <- vector(nsim, mode = 'numeric')     # vector initialization
        for (k in 1:nsim) {
                res <-
                        rbrmotint(weights, ndiscr)       # simulation of the BM integral
                # realization of process
                resnorm <- res - distr[subseq] * res[ndiscr]
                if (alternative == "two.sided") {
                        fres[k] <- max(abs(resnorm))         # maximum of the absolute limit
                } else {
                        fres[k] <- max(resnorm)              # asymptotic random variable
                }
        }                                       # end for
        return(fres / norm_fact)                  # scale the result
}                                       # end function lim.distr



#----------------------------------------------------------------
#                Plot two vectors of p-values
#----------------------------------------------------------------

paired.pvalues <-
        function(pv1,
                 pv2,
                 kwd,
                 nid = 0,
                 title = "paired p-values",
                 xlab = "-log10 p-values",
                 ylab = "-log10 p-values") {
                #   Takes two vectors of p-values pv1 and pv2.
                #   Plots the two vectors (pv1 on the y-axis, pv2 on the x-axis),
                #   -log10 transformed, and represents the gene sets
                #   matching kwd as red triangles. Identifies nid of those.
                #
                #   Usage: paired.pvalues(pvwks,pvgsea,"liver",title="WKS vs. GSEA",
                #            xlab="-log10 p-values GSEA",ylab="-log10 p-values WKS")
                #
                lpv1 <- -log10(pv1)                     # -log10 transformed
                lpv2 <- -log10(pv2[names(pv1)])         # -log10 transformed
                {
                        plot(
                                lpv1,
                                lpv2,
                                # all points
                                xlab = xlab,
                                ylab = ylab,
                                main = title,
                                cex.axis = cexa,
                                cex.lab = cexl,
                                cex.main = cexm,
                                pch = 19,
                                cex = cexp,
                                col = "blue4"
                        )
                }
                abline(0, 1, lty = lty, lwd = lth)             # add bisector
                abline(h = -log10(0.05),
                       lty = lty,
                       lwd = lth)  # horizontal line at pv=0.05
                abline(v = -log10(0.05),
                       lty = lty,
                       lwd = lth)  # vertical line at pv=0.05
                if (!missing(kwd)) {
                        lkwd <- tolower(kwd)                 # switch to lower case
                        nmpv <- tolower(names(pv1))          # switch to lower case
                        indt <-
                                which(grepl(lkwd, nmpv))      # get gene sets matching keyword
                        lpv1t <- lpv1[indt]                  # their p-values for test 1
                        lpv2t <- lpv2[indt]                  # their p-values for test 2
                        {
                                points(
                                        lpv1t,
                                        lpv2t,
                                        # points with biological info
                                        pch = 17,
                                        cex = 1.2,
                                        col = coltri
                                )
                        }       # as red triangles
                        if (nid > 0) {
                                {
                                        identify(
                                                lpv1t,
                                                lpv2t,
                                                n = nid,
                                                # identify nid gene sets
                                                labels = names(lpv1t),
                                                col = colid,
                                                cex = cext
                                        )
                                }
                        }
                }
        }                                       # end function paired.pvalues

#----------------------------------------------------------------
#               Plot cumulated proportions of weights
#----------------------------------------------------------------

cumulated.weights <-
        function(w, db, gene.sets, title = "Cumulated weights") {
                #   Takes a named vector w, a database db and a vector of
                #   gene set names gene.sets.
                #   Plots the realizations of the random quantity used in
                #   WKS test statistic for the given gene sets.
                #   The bisector and the expected theoretical function are
                #   added on the graph.
                #
                #   Usage: dcol <- readRDS("kidney.rds")
                #          C2 <- readRDS("C2.rds");
                #          gene.sets <- c("ACEVEDO_LIVER_CANCER_DN","ACEVEDO_LIVER_CANCER_UP")
                #          cumulated.weights(dcol,C2,gene.sets)
                #
                xlab <- "t"
                ylab <- "Sn(t)"
                weights <-
                        abs(sort(w, decreasing = TRUE))        # sort weights in decreasing order
                gene.list <- names(weights)                    # all genes
                lgw <- length(weights)                         # their number
                stpfunc <- function(gene.set) {
                        ma <-
                                order.matching(gene.list, gene.set)    # indices of common genes
                        ma1 <- ma / lgw                               # proportion
                        vma <-
                                weights[ma]                          # corresponding weights
                        fvma <- cumsum(vma)                         # cumulated weights
                        lg <- length(fvma)
                        fvman <- fvma / fvma[lg]                      # normalize
                        sfun <-
                                stepfun(ma1, c(0, fvman), f = 0)      # make it a step function
                        return(sfun)
                }                                              # end step function
                sfun <-
                        stpfunc(db[[gene.sets[1]]])            # compute first step function
                {
                        plot(
                                sfun,
                                verticals = TRUE,
                                do.points = FALSE,
                                # plot it
                                xlim = c(0, 1),
                                ylim = c(0, 1),
                                col = colstep,
                                xlab = xlab,
                                ylab = ylab,
                                main = title,
                                cex.axis = cexa,
                                cex.main = cexm
                        )
                        }
                len <- length(gene.sets)                       # get number of plots
                if (len > 1) {
                        # if more than one
                        {
                                lapply(gene.sets[2:len],                     # iterate
                                       function(x) {
                                               lines(
                                                       stpfunc(db[[x]]),
                                                       xlim = c(0, 1),
                                                       verticals = TRUE,
                                                       do.points = FALSE,
                                                       col = colstep
                                               )
                                       })
                        }
                }
                abline(0,
                       1,
                       lty = lty,
                       lwd = lth,
                       col = colline)        # add bisector
                distr1 <-
                        cumsum(abs(weights))                 # calculate the integral and
                distr <- distr1 / distr1[lgw]                    # normalize
                {
                        lines(
                                seq(1 / lgw, 1, 1 / lgw),
                                distr,
                                lty = lty,
                                lwd = lth,
                                col = colline
                        )
                }               # add the expected theoretical
        }                                              # end function cumulated.weights


enrichment.plot <- function(w,
                            db,
                            gene.sets,
                            xlab = "genes",
                            ylab = "ES") {
        #   Takes a named vector w, a database db and a vector of
        #   gene set names gene.sets and creates a matrix of
        #   ceiling(length(gene.sets)/2) x 2 plots.
        #   For each gene set, the difference between the cumulated
        #   weights (S_n) and the primitive function (G) is plotted.
        #   At the bottom of the graph, vertical lines are added,
        #   indicating the positions of the genes inside the gene set
        #   in the ranked list.
        #
        #   Usage: dcol <- readRDS("kidney.rds")
        #          C2 <- readRDS("C2.rds");
        #          gene.sets <- c("ACEVEDO_LIVER_CANCER_DN","ACEVEDO_LIVER_CANCER_UP")
        #          enrichment.plot(dcol,C2,gene.sets)
        #          pvwks <- WKS.test(dcol,C2)
        #          gene.sets <- names(pvwks)[1:5]; enrichment.plot(dcol,C2,gene.sets)
        #          gene.sets <- names(tail(pvwks,5)); enrichment.plot(dcol,C2,gene.sets)
        #
        len <- length(gene.sets)                       # number of gene sets
        weights <-
                abs(sort(w, decreasing = TRUE))        # sort weights in decreasing order
        gene.list <- names(weights)                    # all genes
        lgw <- length(weights)                         # their number
        distr1 <-
                cumsum(abs(weights))                 # calculate the integral and
        distr <- distr1 / tail(distr1, 1)                 # normalize
        plotenrich <- function(gene.set, cexmain) {
                name <- gene.set
                gene.set <- db[[gene.set]]
                ma <-
                        order.matching(gene.list, gene.set)    # indices of common genes
                ma1 <- ma / lgw                               # proportion
                vma <-
                        weights[ma]                          # corresponding weights
                fvma <- cumsum(vma)                         # cumulated weights
                lg <- length(fvma)
                fvman <- fvma / tail(fvma, 1)                  # normalize
                trv <-
                        distr[ma]                            # (expected) theoretical values
                X <- rbind(ma1, ma1)
                X <-  as.vector(X)
                X <- c(0, X, 1)                               # ordinates
                Y <- rbind(fvman, fvman)
                Y <- as.vector(Y)
                Y <-
                        c(0, 0, Y)                               # coordinates for the step function
                C <-  rbind(trv, trv)
                C <- as.vector(C)
                C <- c(0, C, 1)
                Y <-
                        Y - C                                    # coordinates for the difference
                scy0 <- min(Y)
                scy1 <- scy0 + (max(Y) - scy0) * sctick
                scyM <- max(Y)
                # plot difference S_n-G
                {
                        plot(
                                X,
                                Y,
                                xlim = c(0, 1),
                                ylim = c(2 * scy0 - scy1, scyM),
                                type = "l",
                                lty = 1,
                                col = colstep,
                                xlab = xlab,
                                ylab = ylab,
                                main = name,
                                cex.lab = cexl,
                                cex.main = cexmain,
                                cex.axis = cexa
                        )
                }
                lines(c(0, 1),
                      c(0, 0),
                      lty = lty,
                      col = colline)    # add the line y=0
                scy0 <- min(Y)
                scy1 <- scy0 + (max(Y) - scy0) * sctick
                X <- rbind(ma1, ma1)
                Y <- rbind(rep(2 * scy0 - scy1, lg), rep(scy0, lg))
                matlines(X, Y, lty = rep(1, lg), col = colline)     # add the vertical lines
        }
        # split graphical window
        if (len > 1) {
                par(mfrow = c(ceiling(len / 2), 2))              # more than one gene set
                par(mar = c(4, 4, 1.5, 1.5))                   # margin parameters
                # apply to all gene sets
                pl <- lapply(gene.sets, function(x) {
                        plotenrich(x, cext)
                })
        } else{
                # only one gene set
                plotenrich(gene.sets, cexm)
        }                                           # end if
}                                              # end function enrichment.plot
