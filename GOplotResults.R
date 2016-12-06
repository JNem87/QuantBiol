GOplotResults <-  function(GOresults, main = "Significant Go Terms Cluster", test = "Fisher elim Test"){
        ggplot( GOresults[[2]], aes(x=pval, y=factor(PATHNAME, levels = GOresults[[2]]$PATHNAME), size=as.factor(hits), color=cov)) +
                geom_point() +
                guides(fill=guide_legend(title = "test")) +
                scale_colour_gradient(low="red", high="blue", name="Coverage", guide = "colourbar") +
                scale_size_discrete(name="# of hits") +
                theme_bw() +
                xlab("-log10(P-Value)") +
                ylab(paste0("Significant Go Terms - ",test)) +
                ggtitle(main)
}
