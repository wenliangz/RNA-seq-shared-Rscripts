plotmultigsea <- function(ranks)
{
    library(fgsea)
    library(ggplot2)
    fgseaRes <- fgsea(pathways = pathways, 
                      stats = rnks[[1]],
                      minSize=15,
                      maxSize=500,
                      nperm=10000)
    
    dim(fgseaRes)
    pathway_top <- fgseaRes[order(padj), ][2:3,]$pathway  # pick the no.2&3 pathway manually after review
    
    plots <- rep(list(list()),length(ranks))
    
    for (i in 1:length(ranks))
    {
        
        for(j in 1:length(pathway_top))
        {
            plt <- plotEnrichment(pathways[[pathway_top[j]]],ranks[[i]]) + 
                labs(title=pathway_top[j]) +
                ylim(-.1,1)
            plots[[i]][[j]] <- plt
            plots[[i]][[j]]
            
        }
        # pdf('pathways_heatmap.pdf')
        # gg_multi <- multiplot(plotlist=plots,cols=2)
        # dev.off()
        
    }
    
    return(plots)
}
