library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
source('multiplot.R')

# read bingo data
path <- 'y:/BlueBirdBio/OneDrive for Business/MyProjects/RNAseq_Mike_mRNATransfection/pathway'
# path <- 'C:/Users/Wen.Zhang/OneDrive for Business/MyProjects/RNAseq_Mike_mRNATransfection/pathway/'
fls <- list.files(path=path,pattern="*bgo$")

read_bingo <- function(path,filename)
{
    df <- str_split(filename,pattern = '\\.')[1]
    df <- read.table(paste(path,filename, sep=''),header = TRUE,sep='\t',skip = 19,fill=TRUE,stringsAsFactors = FALSE) 
    df %<>% mutate(log10corr.p = -log10(corr.p.value))
    
    return(df)
}


bingos <- list()
for (i in 1:length(fls))
{
    # i=1
    bingos[[i]] <- read_bingo(path,fls[i])
    names(bingos)[i] <- str_split(fls[i],pattern = '\\.')[[1]][1]
}



plot_enrichhis <- function(dfs)
{
    
    plots <- list()
    for (i in 1:length(dfs))
    {
        # i=1
        # dfs = bingos
        df <- dfs[[i]]
        n = ifelse(nrow(df) >=10,10,nrow(df))
        df <- dfs[[i]][1:n,]
        df <- df[order(df$log10corr.p),]
        df$Description <- factor(df$Description,levels=df$Description)  
        gg <- ggplot(df[order(df$log10corr.p),],aes(x=Description,y=log10corr.p,fill=Description)) +
        geom_bar(stat='identity') +
        coord_flip() +
        theme(legend.position="none",
              axis.text= element_text(size=16),
              axis.title = element_text(size=16)) +
        labs(x='',y = '-log10corr.p',
             title = names(dfs)[i]
             )  
        gg
        plots[[i]] <- gg
    }
    # pdf('pathways_heatmap.pdf')
    gg_multi <- multiplot(plotlist=plots,cols=1)
    # dev.off()
    return(gg_multi)
}


plot_enrichhis(bingos)

