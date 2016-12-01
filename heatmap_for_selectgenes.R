library(plyr)
library(dplyr)
library(magrittr)

# ======= Read DE gene list  =======

## rawdata is the log2(read counts per million), if need to normalize values to certain control, make sure to use 'substract'.

# project_folder <- 'C:/Users/Wen.Zhang/OneDrive for Business/MyProjects/RNAseq_Mike_mRNATransfection'
project_folder <- 'y:/BlueBirdBio/OneDrive for Business/MyProjects/RNAseq_Mike_mRNATransfection'
cpm <- read.table(file=paste(project_folder,'/EdgeR_EDAanalysis/2016-10-19_cmp.csv',sep=''),header = TRUE,sep=',')
# colnames(cpm)[1] <- 'genes'
cpm %<>%
    mutate(NoTransfect=rowMeans(.[,c(2,3)])) %>%
    mutate(LNP_Only=rowMeans(.[,c(4,5,6)])) %>%
    mutate(PolyIC=rowMeans(.[,c(7,8,9)])) %>%
    mutate(Tri_CCR5=rowMeans(.[,c(10,11,12)])) %>%    
    mutate(HPLC_CCR5=rowMeans(.[,c(13,14,15)])) %>%
    mutate(J2_CCR5=rowMeans(.[,c(16,17)]))
    
colnames(cpm)[1] <- 'genes'    
rownames(cpm) <- cpm$genes
write.csv(cpm,file=paste(project_folder,"/",as.character(Sys.Date()),"_cpm_average.csv",sep=''))


# ====== read DE genes and genes of interest from nanostring ====


library(readxl) # readxl read both modern xlsx and xls files
DEGenes <- read.table(paste(project_folder,'/EdgeR_DEanalysis/2016-10-27_DEgenes_Tri_CCR5_vs_LNP_Only.csv', sep=''),header = TRUE,sep=',')
topgenes <- as.character(DEGenes$genes[1:25])
# 
# Immunology <- read_excel(paste(project_folder,'/sample_info/Nanostring nCounter_Human_Immunology_v2_Panel_Gene_List.xlsx',sep=''),skip=1, sheet='Human Immunology Panel')
# Inflammation <- read_excel(paste(project_folder,'/sample_info/Nanostring LBL-C0269-01_nCounter_Human_Inflammation_v2_Panel_Gene_List.xlsx',sep=''),skip=1,sheet='Human Inflammation Panel')
targeted_path <- read_excel(paste0(project_folder,'/sample_info/targeted_pathways.xlsx'))

pathways <- c('Immunology','Inflammation','RIG-I-like receptor signaling',
              'Cytosolic DNA-sensing', 'TCL signaling' ,'Nucleic Acid Sensors',
              'Anti Viral Response','MYD88-Independent', 'MYD88-Dependent')

# split the df by column



split_df<- function(df,pathways)
{
    n=1
    m=2
    path_dfs = list()
    for (i in 1:length(pathways))
    {
    path <- pathways[i]
    path_dfs[[i]] <- df[,n:m]
    names(path_dfs)[i] = pathways [i]
    path_dfs[[i]] <- path_dfs[[i]][complete.cases(path_dfs[[i]][,1]),]
    path_dfs[[i]]$gene_cat <- pathways[i]
    # path_dfs[[i]]$gene_cat = as.factor(path_dfs[[i]]$gene_cat)
    names(path_dfs[[i]])[1] <- 'genes'
    path_dfs[[i]][is.na(path_dfs[[i]])] <- ''
    n = n+2
    m = m+2  
      
    }
# names(path_dfs) = pathways      
return(path_dfs)
}


pathways_dfs<- split_df(targeted_path,pathways)


# head(Immunology)
# head(Inflammation)
# 
# Immunology$gene_cat <- 'Immunology'
# Inflammation$gene_cat <- 'Inflammation'
# 
# Immunology_genes <- Immunology[,c("Official Symbol",'gene_cat','Alias / Prev Symbol')]
# names(Immunology_genes) <- c("genes",'gene_cat','alias')
# 
# Inflammation_genes <- Inflammation[,c("Official Symbol",'gene_cat','Alias / Prev Symbol')]
# names(Inflammation_genes) <- c("genes",'gene_cat','alias')



# =========== split gene alias and merge ==========

# # note, there are multiple alias for gene names, we need to split into seperate column 

require(stringr)
split_alias <- function(df)
{
    # df <- pathways_dfs[['Immunology']]
    gene_alias <- str_split(string=df$alias, pattern = ',')
    head(gene_alias)
    alias_matrix <- data.frame(Reduce(rbind, gene_alias))
    names(alias_matrix) <- paste('alias',1:ncol(alias_matrix))
    df_cb <- cbind(df,alias_matrix)
    return(list(alias_matrix,df_cb))
}

# merge with DE list, recursively, for all gene alias columns
merge_genes <- function(df1,df2,df2_aliasMx)
{
    geneNames =  c('genes',paste('alias',1:ncol(df2_aliasMx)))
    gene.merged <- merge(x=df2,y=df1, by.x='genes',by.y='genes')
    df.merge <- gene.merged[0,]
    
    for(gene in geneNames)
    {
        merged <- merge(x=df2,y=df1, by.x=gene,by.y='genes')
        df.merge <- rbind(df.merge, merged)
        return(df.merge)
    }
    
}


merge_alias <- function(pathways_dfs)
    
{
    pathways_dfs_merged <- list()
    
    for (i in 1:length(pathways_dfs))
    {
        # i=1
        df <- split_alias(pathways_dfs[[i]])
        pathways_dfs_merged[[i]] <- merge_genes(cpm,df[[2]],df[[1]])
        names(pathways_dfs_merged)[i] <-  names(pathways_dfs)[i]
    }  
    return(pathways_dfs_merged)
}


clean_df <- function(df)
{
    clean_df <- df[,c(1,3,(length(df)-5):length(df))]
    return(clean_df)
}


pathways_dfs_merged <- merge_alias(pathways_dfs)

clean_dfs <- lapply(pathways_dfs_merged,clean_df)


# Immunology_aliasMx<- split_alias(Immunology_genes)[[1]]
# Inflammation_aliasMx <- split_alias(Inflammation_genes)[[1]]
# Immunology_genes <- split_alias(Immunology_genes)[[2]]
# Inflammation_genes <- split_alias(Inflammation_genes)[[2]]
# 
# 
# immunology_merged <- merge_genes(cpm,Immunology_genes,Immunology_aliasMx)
# Inflammation_merged <- merge_genes(cpm,Inflammation_genes,Inflammation_aliasMx)
# cpm_topGenes <- cpm[topgenes,]
# cpm_topGenes$gene_cat <- 'Top Genes'


# immunology_clean <- immunology_merged %>% select(-c(3:26))
# inflammation_clean <- Inflammation_merged %>% select(-c(3:22))
# topGenes_clean <- cpm_topGenes %>% select(-c(2:17))
# topGenes_clean %<>% select(1,length(topGenes_clean),2:(length(topGenes_clean)-1))


# =================== combine and make tidy data =============

require(tidyr)

# df_keep <- immunology_clean[,c(1,2)]
# df_abovemean <- immunology_clean[,-c(1,2)] - rowMeans(immunology_clean[,-c(1,2)])
# df_tocontrol <- data_abovemean - data_abovemean$LNP_Only
# df_melt <- df_abovemean

get_longTable <- function(clean_dfs)
{
    df_melted <- list()
    for(i in 1:length(clean_dfs))
    {
        # i=1
        df = clean_dfs[[i]]
        name = names(clean_dfs)[i]
        df_keep <- df[,c(1,2)]
        df_tocontrol <- df[,-c(1,2)] - df$LNP_Only
        # df_abovemean <- df[,-c(1,2)] - rowMeans(df[,-c(1,2)])
        table_long <- paste0(name,'_tableLong') 
        table_long <- cbind(df_keep,df_tocontrol)
        table_long %<>% gather(conditions,cpm,NoTransfect:J2_CCR5)
        df_melted[[i]] <- table_long 
        names(df_melted)[i] <- name
    }
    
    return(df_melted)
    
}

# clean_dfs <- list(immunology=immunology_clean, inflammation = inflammation_clean,topgene= topGenes_clean)

df_longs <- get_longTable(clean_dfs)

# 
# immunology_long <- cbind(immunology_clean[,c(1,2)],(immunology_clean[,-c(1,2)] - rowMeans(immunology_clean[,-c(1,2)]) ))
# immunology_long %<>% gather(conditions,cpm,NoTransfect:J2_CCR5)
# 
# inflammation_long <- cbind(inflammation_clean[,c(1,2)],(inflammation_clean[,-c(1,2)] - rowMeans(inflammation_clean[,-c(1,2)]) ))
# inflammation_long %<>% gather(conditions,cpm,NoTransfect:J2_CCR5)
# 
# topGene_long <- cbind(topGenes_clean[,c(1,2)],(topGenes_clean[,-c(1,2)] - rowMeans(topGenes_clean[,-c(1,2)]) ))
# topGene_long %<>% gather(conditions,cpm,NoTransfect:J2_CCR5)



# refactoring the order of the conditions

reOrderFactors <- function(dfs)
{
    for(i in 1:length(dfs))
    {
        df <- dfs[[i]]
        df$conditions <- factor(df$conditions, levels=c("NoTransfect", "LNP_Only",
                                                        'HPLC_CCR5','J2_CCR5',"Tri_CCR5",'PolyIC'))
        dfs[[i]] <-df

    }
    
    return(dfs)
}

df_longs<- reOrderFactors(df_longs)


# immunology_long$conditions <- factor(immunology_long$conditions, levels=c("NoTransfect", "LNP_Only",
#                                                              'HPLC_CCR5','J2_CCR5',"Tri_CCR5",'PolyIC'))
# inflammation_long$conditions <- factor(inflammation_long$conditions, levels=c("NoTransfect", "LNP_Only",
#                                                              'HPLC_CCR5','J2_CCR5',"Tri_CCR5",'PolyIC'))
# immune_GOI$conditions <- factor(immune_GOI$conditions, levels=c("NoTransfect", "LNP_Only",
#                                                              'HPLC_CCR5','J2_CCR5',"Tri_CCR5",'PolyIC'))
# topGene_long$conditions <- factor(topGene_long$conditions, levels=c("NoTransfect", "LNP_Only",
#                                                              'HPLC_CCR5','J2_CCR5',"Tri_CCR5",'PolyIC'))

# immune_GOI <- rbind(df_longs$immunology,df_longs$inflammation)
# head(immune_GOI)
# tail(immune_GOI)
# topGene_long <- df_longs$topgene


# visualize data

require(ggplot2)
require(RColorBrewer)
require(ggthemes)
require(cowplot)
source('multiplot.R')

# heatmap for immune genes of interest
# dfs_forheatmap <- list(Immunology=df_longs$immunology,Inflammation=df_longs$inflammation)

## multiplot for individual gene list
heatmap <- function(dfs)
{
    require(ggplot2)
    plots = list()
    for(i in 1:length(dfs))
    {
        hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
        gg <- ggplot(dfs[[i]],aes(x=conditions,y=genes)) 
        gg <- gg + geom_tile(aes(fill=cpm)) +
            facet_wrap(~ gene_cat,ncol=2) +
            scale_fill_gradientn(colours = hm.palette(100),limits=c(-15, 15)) +
            theme_light()+
            theme(
                axis.line=element_blank(),
                axis.text.x=element_text(angle = 45,hjust = 1),
                legend.title = element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()
                ) +
            labs( title = names(dfs)[i])
        plots[[i]] <- gg

        # ggsave(file=paste(names(dfs)[i],'_heatmap.pdf',sep=''),plot=gg)

    }
    pdf('pathways_heatmap.pdf')
    gg_multi <- multiplot(plotlist=plots,cols=length(dfs))
    dev.off()
    return(gg_multi)
}

map <- heatmap(df_longs)


# multifacet for gene list
data <- do.call('rbind',df_longs)
rownames(data) <- NULL

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
gg <- ggplot(data,aes(x=conditions,y=genes)) 
gg <- gg + geom_tile(aes(fill=cpm)) +
    facet_wrap(~ gene_cat,ncol=2,scale='free') +
    scale_fill_gradientn(colours = hm.palette(100),limits=c(-15, 15)) +
    theme_light()+
    theme(
        axis.line=element_blank(),
        axis.text.x=element_text(angle = 45,hjust = 1),
        legend.title = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(1,0.5,0.5,1),'cm')
    ) +
    labs( title = 'Pathways Heatmaps')
gg
ggsave(file='heatmaps.pdf',plot=gg,width = 9, height = 16) 


## combine all the plots using plot_grid or multiplot

# pdf('heatmap.pdf')
# plot_grid(gg_immune, gg_TopGenes, labels=c("A", "B"), ncol = 2, nrow = 1) +
#     theme(plot.margin = unit(c(1,1,0.5,1),'cm'))
# dev.off()

source('multiplot.R')
pdf('heatmap.pdf')
multiplot(gg_immune, gg_TopGenes,cols=2)
dev.off()

