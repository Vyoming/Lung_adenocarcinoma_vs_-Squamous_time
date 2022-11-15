library("tidyverse")
library("AnnotationHub")
library("DESeq2")
library('apeglm')
library('reshape')
library('EnhancedVolcano')
library('fgsea')
library('biomaRt')


## Making a "short cut"

load("~/data/Bio470/deseq2.dds.RData")
dds
dds@colData$stage <- dds@colData$sample
dds@colData$stage <- factor(c("Early", "Late", "Late", "Late", "Early", "Late", "Early", "Late", 'Late', 'Early','Early','Early'))
dds@colData$stage<- factor(dds@colData$stage,levels = c('Early','Late'))

dds@colData$Cancer <- dds@colData$sample
dds@colData$Cancer <- factor(c("adenocarcinoma", "squamous", "adenocarcinoma", "squamous", "adenocarcinoma", "adenocarcinoma", "squamous", "squamous", 'adenocarcinoma', 'squamous','adenocarcinoma','squamous'))
dds@colData$Cancer<- factor(dds@colData$Cancer,levels = c('adenocarcinoma','squamous'))
dds@colData


dds@colData
dds$sample
levels(dds@colData$Intersect)
my_levels <- c("adenocarcinoma_Late", "squamous_Late", "adenocarcinoma_Early",   "squamous_Early"  )
dds@colData$Intersect <- my_levels
dds@colData$Intersect <- factor(dds@colData$Intersect,levels = my_levels)
design(dds) <- formula(~Cancer)
dds <- DESeq(dds)
resultsNames(dds)
dds
{
  results(dds)
  
  rownames(dds) <- gsub("\\.\\d*", "", rownames(dds))
  hub <- AnnotationHub()
  hubid <- "AH7799"
  anno <- hub[[hubid]]
  genemap <- tibble(gene_id=anno$gene_id,
                    symbol=anno$gene_name) %>%
    distinct()
  
  featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    dplyr::select(symbol)
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  
  dds <- nbinomWaldTest(dds)
  res <- results(dds, name='Cancer_squamous_vs_adenocarcinoma')
  res <- lfcShrink(dds, coef="Cancer_squamous_vs_adenocarcinoma", res=res, type='apeglm')
  
  res$gene_name <- mcols(dds)$symbol
  res <- res[order(res$log2FoldChange),]
  res_format <- res %>%
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    as_tibble() %>%
    rename_at(vars(-gene_id, -gene_name), ~ paste0(., ""))
  res_format[complete.cases(res_format),]
  res_format$log2FoldChange
  #res_format <- res_format[order(-res_format$log2FoldChange),]
  write.csv(res_format, 'Cancer_squamous_vs_adenocarcinoma_ranked.csv')
  
  
  Deseq_results <- mcols(dds) 
  Normalized_mean_counts <- assays(dds)[['counts']]
  
  Normalized_mean_counts <- assays(dds)[['rlog']]
  Normalized_mean_counts <- data.frame(Normalized_mean_counts)
  Normalized_mean_counts$Gene_Name <- mcols(dds)$symbol
  Normalized_mean_counts$Zscore <- mcols(dds)$Intersect_Naive_Experimental_vs_Naive_Control
  Normalized_mean_counts[complete.cases(Normalized_mean_counts),]
  write.csv(Normalized_mean_counts, 'Cancer_squamous_vs_adenocarcinoma_counts.csv')
}
sva_E
sva_L <- res_format
#rld <- rlog(dds, blind=FALSE)

#Normalized_mean_counts <- assays(rld)[['rlog']]
EnhancedVolcano(res_format, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = 'Squamous vs Adenocarcinoma', pCutoff = .05, FCcutoff = 1, xlim = c(-7,12), ylim = c(0,10) )
tail(res_format)

EnhancedVolcano(sva_E, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = 'Early squamous vs adenocarcinoma', pCutoff = .05, FCcutoff = 1, xlim = c(-5,7.5), ylim = c(0,6) )
EnhancedVolcano(sva_L, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = 'Late squamous vs adenocarcinoma', pCutoff = .05, FCcutoff = 1, xlim = c(-3,8), ylim = c(0,6) )


Heatmap_data <- data.frame(res_format)
MHC_Genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1')
Klr_Genes <-  c('Klra1','Klra4','Klrk1')
Metabolism_Genes <-  c('Ctsh','Ctse','Lyz2')
Hist_Genes <- c('Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
All_genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1','Klra1','Klra4','Klrk1', 'Ctsh','Ctse','Lyz2','Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
rownames(Heatmap_data) <-  Heatmap_data$gene_name
Heatmap_data <- Heatmap_data[Heatmap_data$gene_name %in% All_genes,]
Heatmap_data <- Heatmap_data[,c('log2FoldChange', 'gene_name')]

GOBP_CHRONIC_INFLAMMATORY_RESPONSE

Heatmap_data <-melt(Heatmap_data)
pathwayColors =rev(viridis::magma(10))
All_genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1','Klra1','Klra4','Klrk1', 'Ctsh','Ctse','Lyz2','Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
Heatmap_data$gene_name <- factor(Heatmap_data$gene_name,levels = All_genes)
ggplot(Heatmap_data, aes(gene_name, variable, fill= value)) + geom_tile() + scale_fill_gradientn(name = "Log2FC", colors = pathwayColors) + theme(panel.background = element_blank())
ggsave(file="Bulk_logfc_scale.pdf",  width=15, height=2, units="in")

#GSEA
res_format
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()
bm

res1 <- inner_join(res_format, bm, by=c("gene_id"="ensembl_gene_id"))
res2 <- res1 %>% 
  dplyr::select(hsapiens_homolog_associated_gene_name, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  summarize(log2FoldChange=mean(log2FoldChange))
res2

ranks <- deframe(res2)
head(ranks, 20)
pathways.hallmark <- gmtPathways("Documents/c5.all.v7.5.1.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgseaResTidier <- fgseaResTidy[which(fgseaResTidy$padj < .05 ),]
ggplot(fgseaResTidier, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


# heatmaps for genes of intrest
gene.list.activated <- c('Cdk5rap1','Lta','H2-D1','Plcb4','Ephb6')
gene.list.naive <- c('Ddr1','H2-T10','H2-D1','Rgs11','Dnahc8')
VP16_bulk <- res_format[res_format$gene_name %in% gene.list.activated,]
VP16_bulk <- VP16_bulk[c('log2FoldChange', 'gene_name' )]
VP16_bulk <-melt(VP16_bulk)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
ggplot(VP16_bulk, aes(gene_name, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.49,1.49)) + coord_flip() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(file="steady_activated_CVP_gene_heatmap.pdf",  width=2, height=3, units="in")

GOMF_MHC_CLASS_I_RECEPTOR_ACTIVITY	
GOBP_CHRONIC_INFLAMMATORY_RESPONSE
GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY
GOMF_CYTOKINE_ACTIVITY
GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION

GEO_MARKERS <- fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj)
Filter_GEO_MARKERS <-  GEO_MARKERS[GEO_MARKERS$pathway %in% c('GOBP_CHRONIC_INFLAMMATORY_RESPONSE', 'GOMF_MHC_CLASS_I_RECEPTOR_ACTIVITY', 'GOMF_CYTOKINE_ACTIVITY',
                                                              'GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION', 'GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY'),]

Filter_GEO_MARKERS$pathway <- factor(Filter_GEO_MARKERS$pathway,levels = c('GOBP_CHRONIC_INFLAMMATORY_RESPONSE', 'GOMF_MHC_CLASS_I_RECEPTOR_ACTIVITY', 'GOMF_CYTOKINE_ACTIVITY',
                                                                           'GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION', 'GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY'))

Filter_GEO_MARKERS <-  GEO_MARKERS[GEO_MARKERS$pathway %in% c('GOBP_T_CELL_MEDIATED_CYTOTOXICITY', 'GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION', 'GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY',
                                                              'GOMF_MHC_CLASS_II_RECEPTOR_ACTIVITY', 'GOBP_REGULATION_OF_FATTY_ACID_OXIDATION'),]

Filter_GEO_MARKERS$pathway <- factor(Filter_GEO_MARKERS$pathway,levels = c('GOBP_T_CELL_MEDIATED_CYTOTOXICITY', 'GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION', 'GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY',
                                                                           'GOMF_MHC_CLASS_II_RECEPTOR_ACTIVITY', 'GOBP_REGULATION_OF_FATTY_ACID_OXIDATION'))



pathwayColors =rev(viridis::plasma(10))
pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")

ggplot(data = Filter_GEO_MARKERS,aes(x=pathway,y=NES )) + theme_cem +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 5), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = pval))+
  geom_point(color = 'black' ,size = 1.4 ) +
  geom_point(aes(color = pval) ,size = 1 ) +
  scale_color_gradientn(colors=pathwayColors, breaks = c(.05, .01, .001, .00001),labels =  c('.05', '.01', '.001', '<.0001'),limits = c(0,.05), trans = scales::boxcox_trans(.25)) +
  scale_y_continuous(limits = c(-2.1,2),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  coord_flip() +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
ggsave(file="GSEA_Bubble_naive_EVC.pdf",  width=5, height=3, units="in")


#heiarchical clustering
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$stage, vsd$Cancer, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

