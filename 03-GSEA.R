#pseudobulk DE to compare MDS and HD ----

library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)

#Load data
dir <- dir
seu <- readRDS(paste0(dir,'seurat_combined.rds'))

#Normalize and scale RNA assay
DefaultAssay(seu) <- "RNA"
seu <- seu %>% NormalizeData() %>% ScaleData()
Idents(seu) <- "cell_annos"


#run pseudobulk
pseudo_seu <- AggregateExpression(seu, assays = "RNA", return.seurat = T, group.by = c("Condition", "Specimen_ID", "cell_annos"))
pseudo_seu@meta.data
pseudo_seu$cell_annos <- gsub(".*_","",rownames(pseudo_seu@meta.data))

tail(Cells(pseudo_seu))


pseudo_seu$celltype.Dataset <- paste(pseudo_seu$cell_annos, pseudo_seu$Condition, sep = "_")
unique(pseudo_seu$celltype.Dataset)
Idents(pseudo_seu) <- "celltype.Dataset"
cell_types <- unique(pseudo_seu$cell_annos)
#MDS vs. HD HSC-HSPCs
for (cell_type in cell_types){
  res <- FindMarkers(object = pseudo_seu, 
                     ident.1 = paste0(cell_type,"_MDS"), 
                     ident.2 = paste0(cell_type,"_HD"),
                     test.use = "DESeq2", min.pct = 0.4)
  res <- res %>%  dplyr::mutate(Threshold = ifelse(p_val_adj < 0.05 & avg_log2FC  > 1.5, "Upregulated_in_MDS", ifelse(p_val_adj < 0.05 & avg_log2FC  < -1.5, "Upregulated_in_Healthy", "FALSE"))) %>% dplyr::arrange(desc(Threshold))
  write.csv(res, paste0(dir,"DE/",cell_type,"_","MDS_vs_Healthy_","pseudobulk.csv"))
}




#GSEA ----

library(fgsea)
library(tidyverse)
library(DT)
library(dplyr)
library(reactome.db)

set.seed(100)
# Load the pathways into a named list
dir <- dir
pathways.hallmark <- gmtPathways(paste0(dir,"h.all.v2024.1.Hs.symbols.gmt"))
pathways.wiki <- gmtPathways(paste0(dir,"c2.cp.wikipathways.v2024.1.Hs.symbols.gmt"))
pathways.reactome <- gmtPathways(paste0(dir,"c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
pathways.kegg <- gmtPathways(paste0(dir,"c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"))
pathways.biocarta <- gmtPathways(paste0(dir,"c2.cp.biocarta.v2024.1.Hs.symbols.gmt"))

cell_de <- list.files(paste0(dir,'DE/'), recursive = T)


#GSEA function
runGSEA <- function(pathway, stats, nperm){
  fgseaRes <- fgsea(pathways=pathway, stats=stats, nperm=nperm)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy <- fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>%
    filter(padj < 0.05)
  
  return(fgseaResTidy)
}


#Calculate gene ranks
for(i in 1:length(cell_de)){
  
  
  cell_type <- gsub('_.*',"", cell_de[i])
  
  res <- read.csv(paste0(dir, 'DE/', cell_de[i]))
  
  res2 <- res %>%
    dplyr::select(X, avg_log2FC) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(X) %>% 
    summarize(stat=mean(avg_log2FC))
  res2
  ranks <- deframe(res2)

  
  #Calulate NES ----
  hallmark <- runGSEA(pathways.hallmark, ranks, 1000)
  wiki <- runGSEA(pathways.wiki, ranks, 1000)
  reactome <- runGSEA(pathways.reactome, ranks, 1000)
  kegg <- runGSEA(pathways.kegg, ranks, 1000)
  biocarta <- runGSEA(pathways.biocarta, ranks, 1000)
  
  
  GSEA_combined <- as.data.frame(rbind(hallmark, wiki, reactome, kegg, biocarta))
  write.csv(GSEA_combined, paste0(dir, 'GSEA/', cell_type, ".csv"))
  
}

#Plotting ----
GSEA <- list.files(paste0(dir,'GSEA/'), recursive = T)



mat1<- read.csv(paste0(dir, "GSEA/", GSEA[9]), row.names = 1)
mat1 <- mat1[c(4,6,9,11, 14, 17, 420),]
mat1$cell_type <- rep("HSPC", length(mat1$pathway))


mat2 <- read.csv(paste0(dir, "GSEA/", GSEA[5]), row.names = 1)
mat2 <- mat2[mat2$pathway%in% mat1$pathway,]
mat2$cell_type <- rep("Monocytes", length(mat2$pathway))


mat3 <- read.csv(paste0(dir, "GSEA/", GSEA[6]), row.names = 1)
mat3 <- mat3[mat3$pathway%in% mat1$pathway,]
mat3$cell_type <- rep("Dendritic Cells", length(mat3$pathway))

mat <- rbind(mat1, mat2, mat3)

mat$cell_type <- factor(mat$cell_type, levels = c("HSPC", "Monocytes", "Dendritic Cells"))
mat$pathway <- factor(mat$pathway, levels=c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                            "HALLMARK_INFLAMMATORY_RESPONSE"  , "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                            "HALLMARK_FATTY_ACID_METABOLISM", "REACTOME_VESICLE_MEDIATED_TRANSPORT", 
                                            "HALLMARK_PROTEIN_SECRETION"))
dotplot_DE <- ggplot() + 
  geom_point(mat, mapping= aes(y=pathway, x=cell_type, size = padj, fill=NES), color='black', shape=21) + 
  scale_fill_gradient2(low = '#5788c9', mid = 'white', high ='#db3e25', midpoint = 0)+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(5, 3)) +
  #facet_wrap(~Comparison) +  
  guides(fill=guide_colorbar(title='NES'), size=guide_legend(title='Adjusted\nP-value', reverse = T)) + 
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position='right',
        axis.text.y=element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = 8, angle=45, vjust=1, hjust=1, color='black'),
        legend.key.size = unit(0.5, "cm"),
        legend.title=element_text(size=9, face='bold'), 
        legend.text=element_text(size=8)) +coord_flip()
dotplot_DE

ggsave(paste0(dir,"gsea_dotplot.svg"), units = "in", height = 4, width = 8)


#SSGSEA----
#adapted from https://rpubs.com/pranali018/SSGSEA
seu <- subset(seu, subset = cell_annos == "HSPC")
pseudo_seu <- AggregateExpression(seu, assays = "RNA", return.seurat = T, group.by = c("Disease", "Specimen_ID"))
pseudo_seu@meta.data

MDSPrimID <- read.csv(paste0(dir,'results/Pawan/LASS0_genes.csv'))[,2]
gene_sets <- as.list(as.data.frame(MDSPrimID))

res <- ssgsea(pseudo_seu@assays$RNA$data, gene_sets, scale = TRUE, norm = FALSE)
res1 <- as.data.frame(t(res)) 

res1$sample <- rownames(res1)
apply(res1$sample, )
res1$Disease <- unlist(lapply(res1$sample,FUN = function(x){gsub(x = x, "_.*","")}))


library(ggpubr)
ggplot(res1, aes( x = sample, y = MDSPrimID, fill = Disease)) +
  geom_bar(stat = "identity") + 
  labs(title = "ssGSEA on MDS scRNA-seq samples") +
  ylab("MDSPrimID ssGSEA score")+
  scale_fill_manual(values = c("grey","firebrick"))+
  xlab(NULL)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 9))

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}


