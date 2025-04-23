library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)
library(monocle3)
library(magrittr)
library(ggplot2)
library(scales)
library(viridis)
library(gridExtra)

set.seed(100)


dir <- dir
seu <- readRDS(paste0(dir,'seurat_combined.rds'))
LASS0 <- read.csv(paste0(dir,'LASS0_genes.csv'))[,2]


umap_df <- cbind(seu@reductions$umap@cell.embeddings, seu@meta.data)

#Subset HSPC/mature myeloid metaclusters
Idents(seu) <- "cell_metaclusters"
seu <- subset(seu, idents = c("HSPC", "Myeloid"))


#Create CDS object for Monocle3 analysis
counts <- seu@assays$RNA$counts
cell_metadata <- seu@meta.data
cds <- new_cell_data_set(counts, cell_metadata = cell_metadata,
                         gene_metadata = data.frame(gene_short_name = rownames(counts),
                                                    row.names = rownames(counts)))

reducedDim(cds, "UMAP", withDimnames=TRUE) <- seu[['umap']]@cell.embeddings
umap_cds <- as.data.frame(cbind(reducedDim(cds, "UMAP", withDimnames=TRUE), colData(cds)))
p1 <- ggplot()+geom_point(data = umap_df, aes(umap_1,umap_2, color = "blank"),  size = 0.5)+
  geom_point(data = umap_cds, aes(umap_1,umap_2, color = cell_metaclusters), size = 0.5)+
  scale_color_manual(name = "Cell Populations",
                     values = c("#A0CBE8","#4E79A7", "grey"),
                     breaks = c("Myeloid", "HSPC", "blank"),
                     labels = c("Myeloid", "HSPC", ""))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 2)))
p1$layers[[1]]$aes_params$alpha <- 0.1                 
p1

#cluster cells
cds <- cluster_cells(cds, reduction_method = "UMAP")

#Construct the graph
cds <- learn_graph(cds, use_partition = F, close_loop = T, learn_graph_control = list("minimal_branch_len" = 20))

#create starting node/cell
cell_ids <- colnames(cds)[cell_metadata$cell_metaclusters ==  "HSPC"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

#compute the trajectory
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE, show_trajectory_graph = T)


traj.coord<- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

umap_cds <- as.data.frame(cbind(reducedDim(cds, "UMAP", withDimnames=TRUE), traj.coord))

p2 <- ggplot()+
  geom_point(data = umap_df, aes(umap_1,umap_2), size = 0.5, color = "grey")+
  geom_point(data = umap_cds, aes(umap_1,umap_2, color = traj.coord), size = 0.5)+
  #scale_color_viridis(name = "Pseudotime")+
  scale_color_viridis(name = "Pseudotime", 
                      limits = c(0,16),
                      breaks = c(0,2,8,14,16),
                      labels = c("","Early","Mid","Late",""))+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()
  )
p2$layers[[1]]$aes_params$alpha <- 0.1                 
p2

p1 / p2

ggsave(filename = paste0(dir,"pseudotime_umaps.png"), dpi = 600, units = "in", height = 6, width = 4)

# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]
endpoints <- endpoints[3]


#Calculate cell weights
cellWeights <- data.frame(endpoints = rep(1, ncol(cds)))
rownames(cellWeights) <- colnames(cds)
pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
                     nrow = ncol(cds), byrow = FALSE)


snow <- BiocParallel::SnowParam(workers = 2, type = "SOCK")

#tradeSeq condition test comparing MDS and HD

sce <- fitGAM(counts = counts,
               pseudotime = pseudotime,
               genes = LASSO,
               parallel = T, BPPARAM = snow,
               condition = factor(colData(cds)$Condition),
               cellWeights = cellWeights)

condRes <- conditionTest(sce, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]
cat(conditionGenes, sep = '\n')


#heatmap of MDS vs. HD condition test
yhatSmooth <- predictSmooth(sce2, gene = conditionGenes, nPoints = 1000, tidy = FALSE)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_MDS <- pheatmap(yhatSmoothScaled[, 1001:2000],
                           cluster_cols = FALSE,
                           color = colorRampPalette((viridis(n = 10)))(100),
                           show_rownames = F, show_colnames = FALSE, main = "MDS", legend = FALSE,
                           silent = TRUE
)

matchingHeatmap_HD <- pheatmap(yhatSmoothScaled[heatSmooth_MDS$tree_row$order, 1:1000],
                               cluster_cols = FALSE, cluster_rows = FALSE,
                               color = colorRampPalette((viridis(n = 10)))(100),
                               show_rownames = T, show_colnames = FALSE, main = "HD",
                               legend = T, silent = TRUE)
pdf(paste0(dir,"tradeseq_heatmap.pdf"), height = 8, width = 8 )
grid.arrange(heatSmooth_MDS[[4]], matchingHeatmap_HD[[4]], ncol = 2)
dev.off()

#display genes per cluster (MDS 1 and 2, HD)
cl <- sort(cutree(heatSmooth_MDS$tree_row, k = 3))
table(cl)
MDSGenes_sub <- names(cl)[cl == "2"]
cat(MDSGenes_sub, sep = '\n')


#gene-wise expression plots across pseudotime
for(i in 1:length(conditionGenes)){
  plotSmoothers(sce2, assays(sce2)$counts, gene =  conditionGenes[i], alpha = 1, border = TRUE, plotLineages = T, curvesCols = c("#C48a47","#F03B20"))+
    scale_color_manual(values = c("#C48a47","#F03B20"),
                       name="Condition",
                       breaks=c("Lineage 1_HD", "Lineage 1_MDS"),
                       labels=c("HD", "MDS"))+
    scale_fill_manual(values = c("#C48a47","#F03B20"),
                      name="Condition",
                      breaks=c("Lineage 1_HD", "Lineage 1_MDS"),
                      labels=c("HD", "MDS"))+
    ggtitle(conditionGenes[i])+
    theme(
      axis.text = element_text(color="black"),
      axis.ticks = element_line(color = "black"),
      axis.line =  element_line(color = "black"),
      legend.text = element_text(size = 10))+
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  ggsave(paste0(dir, conditionGenes[i], ".png"), units = "in", width = 3, height = 3)
}
