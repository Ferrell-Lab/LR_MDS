
library(Seurat)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(readxl)
library(SummarizedExperiment)
library(plyr)
custom_palette <-c('#D0D1E6',
                   '#67A9CF',
                   '#253494',
                   '#FEC44F',
                   '#FE9929',
                   '#EC7014',
                   '#D9F0A3',
                   '#78C679',
                   '#006837',
                   '#BAE4BC',
                   '#D4B9DA',
                   '#DF65B0',
                   '#8856A7',
                   '#810F7C',
                   '#E7298A',
                   '#FC9272',
                   '#EF3B2C',
                   '#99000D'
)

#Read data ----
BMB_rds <- list.files(dir, pattern = "NoNegativeDrops.rds", recursive = T)
Other_rds <- list.files(dir, pattern = "CellFiltered.rds", recursive = T)
#Liu_rds <- Other_rds[c(11,13,14,15,16)]
#Healthy_rds <- Other_rds[c(1:7,12, 17:23)]
#Tom_rounds <- BMB_rds[c(4:10)]
Healthy_rds <- Other_rds[c(1:4,12, 26:28)] #took out patient C2 (sample 3), ck (sample 4), SK1 (sample 14),SK2 (sample 15), S2 (sample 13), M (sample 8) N (sample 9)
all_rds <- c(Healthy_rds, BMB_rds)

seu <- list()
seu <- lapply(paste0(dir,all_rds),readRDS)

#Annotate HD samples
seu[[1]]@meta.data$lane <- 'patient_A'
seu[[2]]@meta.data$lane <- 'patient_P'
seu[[3]]@meta.data$lane <- 'patient_Q'
seu[[4]]@meta.data$lane <- 'patient_S1'
seu[[5]]@meta.data$lane <- 'patient_C1'
seu[[6]]@meta.data$lane <- 'patient_G'
seu[[7]]@meta.data$lane <- 'patient_K'
seu[[8]]@meta.data$lane <- 'patient_L'

#Annotate MDS lanes
seu[[9]]@meta.data$lane <- 'BMBX_1'
seu[[10]]@meta.data$lane <- 'BMBX_2'
seu[[11]]@meta.data$lane <- 'BMBX_3'


seu[[12]]@meta.data$lane <- 'R1'
seu[[13]]@meta.data$lane <- 'R2'
seu[[14]]@meta.data$lane <- 'R3'
seu[[15]]@meta.data$lane <- 'R5'
seu[[16]]@meta.data$lane <- 'R6'
seu[[17]]@meta.data$lane <- 'R8'
seu[[18]]@meta.data$lane <- 'R9'





# merge objects and add sample id
seu <- merge(x = seu[[1]], y = seu[c(2:18)], add.cell.ids = c("H1","H2","H3","H4","H5","H6","H7","H8",
                                                              "X1","X2","X3",
                                                              "R1","R2","R3","R5","R6","R8","R9"))

seu@meta.data$barcodes <- rownames(seu@meta.data) #add barcodes as metadata column to preserve them


#Clean HTO labels
seu$HTO_classification.global[is.na(seu$HTO_classification.global)] <- "Singlet"

seu <- subset(seu, subset = HTO_classification.global == "Singlet")

seu@meta.data$HTO_classification <- gsub(pattern = "-", replacement = "", x =seu@meta.data$HTO_classification)

#Add metadata
seu_metadata <- read_xlsx(paste0(dir,"seu_integrated_metadata.xlsx"))
metadata <- seu@meta.data
metadata <- merge(metadata, seu_metadata, by = c("10X_Lane", "HTO_classification"), all.x = T)
rownames(metadata) <- metadata$barcodes
seu <- AddMetaData(seu, metadata)
dim(seu)

seu <- seu[,!is.na(seu@meta.data$Specimen_ID)] #removes samples from lanes not used in study

Idents(seu) <- "Disease"
seu <- subset(seu, downsample = length(seu$Disease[seu$Condition == "MDS"])) #Equal number of cells in HD and MDS


#Normalization via SCTransform
seu <- SCTransform(seu, method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = T, verbose = F)

seu <- RunPCA(seu, npcs = 50, verbose = F)


#Harmony Integration
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  normalization.method = "SCT", new.reduction = "harmony",
  verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, reduction = "harmony")

res_list <- c(0.1,0.3,0.5,0.7,1.2) #find optimal cluster resolution
seu <- FindClusters(seu, resolution = res_list)

seu <- RunUMAP(seu, dims = 1:30, reduction = "harmony")

seu$clusters <- seu$SCT_snn_res.0.3 #use 0.3 as optimal cluster resolution

#Find top 20 markers per condition
seu <- PrepSCTFindMarkers(seu)
DimPlot(seu, group.by = "SCT_snn_res.0.3",  reduction = "umap", label.box = T, label = T) + scale_color_tableau("Classic 20") 


Idents(seu) <- "cell_annos"
markers <- FindAllMarkers(seu, logfc.threshold = 0.25, cells.use = sample(colnames(seu), size = 10000))
markers %>% group_by(cluster) %>% slice_head(n = 20) %>%
  ungroup() -> top20
top20 <- as.data.frame(top20)
p1 <- DoHeatmap(seu, features = top20$gene) + NoLegend()
p1
write.table(top20, paste0(dir, 'top20markers.csv'))


#Supervised labeling of clusters
cell_annotations <- c("0" = "CD4+ T cells-1",
                      "1" = "Classical Monocytes",
                      "2" = "CD4+ T cells-2",
                      "3" = "NK cells",
                      "4" = "CD8+ T cells",
                      "5" = "B cells",
                      "6" = "HSPC",
                      "7" = "Ery-Early",
                      "8" = "Ery-Late",
                      "9" = "Ery-Mid",
                      "10" = "Ery-Mid",
                      "11" = "pro-B",
                      "12" = "Dendritic cells",
                      "13" = "Exclude",
                      "14" = "LMPP",
                      "15" = "Exclude",
                      "16" = "Plasma cells",
                      "17" = "Megakaryocytes",
                      "18" = "Exclude" #no discernible signature, likely junk
)

#map annotations to clusters
cell_map<- seu@meta.data[seu$clusters %in% names(cell_annotations),] 
seu$cell_annos <- mapvalues(seu$clusters, from = names(cell_annotations), to = cell_annotations)

seu <- subset(seu, subset = cell_annos == "Exclude", invert = T)

#remove extraneuous levels
cell_annotations <- c("0" = "CD4+ T cells-1",
                      "1" = "Classical Monocytes",
                      "2" = "CD4+ T cells-2",
                      "3" = "NK cells",
                      "4" = "CD8+ T cells",
                      "5" = "B cells",
                      "6" = "HSPC",
                      "7" = "Ery-Early",
                      "8" = "Ery-Late",
                      "9" = "Ery-Mid",
                      "10" = "Ery-Mid",
                      "11" = "pro-B",
                      "12" = "Dendritic cells",
                      "13" = "LMPP",
                      "14" = "LMPP",
                      "15" = "Plasma cells",
                      "16" = "Plasma cells",
                      "17" = "Megakaryocytes",
                      "18" = "Megakaryocytes"
)


cell_map<- seu@meta.data[seu$clusters %in% names(cell_annotations),] 
seu$cell_annos <- mapvalues(seu$clusters, from = names(cell_annotations), to = cell_annotations)

#metacluster annotation
cell_annotations <- c("0" = "T_NK",
                      "1" = "Myeloid",
                      "2" = "T_NK",
                      "3" = "T_NK",
                      "4" = "T_NK",
                      "5" = "Erythroid",
                      "6" = "HSPC",
                      "7" = "Erythroid",
                      "8" = "Erythroid",
                      "9" = "Erythroid",
                      "10" = "Erythroid",
                      "11" = "B cells",
                      "12" = "Myeloid",
                      "13" = "HSPC",
                      "14" = "HSPC",
                      "15" = "B cells",
                      "16" = "B cells",
                      "17" = "Erythroid",
                      "18" = "Erythroid")

cell_map<- seu@meta.data[seu$clusters %in% names(cell_annotations),] 
seu$cell_metaclusters <- mapvalues(seu$clusters, from = names(cell_annotations), to = cell_annotations)



saveRDS(seu, paste0(dir, 'seurat_combined.rds'))

write.csv(seu@assays$RNA$counts, file = paste0(dir, "seu_integrated_counts.csv"))



umap_df <- cbind(seu@reductions$umap@cell.embeddings, seu@meta.data)

#UMAP of all cells
p1 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = factor(cell_annos))) +
  scale_color_manual(values = custom_palette) +
  theme_minimal() +
  geom_point(size = 0.5, alpha = 0.8) +
  labs(color = "Cell Populations", caption = paste0("Number of LR-MDS cells: ",table(seu$Disease)[2]), title = "Integrated BMMCs from MDS and HD Patients") +
  #labs(color = "Cell Populations", caption = paste0(dim(seu)[2], " cells")) +
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(), 
        plot.caption = element_text(size = 9, face = "bold"),
        #axis.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 5)))
p1

ggsave(paste0(dir,"umap_integrated.png"),dpi = 300, units = "in", height = 5, width = 6)
# Area + contour
ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white", bins = 15)+
  #geom_density_2d()+
  labs( color = "")+
  scale_fill_continuous(name = "Density", 
                        type = "viridis" ,
                        limits = c(0, 0.015),
                        breaks = c(0, 0.005, 0.010, 0.015),
                        # breaks = c(0.005,0.010),
                        labels = c("","Low","High","")
  ) +
  facet_wrap(~Disease, ncol = 2)+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth= 1),
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        strip.text = element_text(
          size = 12, color = "black"
        ),
        legend.key.size =  unit(0.5, 'lines'))



ggsave(paste0(dir,"cell_densities.svg"), dpi = 300, units = "in", width = 6, height = 3)

ggplot(umap_df_proportions, aes(x = factor(cell_annos, levels = cell_annos_levels), y = proportion, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  #labs(title = "Risk Status") +
  ylab("Proportion of Cell Population")+
  xlab("Cell Population")+
  labs(legend = "Disease")+
  scale_fill_manual(values = c("#BBBBBC","#94221F")) +  # Set custom fill colors# Format y-axis labels as percentages
  #facet_grid(rows = vars(Dataset)) +  # Facet by risk status  theme(legend.position = "bottom")
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),  # Rotate x-axis labels for better readability
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA))

ggsave(paste0(dir,"cell_proportions.svg"), dpi = 300, units = "in", width = 5, height = 4)

