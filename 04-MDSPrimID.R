

library(glmnet)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
set.seed(seed = 100)

dir <- dir
seu <- readRDS('seurat_combined.rds')

#Subset to include only HSPC metacluster
seu <- subset(seu, subset = cell_metaclusters == "HSPC")

#Normalize/Scale RNA assay
DefaultAssay(seu) <- "RNA"
seu <- seu %>% NormalizeData() %>% ScaleData()
pseudo_seu <- AggregateExpression(seu, assays = "RNA", return.seurat = T, group.by = c("Disease","Specimen_ID", "Coded_ID"),  normalization.method = "LogNormalize") #Log normalization


#WGCNA (adapted from https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html) ----
library(WGCNA)
allowWGCNAThreads()

#transpose counts matrix for WGCNA
input_mat  <- t(pseudo_seu@assays$RNA$data)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,           
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


picked_power = 4
temp_cor <- cor
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME$a <- case_when(mME$name == "turquoise" ~ "turquoise",
                   mME$name == "blue" ~ "blue", .default = "black" )


mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90, color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

ggsave(file = paste0(dir,"WGCNA_modules.pdf", units = "in", width = 8, height = 8))


module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

chooseTopHubInEachModule(input_mat, colorh = mergedColors, type = "signed", power = picked_power) 

#save gene modules
write_delim(module_df,
            file = paste0(dir,"gene_modules.txt"),
            delim = "\t")


# pick out a few modules of interest here
modules_of_interest = c( "blue", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest <- pseudo_seu@assays$RNA$data[submod$gene_id,]

write.csv(row.names(expr_of_interest), paste0(dir, "WGCNA_genes.csv"))


# Calculate TOMs on single modules 
# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest[2])

expr_of_interest <- pseudo_seu@assays$RNA$data[submod$gene_id,]

# Only recalculate TOM for modules of interest
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = 4)



# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)



edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) #%>%
# left_join(module_df, join_by(gene1 == gene_id)) %>%
# dplyr::rename(module1 = colors) %>%
# left_join(module_df, join_by(gene2 == gene_id)) %>%
# dplyr::rename(module2 = colors)

write_delim(edge_list,
            file = paste0(dir,'edgelist.tsv'),
            delim = "\t")


#Trim Edge_list

edge_list <- read_delim(paste0(dir,'edgelist.tsv'))


edge_list_trim_MDS <- edge_list %>% dplyr::filter(correlation != 1) %>% 
  dplyr::filter(correlation > 0.4) %>% 
  dplyr::arrange(desc(correlation))
edge_list_trim_MDS
write.csv(unique(edge_list_trim_MDS$gene1), paste0(dir, "WGCNA_genes_MDS.csv"))



#Do the same for HD condition
submod_HD = module_df %>%
  subset(colors %in% modules_of_interest[1])

expr_of_interest <- pseudo_seu@assays$RNA$data[submod_HD$gene_id,]


TOM_HD = TOMsimilarityFromExpr(t(expr_of_interest),
                               power = 4)

# Add gene names to row and columns
row.names(TOM_HD) = row.names(expr_of_interest)
colnames(TOM_HD) = row.names(expr_of_interest)


edge_list = data.frame(TOM_HD) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) 



edge_list_trim_HD <- edge_list %>% dplyr::filter(correlation != 1) %>% 
  dplyr::filter(correlation > 0.4) %>% 
  dplyr::arrange(desc(correlation))
edge_list_trim_HD
write.csv(unique(edge_list_trim_HD)$gene1, paste0(dir, "WGCNA_genes_HD.csv"))



#Visualize hub genes in igraph
library(igraph)

edge_list_mtx <- as.matrix(edge_list_trim)[1:50,1:2] #take top 50 correlations

g <- graph_from_edgelist(edge_list_mtx)
l <- layout.fruchterman.reingold(g)

pdf(file = paste0(dir,"/igraph_WGCNA.pdf"), width = 7, height = 7)

plot(g, layout = l,
     main = "WGCNA Hub Genes",
     edge.color=rep("grey",dim(edge_list_trim)[1]),
     edge.width=edge_list_trim$correlation,
     edge.lty="solid", 
     edge.arrow.size=0, 
     vertex.color = rep("white",dim(edge_list_trim)[1]),
     vertex.frame.color = rep("black",dim(edge_list_trim)[1]),
     vertex.shape = "none",
     vertex.label.color = rep("black",dim(edge_list_trim)[1]),
     vertex.label.cex = rep(1, dim(edge_list_trim)[1]),
     vertex.label.family= "serif")

dev.off()


#Plot genes across all cells (VPS13B used as example)
seu <- readRDS('seurat_combined.rds')

library(viridis)
library(patchwork)

features <- c("VPS13B")
seu_subset <- seu[,seu$Disease == "MDS"]

p_list <- c()
for(feature in features){
  
  p <- FeaturePlot(seu_subset, features = feature, min.cutoff = 0.1)+ 
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    scale_color_viridis(option = "magma", direction = -1)+
    theme(panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_rect(fill = NA),
          axis.ticks = element_blank(),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          plot.title = element_text(size = 14),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          legend.position = "none",
    )+
    ggtitle("MDS")
  p_list[[feature]] <- p
}

seu_subset <- seu[,seu$Disease == "HD"]

r_last <- FeaturePlot(seu_subset, features = features, min.cutoff = 0.1)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_viridis(option = "magma",
                      name = substitute(paste(italic("VPS13B"), " Expression")),
                      breaks = c(0,0.25,0.65,1),
                      labels = c("","Low","","High"),
                      direction = -1)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size =  unit(1, 'lines'),
        legend.title = element_text(size = 12),
        legend.position = "bottom")+
  ggtitle("HD")
r_last

wrap_plots(p_list,  nrow = 1) / wrap_plots(r_last,   nrow = 1)
ggsave(file = paste0(dir,"VPS13B_plot.png"),dpi = 300, width = 4, height = 8)


#LASSO regularization procedure (adapted from https://github.com/firozimtech/NSCLCpred)


library(glmnet)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
set.seed(seed = 100)

dir <- dir

seu <- readRDS('seurat_combined.rds')
seu <- subset(seu, subset = cell_metaclusters == "HSPC")

DefaultAssay(seu) <- "RNA"
seu <- seu %>% NormalizeData() %>% ScaleData()

x <- seu@assays$RNA$data
seu@meta.data$Disease[seu@meta.data$Disease == "MDS"] <- 1
seu@meta.data$Disease[seu@meta.data$Disease == "HD"] <- 0


edge_list_trim_MDS <- read.csv(paste0(dir, "WGCNA_genes_MDS.csv"))[,2]
edge_list_trim_HD <- read.csv(paste0(dir, "WGCNA_genes_HD.csv"))[,2]


expr_of_interest <- c(edge_list_trim_MDS,edge_list_trim_HD)



x <- x[rownames(x) %in% expr_of_interest,]
x <- t(x)

y <- as.data.frame(seu@meta.data$Disease)


dim(y)
dim(x)


n<-nrow(x)
set.seed(100)
train_rows <- sample(1:n, .70*n)
train_rows


xtrain <- x[train_rows,]
xtest <- x[-train_rows,]

ytrain <- as.numeric(y[train_rows,])
ytest <- as.numeric(y[-train_rows,])



fit <- glmnet(xtrain, ytrain, family = "binomial", alpha = 1, standardize = T)
plot(fit)

cv.lasso <- cv.glmnet(xtrain, ytrain, alpha = 1, family = "binomial", type.measure = "auc",  standardize = TRUE) #type.measure = "auc"; type.measure = "deviance",

#plot(cv.lasso)
#AUC
summary(cv.lasso$cvm)

# minimum MSE
min(cv.lasso$cvm)

# lambda.min returns the value of lambda that gives minimum mean cross-validated error. # lambda for this min MSE
cv.lasso$lambda.min

# lambda.1se is the value of lambda that gives the most regularized model such that the cross-validated error is within one standard error of the minimum.
cv.lasso$lambda.1se

# No. of coef at Min MSE
cv.lasso$nzero[cv.lasso$lambda == cv.lasso$lambda.min]

# The best lambda value is stored inside 'cv.lasso$lambda.min'.
# identify the lambda value that produces the lowest test mean squared error (MSE).
best_lambda <- cv.lasso$lambda.min
best_lambda

df_coef_min<- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 5)

# See all contributing variables
coef_output_min<-as.matrix(df_coef_min[df_coef_min[, 1] != 0, ])


GeneName_min<-rownames(coef_output_min)
GeneName_min<-GeneName_min[-1]

# create a function to transform coefficient of glmnet and cvglmnet to data.frame
coeff2dt <- function(fitobject, s) {
  coeffs <- coef(fitobject, s)
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)
  
  # reorder the variables in term of coefficients
  return(coeffs.dt[order(coeffs.dt$coefficient, decreasing = T),])
}

coeff2dt(fitobject = cv.lasso, s = best_lambda) %>% head(20)

library(gt)
library(tidyr)
module_df <- read_delim(paste0(dir, 'gene_modules.txt'))

coeffs.table <- coeff2dt(fitobject = cv.lasso, s = best_lambda)
coeffs.table <- dplyr::left_join(coeffs.table, module_df, join_by('name' == 'gene_id'))
colnames(coeffs.table) <- c("Gene ID", "Lambda", "Module")
coeffs.table
coeffs.table.top10 <- coeffs.table %>% head(10) %>% gt() %>%  opt_table_font( weight = 0.1)%>% fmt_number(columns = c('Lambda'), decimals = 3)


coeffs.table.top10
coeffs.table<-coeffs.table[-1,]

ggplot(data = coeffs.table) +
  geom_col(aes(x = coeffs.table$"Gene ID", y = Lambda, fill = {Lambda > 0})) +
  xlab(label = "") +
  ggtitle((paste0("LASSO Coefficients with ", "lambda.min= ", round(best_lambda, 3)))) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

ggsave(paste0(dir,"lambda_coeff.pdf"), units = "in", width = 8, height = 6)


#final model
x2.train<-as.data.frame(xtrain)
x3.train<-(dplyr::select(x2.train, GeneName_min)) # select the expression data according to important Gene name
x4.train <-as.matrix(x3.train) # for cv.glment, data should be in matrix (not in data frame)
LASSO_model <- glmnet(x4.train, ytrain, alpha = 1, family = "binomial",
                                         lambda = cv.lasso$lambda.min, standardize = TRUE)

#test dataset
x2.test<-as.data.frame(xtest)
x3.test<-(dplyr::select(x2.test, GeneName_min)) # select the expression data according to important Gene name
x4.test <-as.matrix(x3.test) # for cv.glment, data should be in matrix (not in data frame)

predict_probability_selectedgene <- predict(LASSO_model,  type = "response", s=best_lambda, newx=x4.test)

library(gplots)
library(pROC)

rocobj<-roc(ytest,predict_probability_selectedgene)

auc <- round(pROC::auc(ytest, predict_probability_selectedgene),4)


ggroc(rocobj, colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) + theme_bw()

LASSO_names <- colnames(x3.train)
cat(LASSO_names, sep =  "\n")
length(LASSO_names)
LASSO_names
write.csv(LASSO_names, paste0(dir,'LASS0_genes.csv'))


#Validation ----
library(compositions)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(pROC)
library(precrec)

#validation dataset retrieved from  DOI: 10.1080/10428194.2018.1452210 
data_MDS <- read.csv(paste0(dir, 'data/MDS_bulk/GSE111085_MDS_isoform_fpkm_hgnc.txt'), header = T, row.names = 1)

metadata_MDS <- readxl::read_excel(paste0(dir, 'data/MDS_bulk/NIHMS1508337-supplement-Supp3.xlsx'))

#filter by IPSS-R Low and Int-1
metadata_MDS <- metadata_MDS %>% filter(IPSS %in% c("Low", "Int-1"))
metadata_MDS$`Dendrogram branch Fig 1` <- paste0('MDS', metadata_MDS$`Dendrogram branch Fig 1`) 
MDS_IDs <- data_MDS %>% dplyr::select(metadata_MDS$`Dendrogram branch Fig 1`)  %>% colnames()

HD_IDs <- c("MDS33", "MDS44", "MDS45", "MDS46", paste0("MDS", seq(49,67)))

data_MDS <- t(data_MDS) %>% as.data.frame() %>% mutate("ytest" = case_when(colnames(data_MDS) %in% MDS_IDs ~ 1, colnames(data_MDS) %in% HD_IDs ~ 0)) %>% na.omit()
ytest <- data_MDS$ytest
data_MDS <- t(data_MDS[,-dim(data_MDS)[2]])

#remove genes that were not sequenced in bulk dataset from LASSO gene list
excluded_genes <- LASSO_names[!LASSO_names %in% rownames(data_MDS)]


#recalculate model based on subsetted genes
GeneName_min <- GeneName_min[!GeneName_min %in% excluded_genes]

x2.train<-as.data.frame(xtrain)
x3.train<-(dplyr::select(x2.train, GeneName_min))
x4.train <-as.matrix(x3.train)


LASSO_model <- glmnet(x4.train, ytrain, alpha = 1, family = "binomial",
                                         lambda = cv.lasso$lambda.min, standardize = TRUE)


#run validation and calculate ROC curve
x2.test<-as.data.frame(t(data_MDS))
x2.test <- x2.test[,GeneName_min %in% colnames(x2.test)]
x3.test<-(dplyr::select(x2.test, GeneName_min)) 
x4.test <-as.matrix(x3.test)

predict_probability_selectedgene <- predict(LASSO_model,  type = "response", s=best_lambda, newx=x4.test)
rocobj<-roc(ytest,predict_probability_selectedgene)

auc <- round(pROC::auc(ytest, predict_probability_selectedgene),2)


ggroc(rocobj, colour = 'red', size = 2) +
  ggtitle("CD34+ LR-MDS vs. Normal Validation") + 
  labs(subtitle  = paste0('ROC Curve ', '(AUC = ', auc, ')')) + 
  theme_bw() + 
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 9),
        axis.text = element_text(color = "black"))

ggsave(file = paste0(dir,"LASSO_ROC_2.pdf"), units = "in", width = 5, height = 5)


#DE Validation ----
library(DESeq2)

data_MDS <- read.csv(paste0(dir, 'data/MDS_bulk/GSE111085_MDS_isoform_fpkm_hgnc.txt'), header = T, row.names = 1)

metadata_MDS <- readxl::read_excel(paste0(dir, 'data/MDS_bulk/NIHMS1508337-supplement-Supp3.xlsx'))

#filter by IPSS Low
metadata_MDS <- metadata_MDS %>% filter(IPSS %in% c("Low", "Int-1"))

metadata_MDS$`Dendrogram branch Fig 1` <- paste0('MDS', metadata_MDS$`Dendrogram branch Fig 1`) 

MDS_IDs <- data_MDS %>% dplyr::select(metadata_MDS$`Dendrogram branch Fig 1`)  %>% colnames()
HD_IDs <- c("MDS33", "MDS44", "MDS45", "MDS46", paste0("MDS", seq(49,67)))

data_MDS <- t(data_MDS) %>% as.data.frame() %>% mutate("ytest" = case_when(colnames(data_MDS) %in% MDS_IDs ~ 1, colnames(data_MDS) %in% HD_IDs ~ 0))%>% na.omit()
data_MDS <- t(data_MDS[,-dim(data_MDS)[2]])

coldata <- data.frame(ID = colnames(data_MDS)) %>%  mutate("condition" = case_when(colnames(data_MDS) %in% MDS_IDs ~ "MDS", colnames(data_MDS) %in% HD_IDs ~ "HD"))%>% na.omit()


dds <- DESeqDataSetFromMatrix(countData = round(data_MDS),
                              colData = coldata,
                              design = ~ condition)
dds

dds <- DESeq(dds)
res <- results(dds)

# Convert to data frame and identify significant genes
res_df <- as.data.frame(res)
res_df$significant <- case_when(res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange > 0.5 ~ "LR_MDS",res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange < -0.5 ~ "HD", res_df$padj > 0.05 ~ "No")


# Filter for significant genes that are also in your gene_labels vector
res_df_sig <- res_df[res_df$significant  %in% c("HD","LR_MDS"),]
write_delim(res_df_sig,  file = paste0(dir,"validation_DEG.txt"), delim = "\t")


label_names <- LASSO_names[LASSO_names %in% rownames(res_df_sig)]
res_df$label <- ifelse(rownames(res_df) %in% label_names, rownames(res_df), "")

# Create volcano plot with labeled significant genes
library(ggrepel)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "LR_MDS" = "red", "HD" = "darkgreen")) +
  #labs(title = "LR-MDS vs. HD CD34+ Cells Bulk RNA-seq", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA)) +
  geom_label_repel(data = res_df, aes(label = label),max.overlaps= Inf,
                   vjust = 1, nudge_x  = ifelse(res_df$significant == "LR_MDS", 1, -1), 
                   nudge_y  = ifelse(res_df$significant == "LR_MDS", 0, 1), 
                   arrow =arrow(type = "open", angle = 15, length = unit(0.1, "inches")), 
                   size = 4, color = "black")+
  geom_segment(aes(x = 2, y = 15, xend = 8, yend = 15),
               arrow = arrow(type = "closed"), color = "black", size = 1.5)+
  geom_segment(aes(x = -2, y = 15, xend = -8, yend = 15),
               arrow = arrow(type = "closed"), color = "black", size = 1.5)+
  #labs(subtitle = "CD34+ cells from Im et al., 2019")+
  annotate("text", label = "HD", x = -5, y = 16)+
  annotate("text", label = "LR-MDS", x = 5, y = 16)+
  theme(axis.text=element_text(size=11, color='black'))

ggsave(filename = paste0(dir,"mds_validation.png"), dpi = 600, height = 3, width = 4)

#SSGSEA ----
library(matrixStats)
library(data.table)
dir <- "C:/Users/bhatp3/OneDrive - Vanderbilt/2023_MDS_BM_scRNAseq/"

data_MDS <- read.csv(paste0(dir, 'data/MDS_bulk/GSE111085_MDS_isoform_fpkm_hgnc.txt'), header = T, row.names = 1)
metadata_MDS <- readxl::read_excel(paste0(dir, 'data/MDS_bulk/NIHMS1508337-supplement-Supp3.xlsx'))
LASSO_names <- read.csv(paste0(dir,'results/Pawan/LASS0_genes.csv'))[,2]
LASSO_names <- LASSO_names[!LASSO_names == "EEF1G"] #remove genes with HD coefficient

#split dataset into IPSS Low vs. high
metadata_MDS_low <- metadata_MDS %>% filter(IPSS %in% c("Low", "Int-1"))
metadata_MDS_high <- metadata_MDS %>% filter(IPSS %in% c("High", "Int-2"))

metadata_MDS_low$`Dendrogram branch Fig 1` <- paste0('MDS', metadata_MDS_low$`Dendrogram branch Fig 1`) 
metadata_MDS_high$`Dendrogram branch Fig 1` <- paste0('MDS', metadata_MDS_high$`Dendrogram branch Fig 1`)
MDS_IDs_low <- data_MDS %>% dplyr::select(metadata_MDS_low$`Dendrogram branch Fig 1`)  %>% colnames()
MDS_IDs_high <- data_MDS %>% dplyr::select(metadata_MDS_high$`Dendrogram branch Fig 1`)  %>% colnames()
HD_IDs <- c("MDS33", "MDS44", "MDS45", "MDS46", paste0("MDS", seq(49,67)))

data_MDS <- as.matrix(data_MDS)
MDSPrimID <- LASSO_names
gene_sets <- as.list(as.data.frame(MDSPrimID))

res <- ssgsea(data_MDS, gene_sets, scale = TRUE, norm = FALSE)
res1 <- as.data.frame(t(res)) 

res1$sample <- rownames(res1)
res1 <- res1 %>% mutate('IPSS-R Risk' = case_when(
                           res1$sample %in% MDS_IDs_high ~ "Higher risk", 
                           res1$sample %in% MDS_IDs_low ~ "Lower risk")) %>% na.omit()

library(ggpubr)
ggplot(res1, aes(x = res1$'IPSS-R Risk', y = MDSPrimID, color = res1$'IPSS-R Risk')) +
  geom_boxplot() + 
  geom_point(position = 'jitter')+
  labs(title = "ssGSEA on MDS bulk RNA-seq samples", caption = "Higher risk n = 18, Lower risk n = 26") +
  ylab("MDSPrimID ssGSEA score")+
  xlab('IPSS-R Risk')+
  stat_compare_means()+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text = element_text(color = 'black'),
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


