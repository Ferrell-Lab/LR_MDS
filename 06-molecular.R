
#Cohort molecular genetics plots ----
library(UpSetR)
library(data.table)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
dir <- dir
mutations <- data.table(read.csv(paste0(dir, "/doc/metadata/upset_v2.csv")))

upset(mutations, nsets = 17, number.angles = 30, point.size = 3.5, line.size = 2,
      mainbar.y.label = "Mutation Intersections", sets.x.label = "Samples Per Mutation",
      text.scale = c(1.3, 1.3, 1, 1, 2, 1.5))




mutation_type = list(
  'DNA Methylation' = c("TET2","DNMT3A"),
  'Signaling' = c("JAK2"),
  'Chromatin Modification' = c("ASXL1"),
  'RNA Splicing' = c("SF3B1","SRSF2"),
  'Tumor Suppressors' = c("TP53","CUX1")
)


mutations <- reshape2::melt(mutations)

mutations <- mutations %>%
  mutate(across(variable, ~ case_when
                (variable == "TET2" ~ 'DNA Methylation',
                  variable == "DNMT3A" ~ 'DNA Methylation',
                  variable == "JAK2" ~ 'Signaling',
                  variable == "ASXL1" ~ 'Chromatin Modification',
                  variable == "SF3B1" ~ 'RNA Splicing',
                  variable == "SRSF2" ~ 'RNA Splicing',
                  variable == "TP53" ~ 'Tumor Suppressors',
                  variable == "CUX1" ~ 'Tumor Suppressors')))
mutation_class <- pivot_wider(mutations, names_from = variable, values_fn = sum) %>% t()
mutation_class



library(scales)
library(ggthemes)
colnames(mutations) <- c("Samples", "Category", "value")

s1 <- ggplot(mutations, aes(x = Category, y = value, fill = Category)) + geom_col() + 
  xlab("Mutation Class")+
  ylab("# of Samples")+
  scale_fill_tableau(palette = 'Classic Cyclic')+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.text.y = element_blank(),
        axis.text = element_text(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
s1


#Clinical data table
library(table1)
library(readxl)
library(ggpubr)
dir <- dir
clinical <- read_xlsx(paste0(dir, "doc/metadata/clinical_data_table.xlsx"))

label(clinical$Myeloid) <- "Myeloid Dysplasia"
label(clinical$Erythroid) <- "Erythroid Dysplasia"
label(clinical$Megakaryocytic) <- "Megakaryocytic Dysplasia"
clinical$`Multilineage_Dysplasia` <- factor(clinical$`Multilineage_Dysplasia`, levels = c(1:3), labels = c("One","Two","Three" ))
label(clinical$`Multilineage_Dysplasia`) <- "Multilineage Dysplasia"
label(clinical$`IPSS_R`) <- "IPSS-R"
label(clinical$`IPSS_M`) <- "IPSS-M"
units(clinical$Age)       <- "years"

footnote <- "ᵃ At Biopsy from Bone Marrow Aspirate"
pdf(paste0(dir,"MDS_clin_table.pdf"), width = 8.5, height = 11)
p1 <- table1(~ Sex + Age + Myeloid + Erythroid + Megakaryocytic + Multilineage_Dysplasia + IPSS_R + IPSS_M,
       footnote = footnote, overall=c(left="Total MDS"), data=clinical)
ggsave(p1,filename= paste0(dir,"MDS_clin_table.pdf"), units = "in", height = 11, width = 8.5) 

dev.off()



