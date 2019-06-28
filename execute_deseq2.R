rm(list=ls())
library("DESeq2")
library("dplyr")
library("ggplot2")

options(java.parameters = "- Xmx1024m")
library(xlsx)

setwd(file.path("/Users", "gonzaleza", "rres", "Yongfang_grainProtein"))

#Counts from plate1
#counts <- read.table("counts.tab", row.names = 1, header=F, sep="\t", comment="", as.is=T)

#All samples with technical replicates
#counts <- read.table("all_sample_counts.tab", row.names = 1, header=F, sep="\t", comment="", as.is=T)
#sampleInfo <- cbind(rbind(sampleInfo, sampleInfo), technical_rep=c(rep("run2", 108), rep("run1", 108))) %>% arrange(tube_no)
#sampleInfo$group <- as.factor(paste(sampleInfo$variety, sampleInfo$N, sampleInfo$stage, sampleInfo$technical_rep, sep = "_"))

#All samples after merging technical replicates
counts <- read.table("all_samples_tech_reps_merged_counts.tab", row.names = 1, header=F, sep="\t", comment="", as.is=T)

#Read sample information from file
sampleInfo <- read.table("sample_info.txt", colClasses = c("numeric","factor", "factor", "factor"), header = T)
sampleInfo$group <- as.factor(paste(sampleInfo$variety, sampleInfo$N, sampleInfo$stage, sep = "_"))
#sampleInfo$group_noN <- as.factor(paste(sampleInfo$variety, sampleInfo$stage, sep = "_"))


  
#PCA for all the analyses
ddsMat <- DESeqDataSetFromMatrix(countData = counts, sampleInfo, ~ 0 + group)
#Only consider genes that have more than 10 counts in any of the samples
#dds <- ddsMat[rowSums(counts(ddsMat)) > 1,]
dds <- ddsMat[apply(counts(ddsMat), 1, max) > 10,]
rld <- vst(dds)
x11();plotPCA(rld, intgroup = "group")
#x11();plotPCA(rld, intgroup = "technical_rep")
x11();plotPCA(rld, intgroup = "variety")
x11();plotPCA(rld, intgroup = "N")
x11();plotPCA(rld, intgroup = "stage")

dds <- DESeq(dds)
#load("dds.rdata")#Load precomputed dds
resultsNames(dds)

##############################################################
#Identify genes DE between Skyfall and Gallant
##############################################################
#Compare the mean expression of Skyfall versus the mean expression of Gallant
#Create numeric contrast by extracting the coefficients related to Skyfall and Gallant
contrasts_Sky_v_Gal <- strsplit(resultsNames(dds), "_") %>%
                        unlist() %>%
                        grep("group", ., value = T) %>%
                        substr(.,6,nchar(.)) %>%
                        sapply(function(coeff) if(coeff=="Gallant"){-1/9}else if(coeff=="Skyfall"){1/9}else(0))
res_Sky_v_Gal <- results(dds, contrast = contrasts_Sky_v_Gal, alpha=0.05)
summary(res_Sky_v_Gal)
degs_Sky_v_Gal <- subset(res_Sky_v_Gal, padj < 0.05)
dim(degs_Sky_v_Gal)
degs_Sky_v_Gal <- degs_Sky_v_Gal[order(degs_Sky_v_Gal$padj),]
write.xlsx(degs_Sky_v_Gal,
           "variety_degs.xlsx",
           sheetName = "degs_Sky_v_Gal",
           row.names = T,
           append = F)

##############################################################
#Identify genes DE between Skyfall and Gallant in 100N
##############################################################
#Compare the mean expression of Skyfall versus the mean expression of Gallant
#Create numeric contrast by extracting the coefficients related to Skyfall and Gallant
contrasts_Sky_v_Gal_100N <- strsplit(resultsNames(dds), "N_") %>%
  unlist() %>%
  grep("group", ., value = T) %>%
  #substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="groupGallant_100"){-1/3}else if(coeff=="groupSkyfall_100"){1/3}else(0))
res_Sky_v_Gal_100N <- results(dds, contrast = contrasts_Sky_v_Gal_100N, alpha=0.05)
summary(res_Sky_v_Gal_100N)
degs_Sky_v_Gal_100N <- subset(res_Sky_v_Gal_100N, padj < 0.05)
dim(degs_Sky_v_Gal_100N)
degs_Sky_v_Gal_100N <- degs_Sky_v_Gal_100N[order(degs_Sky_v_Gal_100N$padj),]
write.xlsx(degs_Sky_v_Gal_100N,
           "variety_degs.xlsx",
           sheetName = "degs_Sky_v_Gal_100N",
           row.names = T,
           append = T)

#Compute overlap between overall DEGs and 100N DEGs
common_degs_Sky_v_Gal <- rownames(degs_Sky_v_Gal_100N) %>% intersect(rownames(degs_Sky_v_Gal))
diff_degs_Sky_v_Gal <- rownames(degs_Sky_v_Gal_100N) %>% setdiff(rownames(degs_Sky_v_Gal))

##############################################################
#Identify genes DE between Crusoe and all other varieties
##############################################################
#Compare the mean expression of Crusoe against the mean expression of Skyfall,Gallant and Solstice
#Create numeric contrast by extracting the coefficients related to Skyfall and Gallant
contrasts_Crusoe_v_allOthers <- strsplit(resultsNames(dds), "_") %>%
  unlist() %>%
  grep("group", ., value = T) %>%
  substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="Crusoe"){1/9}else if(coeff=="Gallant"){-1/27}else if(coeff=="Skyfall"){-1/27}else(-1/27))
res_Crusoe_v_allOthers <- results(dds, contrast = contrasts_Crusoe_v_allOthers, alpha=0.05)
summary(res_Crusoe_v_allOthers)
degs_Crusoe_v_allOthers <- subset(res_Crusoe_v_allOthers, padj < 0.05)
dim(degs_Crusoe_v_allOthers)
degs_Crusoe_v_allOthers <- degs_Crusoe_v_allOthers[order(degs_Crusoe_v_allOthers$padj),]
write.xlsx(degs_Crusoe_v_allOthers,
           "variety_degs.xlsx",
           sheetName = "degs_Crusoe_v_allOthers",
           row.names = T,
           append = T)

##################################################################
#Identify genes DE between Crusoe and all other varieties in 100N
##################################################################
#Compare the mean expression of Crusoe against the mean expression of Skyfall,Gallant and Solstice
#Create numeric contrast by extracting the coefficients related to Skyfall and Gallant
contrasts_Crusoe_v_allOthers_100N <- strsplit(resultsNames(dds), "N_") %>%
  unlist() %>%
  grep("group", ., value = T) %>%
  #substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="groupCrusoe_100"){1/3}else if(coeff=="groupGallant_100"){-1/9}else if(coeff=="groupSkyfall_100"){-1/9}else if(coeff=="groupSolstice_100"){-1/9}else(0))
res_Crusoe_v_allOthers_100N <- results(dds, contrast = contrasts_Crusoe_v_allOthers_100N, alpha=0.05)
summary(res_Crusoe_v_allOthers_100N)
degs_Crusoe_v_allOthers_100N <- subset(res_Crusoe_v_allOthers_100N, padj < 0.05)
dim(degs_Crusoe_v_allOthers_100N)
degs_Crusoe_v_allOthers_100N <- degs_Crusoe_v_allOthers_100N[order(degs_Crusoe_v_allOthers_100N$padj),]
write.xlsx(degs_Crusoe_v_allOthers_100N,
           "variety_degs.xlsx",
           sheetName = "degs_Crusoe_v_allOthers_100N",
           row.names = T,
           append = T)

#Compute overlap between overall DEGs and 100N DEGs
common_degs_Crusoe_v_allOthers <- rownames(degs_Crusoe_v_allOthers_100N) %>% intersect(rownames(degs_Crusoe_v_allOthers))
diff_degs_Crusoe_v_allOthers <- rownames(degs_Crusoe_v_allOthers_100N) %>% setdiff(rownames(degs_Crusoe_v_allOthers))

##############################################################
#Identify genes with N effect
##############################################################
#Use a full model with N as a factor and a reduced model without it
#dds_noN <- DESeq(dds, test = "LRT", full = model.matrix(~0+group, sampleInfo), reduced = model.matrix(~0+group_noN, sampleInfo))
#resultsNames(dds_noN)
#res_N <- results(dds_noN)
#summary(res_N)
#degs_N <- subset(res_N, padj < 0.05)
#dim(degs_N)
#degs_N <- degs_N[order(degs_N$padj),]

#100N vs 200N across all varieties and time points
contrasts_200N_v_100N <- strsplit(resultsNames(dds), "_") %>%
  unlist() %>%
  grep("N", ., value = T) %>%
  #substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="100N"){-1/12}else if(coeff=="200N"){1/12}else(0))
res_200N_v_100N <- results(dds, contrast = contrasts_200N_v_100N, alpha=0.05)
summary(res_200N_v_100N)
degs_200N_v_100N <- subset(res_200N_v_100N, padj < 0.05)
dim(degs_200N_v_100N)
degs_200N_v_100N <- degs_200N_v_100N[order(degs_200N_v_100N$padj),]
write.xlsx(degs_200N_v_100N,
           "N_degs.xlsx",
           sheetName = "degs_200N_v_100N",
           row.names = T,
           append = F)

#100N vs 350N across all varieties and time points
contrasts_350N_v_100N <- strsplit(resultsNames(dds), "_") %>%
  unlist() %>%
  grep("N", ., value = T) %>%
  #substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="100N"){-1/12}else if(coeff=="350N"){1/12}else(0))
res_350N_v_100N <- results(dds, contrast = contrasts_350N_v_100N, alpha=0.05)
summary(res_350N_v_100N)
degs_350N_v_100N <- subset(res_350N_v_100N, padj < 0.05)
dim(degs_350N_v_100N)
degs_350N_v_100N <- degs_350N_v_100N[order(degs_350N_v_100N$padj),]
write.xlsx(degs_350N_v_100N,
           "N_degs.xlsx",
           sheetName = "degs_350N_v_100N",
           row.names = T,
           append = T)

#200N vs 350N across all varieties and time points
contrasts_350N_v_200N <- strsplit(resultsNames(dds), "_") %>%
  unlist() %>%
  grep("N", ., value = T) %>%
  #substr(.,6,nchar(.)) %>%
  sapply(function(coeff) if(coeff=="200N"){-1/12}else if(coeff=="350N"){1/12}else(0))
res_350N_v_200N <- results(dds, contrast = contrasts_350N_v_200N, alpha=0.05)
summary(res_350N_v_200N)
degs_350N_v_200N <- subset(res_350N_v_200N, padj < 0.05)
dim(degs_350N_v_200N)
degs_350N_v_200N <- degs_350N_v_200N[order(degs_350N_v_200N$padj),]
write.xlsx(degs_350N_v_200N,
           "N_degs.xlsx",
           sheetName = "degs_350N_v_200N",
           row.names = T,
           append = T)

#Compute overlap of DEGS between N levels
degs_v100N <- rownames(degs_200N_v_100N) %>% intersect(rownames(degs_350N_v_100N))
#write.xlsx(cbind(degs_200N_v_100N[degs_v100N,],degs_350N_v_100N[degs_v100N,]),
#           "N_degs.xlsx",
#           sheetName = "degs_350N_and_200N_v_100N",
#           row.names = T,
#           append = T)
degs_all_N <- rownames(degs_200N_v_100N) %>% intersect(rownames(degs_350N_v_100N)) %>% intersect(rownames(degs_350N_v_200N))
write.xlsx(cbind(degs_200N_v_100N[degs_all_N,],degs_350N_v_100N[degs_all_N,], degs_350N_v_200N[degs_all_N,]),
           "N_degs.xlsx",
           sheetName = "degs_all_N",
           row.names = T,
           append = T)

gene <- rownames(degs_Sky_v_Gal_100N)[1]


#############################################################
#Study expression pattern of storage proteins
#############################################################
library("biomaRt")

#Storage proteins will be retrieved using the GO term "Nutrient reservoir activity" (GO:0045735) which,
#according to GO "can be used in place of the obsolete term 'storage protein; GO:0005187'"
plants_ensembl <- useMart("plants_mart", dataset = "taestivum_eg_gene" ,host="http://plants.ensembl.org")
storage_protein_ids <- getBM(attributes=c("ensembl_gene_id"), filters = "go", values = "GO:0045735", plants_ensembl) %>%
  pull(ensembl_gene_id)

#According to Yongfang many genes annotated with GO:0045735 are not storage proteins. New attempt done
#with InterPro domains IPR001954 and IPR036312:
#IPR001954: Gliadin/LMW glutenin
lmw_glutenin_ids <- getBM(attributes=c("ensembl_gene_id"), filters = "interpro", values = "IPR001954", plants_ensembl) %>%
  pull(ensembl_gene_id)

#IPR036312: Bifunctional inhibitor/plant lipid transfer protein/seed storage helical domain superfamily
superfamily_ids <- getBM(attributes=c("ensembl_gene_id"), filters = "interpro", values = "IPR036312", plants_ensembl) %>%
  pull(ensembl_gene_id)

#Make heatmap of the expressed storage proteins
library("pheatmap")
mat  <- assay(rld)[rownames(rld) %in% superfamily_ids,]
rownames(mat) <- NULL
colnames(mat) <- sampleInfo$group %>% as.vector() %>% make.unique(sep="_")
mat  <- mat - rowMeans(mat)
x11();clust <- pheatmap(mat, scale = "row", cutree_rows = 6)
cluster_membership <- cutree(clust$tree_row, k=6)

library(gridExtra)
x11(); grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, top="Expression pattern of genes in cluster 5")



plot_expression <- function(gene_id, dds){
  deg_expression <- cbind(sampleInfo, norm_counts=counts(dds, normalized=T)[gene_id,]+0.5)
  
  #x11()
  plot <- 
  ggplot(deg_expression, aes(x=stage, y=norm_counts, group=N, col=N)) +
    geom_point() +
    stat_summary(fun.y = "mean", geom = "line", size=1) +
    expand_limits(y=1) +
    scale_y_log10() +
    scale_colour_manual(values = c("blue", "green", "red")) +
    facet_grid(~variety) +
    ggtitle(gene_id) +
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
  #print(plot)
}
  