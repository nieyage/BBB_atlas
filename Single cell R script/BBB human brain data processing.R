###BBB human brain data processing 
#DencontX #
## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
AD_Ctx_1<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Ctx_1"),project="AD_Ctx_1",assay ="RNA")
AD_Ctx_2<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Ctx_2"),project="AD_Ctx_2",assay ="RNA")
AD_Ctx_3<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Ctx_3"),project="AD_Ctx_3",assay ="RNA")
AD_Ctx_4<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Ctx_4"),project="AD_Ctx_4",assay ="RNA")
AD_Hpc_1<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_1"),project="AD_Hpc_1",assay ="RNA")
AD_Hpc_2<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_2"),project="AD_Hpc_2",assay ="RNA")
AD_Hpc_3<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_3"),project="AD_Hpc_3",assay ="RNA")
AD_Hpc_4<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_4"),project="AD_Hpc_4",assay ="RNA")
AD_Hpc_5<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_5"),project="AD_Hpc_5",assay ="RNA")
AD_Hpc_6<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_6"),project="AD_Hpc_6",assay ="RNA")
AD_Hpc_7<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_7"),project="AD_Hpc_7",assay ="RNA")
AD_Hpc_8<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_8"),project="AD_Hpc_8",assay ="RNA")
AD_Hpc_9<-  CreateSeuratObject(counts = Read10X(data.dir = "./AD_Hpc_9"),project="AD_Hpc_9",assay ="RNA")
NC_Ctx_1<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Ctx_1"),project="NC_Ctx_1",assay ="RNA")
NC_Ctx_2<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Ctx_2"),project="NC_Ctx_2",assay ="RNA")
NC_Ctx_3<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Ctx_3"),project="NC_Ctx_3",assay ="RNA")
NC_Ctx_4<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Ctx_4"),project="NC_Ctx_4",assay ="RNA")
NC_Hpc_1<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_1"),project="NC_Hpc_1",assay ="RNA")
NC_Hpc_2<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_2"),project="NC_Hpc_2",assay ="RNA")
NC_Hpc_3<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_3"),project="NC_Hpc_3",assay ="RNA")
NC_Hpc_4<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_4"),project="NC_Hpc_4",assay ="RNA")
NC_Hpc_5<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_5"),project="NC_Hpc_5",assay ="RNA")
NC_Hpc_6<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_6"),project="NC_Hpc_6",assay ="RNA")
NC_Hpc_7<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_7"),project="NC_Hpc_7",assay ="RNA")
NC_Hpc_8<-  CreateSeuratObject(counts = Read10X(data.dir = "./NC_Hpc_8"),project="NC_Hpc_8",assay ="RNA")

objList <- list(AD_Ctx_1,AD_Ctx_2,AD_Ctx_3,AD_Ctx_4,AD_Hpc_1,AD_Hpc_2,AD_Hpc_3,AD_Hpc_4,AD_Hpc_5,AD_Hpc_6,AD_Hpc_7,AD_Hpc_8,
AD_Hpc_9,NC_Ctx_1,NC_Ctx_2,NC_Ctx_3,NC_Ctx_4,NC_Hpc_1,NC_Hpc_2,NC_Hpc_3,NC_Hpc_4,NC_Hpc_5,NC_Hpc_6,NC_Hpc_7,NC_Hpc_8)

for (i in seq_len(length(objList))) {
		 objList[[i]][["percent.mt"]] <- PercentageFeatureSet(objList[[i]], pattern = "^MT-")
      }

pdf("./results/QC_plot.pdf")

for (i in seq_len(length(objList))){
	p1<-VlnPlot(objList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	p1<-p1&ggtitle (objList[[i]]@project.name)
	print(p1);
}

dev.off()

# Retain cells that express not more than 4000 genes (remove potential homotypic doublets ?) and had mitochondrial content <15%
for (i in seq_len(length(objList))){
	objList[[i]] = subset(objList[[i]],nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt < 5 )
}

merged_obj <- merge(x=objList[[1]],y=c(objList[c(2:25)]),add.cell.ids = c("AD_Ctx_1","AD_Ctx_2","AD_Ctx_3","AD_Ctx_4",
	"AD_Hpc_1","AD_Hpc_2","AD_Hpc_3","AD_Hpc_4","AD_Hpc_5","AD_Hpc_6","AD_Hpc_7","AD_Hpc_8","AD_Hpc_9",
"NC_Ctx_1","NC_Ctx_2","NC_Ctx_3","NC_Ctx_4",
"NC_Hpc_1","NC_Hpc_2","NC_Hpc_3","NC_Hpc_4","NC_Hpc_5","NC_Hpc_6","NC_Hpc_7","NC_Hpc_8"),project = "BBB")
Idents(merged_obj) <- gsub("_................-1", "", colnames(merged_obj))

#split the combined object into a list, with each dataset as an element
BBB_EC.list <- SplitObject(merged_obj,split.by = "ident")
BBB_EC.list <- lapply(X = BBB_EC.list, FUN = function(x) {
    x <- NormalizeData(x,verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})



features <- SelectIntegrationFeatures(object.list = BBB_EC.list)
BBB_EC.list <- lapply(X = BBB_EC.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = BBB_EC.list, reference = c(1,5,14,18), reduction = "rpca",dims = 1:50)
BBB_EC.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
##switch to integrated assay
DefaultAssay(BBB_EC.integrated) <- "integrated"

BBB_EC.integrated <- ScaleData(BBB_EC.integrated, verbose = FALSE)
BBB_EC.integrated <-    RunPCA(BBB_EC.integrated, verbose = FALSE)
# Determine the ‘dimensionality’ of the dataset
# JackStraw 
BBB_EC.integrated <- JackStraw(BBB_EC.integrated, num.replicate = 100, dims =50)
BBB_EC.integrated  <- ScoreJackStraw(BBB_EC.integrated, dims = 1:50)
pdf("/md01/nieyg/project/BBB/data/Nature_GSE163577_RAW/results/seclect_pca_allcelltype2.pdf")
JackStrawPlot(BBB_EC.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(BBB_EC.integrated,ndims=50)
dev.off()

saveRDS(BBB_EC.integrated, file = "Human_brain_allcelltype-changecutoff.rds")
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
set.seed(1234)
BBB_EC.integrated<-readRDS("Human_brain_allcelltype-changecutoff.rds")
DefaultAssay(BBB_EC.integrated)<-"integrated";
BBB_EC.integrated <- FindNeighbors(object = BBB_EC.integrated, dims = 1:30)
BBB_EC.integrated <- FindClusters(object = BBB_EC.integrated, resolution = 0.2)
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:30)
table(BBB_EC.integrated@meta.data$orig.ident,BBB_EC.integrated$seurat_clusters)

saveRDS(BBB_EC.integrated, file = "Human_brain_allcelltype-changecutoff.rds")


# #downsample 
# down<-sample(colnames(BBB_EC.integrated),100000)
# BBB_EC.integrated<-subset(BBB_EC.integrated,cell=down)
# 
# 
# DefaultAssay(BBB_EC.integrated)<-"RNA";
# BBB_EC.integrated <- NormalizeData(BBB_EC.integrated)


DefaultAssay(BBB_EC.integrated)<-"integrated";
BBB_EC.integrated <- FindNeighbors(object = BBB_EC.integrated, dims = 1:30)
BBB_EC.integrated <- FindClusters(object = BBB_EC.integrated, resolution = 0.5)
BBB_EC.integrated <- RunUMAP(object = BBB_EC.integrated, dims = 1:30)
table(BBB_EC.integrated$seurat_clusters)

pdf("./results/Allcelltype_cluster_UMAP.pdf")
DimPlot(object = BBB_EC.integrated, group.by = "orig.ident",label = F)
DimPlot(object = BBB_EC.integrated, reduction = "umap",label = T)
dev.off()

#cluster annotation
DefaultAssay(BBB_EC.integrated)<-"RNA";
pdf("./results/cluster_annotation-allcelltype.pdf",width=10,height=10)

features <- c( "FLT1","CLDN5",#Endothelial
	"VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1",#Tip-like
	"LAMA2","ATP1A2",#Pericyte
	"SLIT3","RCAN2",#SMC
	"ABCA9","FBLN1",#Peri.fibroblast
	"SLC47A1",#Meningeal fibroblast
	"GFAP",#Astrocyte
	"ST18",#Oligodendrocyte
	"DSCAM",#OPC
	"DOCK8",#Microglia
	"RBFOX1",#Neuron
	"CFAP299",#Ependymal
	"PTPRC","CD247"#T cell
	)
DotPlot(BBB_EC.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

saveRDS(BBB_EC.integrated, file = "Human_brain_allcelltype-first.rds")




#  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#  library(DoubletFinder)
#  sweep.res.list_kidney <- paramSweep_v3(BBB_EC.integrated, PCs = 1:30, sct = T)
#  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
#  bcmvn_kidney <- find.pK(sweep.stats_kidney)
#  pK_BBB<-bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)] %>% as.character() %>% as.numeric()
#  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#  homotypic.prop <- modelHomotypic(BBB_EC.integrated$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
#  DoubletRate=ncol(BBB_EC.integrated)*8*1e-6
#  nExp_poi <- round(DoubletRate*ncol(BBB_EC.integrated))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#  
#  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#  #seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#  BBB_EC.integrated_DoubletFinder <- doubletFinder_v3(BBB_EC.integrated, PCs = 1:30, 
#  	pN = 0.25, pK = pK_BBB, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#  
#  pdf("/md01/nieyg/project/BBB/data/Nature_GSE163577_RAW/results/Allcelltype_cluster_DoubletFinder.pdf")
#  
#  DimPlot(BBB_EC.integrated_DoubletFinder,reduction="umap",group.by="DF.classifications_0.25_0.03_69063")
#  DimPlot(BBB_EC.integrated_DoubletFinder,reduction="umap",group.by="seurat_clusters",label = F)
#  
#  dev.off()


#filter_BBB<-subset(BBB_EC.integrated_DoubletFinder,DF.classifications_0.25_0.03_69063=="Singlet")
#DefaultAssay(filter_BBB)<-"integrated";
#filter_BBB <- FindNeighbors(object = filter_BBB, dims = 1:30)
#filter_BBB <- FindClusters(object = filter_BBB, resolution = 0.2)
#filter_BBB <- RunUMAP(object = filter_BBB, dims = 1:30)

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
set.seed(1234)
BBB_EC.integrated<-readRDS("Human_brain_allcelltype-first.rds")

filter_BBB<-subset(BBB_EC.integrated,seurat_clusters%in%as.character(c(0:11,13,15:17,20:22)))
DefaultAssay(filter_BBB)<-"integrated";
filter_BBB <-    RunPCA(filter_BBB, verbose = FALSE)
filter_BBB <- FindNeighbors(object = filter_BBB, dims = 1:30)
filter_BBB <- FindClusters(object = filter_BBB, resolution = 0.5)
filter_BBB <- RunUMAP(object = filter_BBB, dims = 1:30)

pdf("./results/Allcelltype_cluster_UMAP-filter1.pdf")
DimPlot(object = filter_BBB, group.by = "orig.ident",label = F)
DimPlot(object = filter_BBB, reduction = "umap",label = T)
dev.off()
#saveRDS(BBB_EC.integrated, file = "Human_brain_allcelltype-changecutoff2.rds")

#cluster annotation
DefaultAssay(filter_BBB)<-"RNA";
pdf("./results/cluster_annotation-allcelltype-filter1.pdf",width=10,height=8)

features <- c( "FLT1","CLDN5",#Endothelial
	"VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1",#Tip-like
	"LAMA2","ATP1A2",#Pericyte
	"SLIT3","RCAN2",#SMC
	"ABCA9","FBLN1",#Peri.fibroblast
	"SLC47A1",#Meningeal fibroblast
	"GFAP",#Astrocyte
	"ST18",#Oligodendrocyte
	"DSCAM",#OPC
	"DOCK8",#Microglia
	"RBFOX1",#Neuron
	"CFAP299",#Ependymal
	"PTPRC","CD247"#T cell
	)
DotPlot(filter_BBB, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()


filter_BBB<-subset(filter_BBB,seurat_clusters%in%as.character(c(0:7,9:11,13:19)))
DefaultAssay(filter_BBB)<-"integrated";
filter_BBB <-    RunPCA(filter_BBB, verbose = FALSE)
filter_BBB <- FindNeighbors(object = filter_BBB, dims = 1:30)
filter_BBB <- FindClusters(object = filter_BBB, resolution = 0.5)
filter_BBB <- RunUMAP(object = filter_BBB, dims = 1:30)


pdf("./results/Allcelltype_cluster_UMAP-filter2.pdf")
DimPlot(object = filter_BBB, group.by = "orig.ident",label = F)
DimPlot(object = filter_BBB, reduction = "umap",label = T)
dev.off()
#saveRDS(BBB_EC.integrated, file = "Human_brain_allcelltype-changecutoff2.rds")

#cluster annotation
DefaultAssay(filter_BBB)<-"RNA";
pdf("./results/cluster_annotation-allcelltype-filter2.pdf",width=10,height=8)

features <- c( "FLT1","CLDN5",#Endothelial
	"VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1",#Tip-like
	"LAMA2","ATP1A2",#Pericyte
	"SLIT3","RCAN2",#SMC
	"ABCA9","FBLN1",#Peri.fibroblast
	"SLC47A1",#Meningeal fibroblast
	"GFAP",#Astrocyte
	"ST18",#Oligodendrocyte
	"DSCAM",#OPC
	"DOCK8",#Microglia
	"RBFOX1",#Neuron
	"CFAP299",#Ependymal
	"PTPRC","CD247"#T cell
	)
DotPlot(filter_BBB, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()


#####add celltype information in umap#######
Idents(filter_BBB)<-filter_BBB$seurat_clusters
new.cluster.ids <- c("Oligo.","EC","Pericyte","Astrocyte","Astrocyte",
	"EC","EC","SMC","Microglia","Oligo.","Oligo.","Peri.fibroblast",
	"OPC","Neuron","Neuron","Ependymal","T cell","Pericyte")

names(new.cluster.ids) <- levels(filter_BBB)
filter_BBB <- RenameIdents(filter_BBB, new.cluster.ids)


saveRDS(filter_BBB, file = "BBB_filter_allcelltype-last.rds")


library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
set.seed(1234)
filter_BBB<-readRDS("BBB_filter_allcelltype-last.rds")
DefaultAssay(filter_BBB)<-"integrated";
EC_cells <- subset(filter_BBB,seurat_clusters%in%as.character(c(1,5,7)))
EC_cells.integrated <- RunPCA(EC_cells, npcs = 50, verbose = FALSE)

DefaultAssay(EC_cells.integrated)<-"integrated";
EC_cells.integrated <- FindNeighbors(object = EC_cells.integrated, dims = 1:20)
EC_cells.integrated <- FindClusters(object = EC_cells.integrated, resolution = 0.1)
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:20)

pdf("./results/EC_cluster_Umap.pdf")
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, reduction = "umap",label = T)
dev.off()

Clusters
were defined by established zonation markers, such as arterial VEGFC
and ALPL; capillary MFSD2A and SLC7A5; and venous IL1R1 and NR2F2

 ‘tip’ cells (for example,
PLAUR and LAMB1) as well as ‘proteostatic’ heat shock proteins.

#EC cluster annotation
DefaultAssay(EC_cells.integrated)<-"RNA";
pdf("./results/cluster_annotation-EC.pdf",width=8,height=8)

features <- c( "VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1"#Tip-like
	)
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()




EC_cells.integrated2 <- subset(EC_cells.integrated,seurat_clusters%in%as.character(c(0,1,2,6)))

DefaultAssay(EC_cells.integrated2)<-"integrated";
EC_cells.integrated2 <- FindNeighbors(object = EC_cells.integrated2, dims = 1:20)
EC_cells.integrated2 <- FindClusters(object = EC_cells.integrated2, resolution = 0.2)
EC_cells.integrated2 <- RunUMAP(object = EC_cells.integrated2, dims = 1:20)

pdf("./results/EC_cluster_Umap2.pdf")
DimPlot(object = EC_cells.integrated2, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated2, reduction = "umap",label = T)
dev.off()

DefaultAssay(EC_cells.integrated2)<-"RNA";
pdf("./results/cluster_annotation-EC2.pdf",width=8,height=8)

features <- c( "FLT1","CLDN5",#Endothelial
    "VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1"#Tip-like
	)
DotPlot(EC_cells.integrated2, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()


EC_cells.integrated3 <- subset(EC_cells.integrated2,seurat_clusters%in%as.character(c(0,1,2,6)))

DefaultAssay(EC_cells.integrated3)<-"integrated";
EC_cells.integrated3 <- FindNeighbors(object = EC_cells.integrated3, dims = 1:20)
EC_cells.integrated3 <- FindClusters(object = EC_cells.integrated3, resolution = 0.1)
EC_cells.integrated3 <- RunUMAP(object = EC_cells.integrated3, dims = 1:20)

pdf("./results/EC_cluster_Umap3.pdf")
DimPlot(object = EC_cells.integrated3, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated3, reduction = "umap",label = T)
dev.off()

DefaultAssay(EC_cells.integrated3)<-"RNA";
pdf("./results/cluster_annotation-EC3.pdf",width=8,height=8)

features <- c( "FLT1","CLDN5",#Endothelial
    "VEGFC","ALPL",#A
	"MFSD2A","SLC7A5",#C
	"IL1R1","NR2F2",#V
	"PLAUR","LAMB1"#Tip-like
	)
DotPlot(EC_cells.integrated3, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()



Idents(EC_cells.integrated3)<-EC_cells.integrated3$seurat_clusters
new.cluster.ids <- c("Capillary","Venous","Arterial","Tip-like")

names(new.cluster.ids) <- levels(EC_cells.integrated3)
EC_cells.integrated3 <- RenameIdents(EC_cells.integrated3, new.cluster.ids)

pdf("./results/ECcell_add_celltype_umap.pdf")
DimPlot(EC_cells.integrated3, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "EC_cells")
dev.off()

sampleinfo<-(EC_cells.integrated3$orig.ident)
group<-gsub("_.*","",sampleinfo)
str(group)
EC_cells.integrated3$group<-group
saveRDS(EC_cells.integrated3, file = "EC_cells.integrated-last.rds")

#Figure S7 
allcelltype<-filter_BBB
allcelltype<-readRDS("BBB_filter_allcelltype-last.rds")


pdf("./results/allcelltype_add_celltype_umap.pdf")
DimPlot(allcelltype, reduction = "umap", raster=FALSE, label = TRUE, pt.size = 0.5) 
dev.off()

pdf("./results/allcelltype_CLDN5_umap.pdf")
FeaturePlot(allcelltype, reduction = 'umap',raster=FALSE, features = "CLDN5", ncol = 1)
dev.off()

EC_cells.integrated3$subtype<-Idents(EC_cells.integrated3)
EC_cells.integrated3$subtype<-factor(EC_cells.integrated3$subtype,levels=c("Arterial","Capillary","Venous","Tip-like"))
Idents(EC_cells.integrated3)<-EC_cells.integrated3$subtype

pdf("./results/EC_3gene_Vln.pdf")
VlnPlot(EC_cells.integrated3, features = c("CLDN5","ALPL","IL1R1"), ncol = 1, pt.size = 0) 
dev.off();

pdf("./results/ECcell_add_celltype_umap.pdf",width=6,height=5)
DimPlot(EC_cells.integrated3, reduction = "umap", label = TRUE, pt.size = 0.5) 
dev.off()

R221<-read.table("./R221.txt",header=T)
R221<- R221$ genename.221
# R221_human<-toupper(R221)
which(R221_human%in% rownames(EC_cells.integrated3))

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
head(mouse_human_genes)
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}
human_R221<-convert_mouse_to_human(R221)
which(human_R221%in% rownames(EC_cells.integrated3))

# select the normal control group 
EC_normal<-subset(EC_cells.integrated3,group%in%"NC")
EC_normal<-ScaleData(EC_normal);

p<-DotPlot(EC_normal, features = human_R221[which(human_R221%in% rownames(EC_normal))]) +  
xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data_Capillary<-dotplot_data[which(dotplot_data$id%in%"Capillary"),]
write.csv(dotplot_data_Capillary,"Normal_control_Capillary_R221gene.csv")

###
data<-dotplot_data_Capillary
data$cutoff = ifelse(data$pct.exp > 0.1, "pct.exp > 0.1","pct.exp < 0.1")
pdf("./results/Normal_control_Capillary_R221gene.pdf")
ggplot(data=data, aes(x=avg.exp, y=pct.exp, color=cutoff)) +
  geom_point(alpha=0.8,size=1) +  #以点图形式呈现
  geom_hline(yintercept = 0.1,lty=4,col="black",lwd=0.8) +
  theme_bw(base_size=15)+   #去除灰色背景并设置字体大小
  scale_colour_manual(values = c('grey',"blue")) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_text(size=15,hjust = 0.5))

dev.off()


#### Figure 9 
EC_cells.integrated3

sampleinfo<-(EC_cells.integrated3$orig.ident)
sample_info<-gsub("_.$","",sampleinfo)
treat<-gsub("^.._","",sample_info)
EC_cells.integrated3$condition<-treat

EC_cells_Ctx<-subset(EC_cells.integrated3,condition%in%"Ctx")
Idents(EC_cells_Ctx) <- EC_cells_Ctx$group
EC.markers <- FindMarkers(EC_cells_Ctx,ident.1 = "AD",ident.2 = "NC", min.pct = 0.1, logfc.threshold = 0)
	mydata<-EC.markers
	mydata$Condition=ifelse(mydata$avg_log2FC>=0.5 & mydata$p_val_adj<=0.05,"up",ifelse(mydata$avg_log2FC<=-0.5 & mydata$p_val_adj<=0.05,"down","normal"))
write.csv(mydata,"./results/Ctx_Capillary_ADvsNC.csv")

pdf("./results/Ctx_Capillary_ADvsNC.pdf")
ggplot(data=mydata, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+
              ggtitle("Ctx AD vs NC")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('up'='red','down'='blue','normal'='gray'));

dev.off()

EC_cells_Hpc<-subset(EC_cells.integrated3,condition%in%"Hpc")
Idents(EC_cells_Hpc) <- EC_cells_Hpc$group;
EC.markers <- FindMarkers(EC_cells_Hpc,ident.1 = "AD",ident.2 = "NC", min.pct = 0.1, logfc.threshold = 0)
mydata<-EC.markers
mydata$Condition=ifelse(mydata$avg_log2FC>=0.5 & mydata$p_val_adj<=0.05,"up",ifelse(mydata$avg_log2FC<=-0.5 & mydata$p_val_adj<=0.05,"down","normal"))

pdf("./results/Hpc_Capillary_ADvsNC.pdf")
ggplot(data=mydata, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=Condition)) + 
              geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+
              ggtitle("Ctx AD vs NC")+theme_bw()+theme(panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
              geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
              scale_color_manual(values=c('up'='red','down'='blue','normal'='gray'));
dev.off()

write.csv(mydata,"./results/Hpc_Capillary_ADvsNC.csv")



