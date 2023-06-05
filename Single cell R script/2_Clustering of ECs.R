## Load required packages
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

BBB_EC.integrated<-readRDS("BBB_YMO.rds")
DefaultAssay(BBB_EC.integrated) <- "RNA" 
EC_cells <- subset(BBB_EC.integrated,idents=c("C_V_1","C_A_1","Capillary_1","C_V_2","C_A_2","Venous","Interferon","Arterial_1","Capillary_2",
  "Arterial_2","Capillary_3","Choroid_plexus"))

EC_cells <- DietSeurat(
  object=EC_cells,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
)

## split the integrated object into a list, with each dataset as an element
EC_cells.ls <- SplitObject(EC_cells,split.by = "orig.ident")
for (i in 1:length(EC_cells.ls)){
  EC_cells.ls[[i]] <- FindVariableFeatures(EC_cells.ls[[i]],selection.method = "vst", nfeatures = 3000)
}
## identify anchors using the FindIntegrationAnchors function
EC_cells.anchors <- FindIntegrationAnchors(object.list = EC_cells.ls,anchor.features = 3000,dims = 1:30)
## pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
EC_cells.integrated <- IntegrateData(anchorset = EC_cells.anchors, dims = 1:30,features.to.integrate = rownames(EC_cells.ls[[1]]))
## switch to integrated assay

DefaultAssay(EC_cells.integrated) <- "integrated"
## scale and center features in the dataset
EC_cells.integrated <- ScaleData(EC_cells.integrated, features =rownames(EC_cells.integrated))
## Perform linear dimensional reduction
EC_cells.integrated <- RunPCA(EC_cells.integrated, npcs = 50, verbose = FALSE)
## Determine the ‘dimensionality’ of the dataset
EC_cells.integrated <- JackStraw(EC_cells.integrated, num.replicate = 100, dims =50)
EC_cells.integrated <- ScoreJackStraw(EC_cells.integrated, dims = 1:50)
pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/OnlyEC_selectPC.pdf")
JackStrawPlot(EC_cells.integrated, dims = 1:50)
ElbowPlot(EC_cells.integrated,ndims=50)
dev.off()

EC_cells.integrated<-BuildClusterTree(EC_cells.integrated)
Tool(object = EC_cells.integrated, slot = 'BuildClusterTree')
PlotClusterTree(EC_cells.integrated)
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/EC_cluster_Umap.pdf")
EC_cells.integrated <- RunUMAP(object = EC_cells.integrated, dims = 1:22)
DimPlot(object = EC_cells.integrated, group.by = "orig.ident",label = T)
DimPlot(object = EC_cells.integrated, split.by = "orig.ident",label = T)
dev.off();

Umap = EC_cells.integrated@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx =EC_cells.integrated@meta.data$orig.ident)
pdf(file = "/md01/nieyg/project/BBB/YMO_results/ECplot/Umap_onlyYoung_OnlyOld.pdf", width = 8, height = 8)
p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c("Young" = "deepskyblue", "Old" = "grey"))   
p+labs(title="Young cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))

p<-ggplot(Umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
      geom_point(size = 0.2, alpha = 1) + 
      scale_color_manual(values=c( "Old"= "red", "Young" = "grey"))

p+labs(title="Aged cells distribution")+theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()


#cluster annotation

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_onlyEC_cluster_annotation.pdf",width=12,height=8)

features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
	"Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
	"Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",####1
	"Tfrc", "Car4",	"Itm2a","Chn2",#####C-V
	"Lcn2","Slc38a5","Nr2f2", "Sox12",##Venous
	"Tbx1", "Tmbs10","Icam1","Vcam1", "Vwf", "P2ry1", ####A/V
	"Plvap", "Plpp3","Esm1",####choroid plexus
	"Pdgfrb", "Cspg4","Kcnj8",####Pericyte
	"Isg15", "Ifit1","Ifit3","Ifit3b"####Interferon
	)
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))

dev.off()
new.cluster.ids <- c("Capillary_1","C_V_1","Interferon","C_V_2","C_A","Capillary_2","Arterial","Capillary_3",
  "Capillary_4","Venous","Choroid_plexus")
names(new.cluster.ids) <- levels(EC_cells.integrated)
EC_cells.integrated <- RenameIdents(EC_cells.integrated, new.cluster.ids)
EC_cells.integrated$celltype<-EC_cells.integrated@active.ident

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_EC_cells_annotation_UMAP.pdf")
DimPlot(EC_cells.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "Brain ECs subtypes")
dev.off();

pdf("/md01/nieyg/project/BBB/YMO_results/ECplot/BBB_onlyEC_subtype_dotplot.pdf",width=12,height=8)
features <- c( "Fbln5", "Cytl1","Mgp","S100a6", "Azin1","Pi16","Bmx", "Vegfc","Fbln2","Gkn3","Hey1","Edn3", ######Arterial
  "Tgfb2", "Glul","Slc26a10","Lypd1",###C-A
  "Ddc","Mfsd2a", "Cxcl12","Spock2","Rgcc",
  "Tfrc", "Car4", "Itm2a","Chn2",#####C-V
  "Lcn2","Slc38a5","Nr2f2", "Sox12","Icam1","Vcam1", "Vwf",##Venous
  "Plvap", "Plpp3","Esm1",####choroid plexus
  "Isg15", "Ifit1","Ifit3","Ifit3b"####Interferon
  )
DotPlot(EC_cells.integrated, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

saveRDS(EC_cells.integrated,"YMO_onlyEC_integrated.rds")
