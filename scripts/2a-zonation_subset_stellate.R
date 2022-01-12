library(Seurat)
library(tidyverse)

#############

####### ---- Prepare Seurat object to work with locally -------------------
metadata = readRDS("data/metaData_bigMouse_Robin.rds") ## use this metadata for the analysis on the PRISM cluster
counts = readRDS("data/rawCounts_bigMouse.rds")


metadata %>% filter(annot== "Stellate cells") %>% pull(digest) %>% table()
cells_oi = metadata %>% filter(annot %in% "Stellate cells" & digest != "GentleMacs") %>% pull(cell) %>% unique()

metadata = metadata[cells_oi,]
counts = counts[,cells_oi]

###

seurat_obj = CreateSeuratObject(counts = counts, project = "mouseLiverSCA", min.cells = 5, min.features = 500, meta.data = metadata)
rm(counts)
rm(metadata)

# Perform integration
liver.list <- SplitObject(seurat_obj, split.by = "digest")
rm(seurat_obj)
liver.list <- lapply(X = liver.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = liver.list) # nr of features cannot be higher...if not higher in previous thing!
liver.list <- PrepSCTIntegration(object.list = liver.list, anchor.features = features)
liver_all.anchors <- FindIntegrationAnchors(object.list = liver.list, normalization.method = "SCT", 
                                            anchor.features = features)
seurat_obj_stellate <- IntegrateData(anchorset = liver_all.anchors, normalization.method = "SCT")
seurat_obj_stellate <- RunPCA(seurat_obj_stellate, verbose = FALSE)

DimPlot(seurat_obj_stellate, group.by = "digest")


# Define which PCs most strongly correlated to the zonation
cor(seurat_obj_stellate@reductions$pca@cell.embeddings, as.matrix(Matrix::t(seurat_obj_stellate@assays$SCT@data[c("Ngfr","Igfbp3","Dach1"),]))) %>% 
  rowMeans() %>% abs() %>% sort(decreasing = T) 

dims_oi = cor(seurat_obj_stellate@reductions$pca@cell.embeddings, as.matrix(Matrix::t(seurat_obj_stellate@assays$SCT@data[c("Ngfr","Igfbp3","Dach1"),]))) %>% 
  rowMeans() %>% abs() %>% sort(decreasing = T)  %>% .[. > 0.09] %>% names() %>% gsub("PC_","",.) %>% as.double()

# Run UMAP and clustering based on those PCs
seurat_obj_stellate <- RunUMAP(seurat_obj_stellate, reduction = "pca", dims = dims_oi)
DimPlot(seurat_obj_stellate, split.by = "digest", pt.size = 0.85)

seurat_obj_stellate = seurat_obj_stellate %>% 
  FindNeighbors(reduction = "pca", dims = dims_oi) %>% 
  FindClusters(resolution = 0.2) %>% identity()

DimPlot(seurat_obj_stellate, split.by = "digest", pt.size = 0.85)

DefaultAssay(seurat_obj_stellate) = "SCT"

FeaturePlot(seurat_obj_stellate, c("Ngfr","Igfbp3","Dach1","Il34","Rgs4", "Spon2"), pt.size = 0.85)
DimPlot(seurat_obj_stellate, pt.size = 0.85)
VlnPlot(seurat_obj_stellate, c("Ngfr","Igfbp3","Dach1","Il34","Rgs4", "Spon2"))

# it seems that only cluster 0 is the central vein zone
seurat_obj_stellate@meta.data$zonation = "portal"
seurat_obj_stellate@meta.data$zonation[seurat_obj_stellate@meta.data$seurat_clusters == 0] = "central"

seurat_obj_stellate@meta.data$celltype = paste(seurat_obj_stellate@meta.data$annot, seurat_obj_stellate@meta.data$zonation, sep = "_")
seurat_obj_stellate = SetIdent(seurat_obj_stellate, value = "celltype")
DimPlot(seurat_obj_stellate, reduction = "umap", pt.size = 1)
VlnPlot(seurat_obj_stellate, c("Ngfr","Igfbp3","Dach1","Il34","Rgs4", "Spon2"))

allmarkers_stellate = FindAllMarkers(seurat_obj_stellate, only.pos = TRUE, logfc.threshold = 0.15)
allmarkers_stellate %>% filter(gene %in% c("Ngfr","Igfbp3","Dach1","Il34","Rgs4", "Spon2"))

allmarkers_stellate %>% group_by(cluster) %>% top_n(25, avg_log2FC) %>% top_n(6, -p_val_adj)

FeaturePlot(seurat_obj_stellate, c("Lsamp"), pt.size = 0.85)
FeaturePlot(seurat_obj_stellate, c("Adamtsl2"), pt.size = 0.85)
FeaturePlot(seurat_obj_stellate, c("Ngfr"), pt.size = 0.85)
FeaturePlot(seurat_obj_stellate, c("Igfbp3"), pt.size = 0.85)
FeaturePlot(seurat_obj_stellate, c("Dach1"), pt.size = 0.85)

### now export the zonation information

seurat_obj_stellate@meta.data %>% select(cell, digest, annot, celltype, zonation) %>% distinct() %>% write_csv2("output/stellate_zonation.csv")
seurat_obj_stellate@meta.data %>% select(cell, digest, annot, celltype, zonation) %>% distinct() %>% saveRDS("output/stellate_zonation_ids.rds")

### now get the zonation DE of all genes

DE_table_zonation = FindMarkers(object = seurat_obj_stellate, ident.1 = "Stellate cells_portal", ident.2 = "Stellate cells_central", min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = "SCT") %>% rownames_to_column("gene") %>% as_tibble()
# DE_table_zonation  = DE_table_zonation %>% mutate(sender = "Stellate cells", significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct)

DE_table_zonation %>% saveRDS("output/stellate_DE_table_zonation.rds")


