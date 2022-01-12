library(Seurat)
library(tidyverse)
library(qsub)

path_metadata = "/group/irc/personal/robinb/Liver_Visium/metaData_bigMouse_Robin.rds"
path_counts = "/group/irc/personal/robinb/Liver_Visium/rawCounts_bigMouse.rds"
path_output = "/group/irc/personal/robinb/Liver_Visium/output/"

qsub_config = create_qsub_config(
  remote = "robinb@prism.psb.ugent.be:7777",
  local_tmp_path = "/home/robin/r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/robinb/r2gridengine",
  modules = "R/x86_64/4.0.3",
  memory = "180G",
  wait = FALSE,
  remove_tmp_folder = FALSE,
  name = "Liver-SCA",
  max_wall_time = "500:00:00",
  stop_on_error = TRUE,
  num_cores = 1
)

process_seurat = function(i, path_counts, path_metadata, path_output){

  library(dplyr)
  library(tibble)
  library(Seurat)

  # Read in scRNAseq data and perform SCT integration
  print("read in metadata and counts")
  metadata = readRDS(path_metadata)
  rownames(metadata) = metadata$cell

  cells_oi = metadata %>% filter(annot %in% "LSECs") %>% pull(cell) %>% unique()

  metadata = metadata[cells_oi,]

  counts = readRDS(path_counts)
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
  seurat_obj_LSEC <- IntegrateData(anchorset = liver_all.anchors, normalization.method = "SCT")
  seurat_obj_LSEC <- RunPCA(seurat_obj_LSEC, verbose = FALSE)
  DimPlot(seurat_obj_LSEC, group.by = "digest")

   # Define which PCs most strongly correlated to the zonation
  cor(seurat_obj_LSEC@reductions$pca@cell.embeddings, as.matrix(Matrix::t(seurat_obj_LSEC@assays$SCT@data[c("Mecom","Msr1","Efnb2"),]))) %>%
    rowMeans() %>% abs() %>% sort(decreasing = T)

  # dims_oi = cor(seurat_obj_LSEC@reductions$pca@cell.embeddings, as.matrix(Matrix::t(seurat_obj_LSEC@assays$SCT@data[c("Mecom","Msr1","Efnb2"),]))) %>%
  #   rowMeans() %>% abs() %>% sort(decreasing = T)  %>% .[. > 0.10] %>% names() %>% gsub("PC_","",.) %>% as.double()

  dims_oi = cor(seurat_obj_LSEC@reductions$pca@cell.embeddings, as.matrix(Matrix::t(seurat_obj_LSEC@assays$SCT@data[c("Mecom","Msr1","Efnb2"),]))) %>%
    rowMeans() %>% abs() %>% sort(decreasing = T)  %>% .[. > 0.15] %>% names() %>% gsub("PC_","",.) %>% as.double()
  # Run UMAP and clustering based on those PCs
  seurat_obj_LSEC <- RunUMAP(seurat_obj_LSEC, reduction = "pca", dims = dims_oi)
  DimPlot(seurat_obj_LSEC, split.by = "digest", pt.size = 0.85)


  DefaultAssay(seurat_obj_LSEC) = "integrated"

  seurat_obj_LSEC = seurat_obj_LSEC %>%
    FindNeighbors(reduction = "pca", dims = dims_oi) %>%
    FindClusters(resolution = 0.10) %>% identity()

  DefaultAssay(seurat_obj_LSEC) = "SCT"

  allmarkers_LSEC = FindAllMarkers(seurat_obj_LSEC, only.pos = TRUE, logfc.threshold = 0.25)
  allmarkers_LSEC %>% filter(gene %in% c("Mecom","Msr1","Efnb2") & p_val_adj < 0.05)
  cluster_portal = allmarkers_LSEC %>% filter(gene %in% c("Mecom","Msr1","Efnb2") & p_val_adj < 0.05) %>% pull(cluster) %>% unique()

  left_over_idents = seurat_obj_LSEC %>% Idents() %>% setdiff(cluster_portal)

  if(length(left_over_idents) > 1){
    allmarkers_LSEC_filtered = FindAllMarkers(seurat_obj_LSEC %>% subset(idents = left_over_idents) , only.pos = TRUE, logfc.threshold = 0.25)
    allmarkers_LSEC_filtered %>% filter(gene %in% c("Mecom","Msr1","Efnb2") & p_val_adj < 0.05)

    cluster_portal = union(cluster_portal, allmarkers_LSEC_filtered %>% filter(gene %in% c("Mecom","Msr1","Efnb2") & p_val_adj < 0.05) %>% pull(cluster) %>% unique())
  }

  cells_portal = seurat_obj_LSEC %>% subset(idents = cluster_portal) %>% Cells()
  cells_central = seurat_obj_LSEC %>% Cells() %>% setdiff(cells_portal)

  zonation = seurat_obj_LSEC@meta.data$cell %in% cells_portal
  zonation[zonation == TRUE] = "portal"
  zonation[zonation == FALSE] = "central"

  seurat_obj_LSEC@meta.data$zonation = zonation

  seurat_obj_LSEC@meta.data$celltype = paste(seurat_obj_LSEC@meta.data$annot, seurat_obj_LSEC@meta.data$zonation, sep = "_")
  seurat_obj_LSEC = SetIdent(seurat_obj_LSEC, value = "celltype")

  allmarkers_LSEC_final = FindAllMarkers(seurat_obj_LSEC, only.pos = TRUE, logfc.threshold = 0.25)

   ### now get the zonation DE of all genes

  DE_table_zonation = FindMarkers(object = seurat_obj_LSEC, ident.1 = "LSECs_portal", ident.2 = "LSECs_central", min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = "SCT") %>% rownames_to_column("gene") %>% as_tibble()

  return(list(zonation_metadata = seurat_obj_LSEC@meta.data %>% select(cell, digest, annot, celltype, zonation) %>% distinct(), DE_table_zonation = DE_table_zonation))

}

set.seed(1)

job_liver_sca = qsub_lapply(X = 1, FUN = process_seurat,
                            object_envir = environment(process_seurat),
                            qsub_config = qsub_config,
                            qsub_environment = NULL,
                            qsub_packages = NULL, path_counts, path_metadata, path_output)


saveRDS(job_liver_sca, "output/job_liver_sca_lsecs.rds")


job_liver_sca = readRDS("output/job_liver_sca_lsecs.rds")
output_liver_de = qsub_retrieve(job_liver_sca) %>% .[[1]]

saveRDS(output_liver_de, "output/output_liver_zonation_lsecs.rds")

output_liver_de$DE_table_zonation %>% arrange(-avg_log2FC)
output_liver_de$DE_table_zonation %>% arrange(avg_log2FC)


