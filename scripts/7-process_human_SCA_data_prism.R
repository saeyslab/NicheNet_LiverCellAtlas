library(Seurat)
library(tidyverse)
library(qsub)

path_metadata = "/group/irc/personal/robinb/Liver_Visium/metaData_bigHuman_Robin.rds"
path_counts = "/group/irc/personal/robinb/Liver_Visium/rawCounts_bigHuman.rds"
path_output = "/group/irc/personal/robinb/Liver_Visium/output/"

qsub_config = create_qsub_config(
  remote = "robinb@prism.psb.ugent.be:7777",
  local_tmp_path = "/home/robin/r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/robinb/r2gridengine",
  modules = "R/x86_64/4.0.3",
  memory = "215G",
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
  # metadata = metadata %>% rownames_to_column("cell")
  rownames(metadata) = metadata$cell
  
  counts = readRDS(path_counts)
  counts = counts[,metadata$cell]
  
  print("make scrnaseq seurat obj")
  
  seurat_obj = CreateSeuratObject(counts = counts, project = "mouseLiverSCA", min.cells = 5, min.features = 500, meta.data = metadata)
  rm(counts)
  rm(metadata)
  
  print(table(seurat_obj@meta.data$annot, seurat_obj@meta.data$digest))
  
  
  
  print("perform sct integration")
  
  liver.list <- SplitObject(seurat_obj, split.by = "type")
  rm(seurat_obj)
  liver.list <- lapply(X = liver.list, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = liver.list, nfeatures = 5000)
  liver.list <- PrepSCTIntegration(object.list = liver.list, anchor.features = features)
  liver_all.anchors <- FindIntegrationAnchors(object.list = liver.list, normalization.method = "SCT", 
                                              anchor.features = features)
  seurat_obj_scrnaseq <- IntegrateData(anchorset = liver_all.anchors, normalization.method = "SCT")
  
  print("run DR")
  
  seurat_obj_scrnaseq <- RunPCA(seurat_obj_scrnaseq, verbose = FALSE)
  seurat_obj_scrnaseq <- RunUMAP(seurat_obj_scrnaseq, reduction = "pca", dims = 1:30)
  p1 <- DimPlot(seurat_obj_scrnaseq, reduction = "umap", group.by = "digest")
  p2 <- DimPlot(seurat_obj_scrnaseq, reduction = "umap", group.by = "annot", label = TRUE, 
                repel = TRUE)
  plot = p1 + p2
  
  rm(liver_all.anchors)
  
  seurat_obj_scrnaseq = seurat_obj_scrnaseq %>% SetIdent(value = "annot")
  seurat_obj_scrnaseq@meta.data$celltype = seurat_obj_scrnaseq@meta.data$annot

  if(DefaultAssay(seurat_obj_scrnaseq) != "SCT"){
    DefaultAssay(seurat_obj_scrnaseq) = "SCT"
  }
  seurat_obj_scrnaseq = seurat_obj_scrnaseq %>% SetIdent(value = "annot")
  seurat_obj_scrnaseq@meta.data$celltype = seurat_obj_scrnaseq@meta.data$annot
  seurat_obj_scrnaseq@meta.data$celltype[seurat_obj_scrnaseq@meta.data$celltype == "immLAMs"] = "LAMs"
  seurat_obj_scrnaseq@meta.data$celltype[seurat_obj_scrnaseq@meta.data$celltype == "matLAMs"] = "LAMs"
  seurat_obj_scrnaseq = seurat_obj_scrnaseq %>% SetIdent(value = "celltype")
  
  seurat_obj_scrnaseq %>% saveRDS(paste0(path_output, "seurat_obj_scrnaseq_human.rds"))
  

  a = NULL
  return(a)
  
}

set.seed(1)

job_liver_sca = qsub_lapply(X = 1, FUN = process_seurat,
                                   object_envir = environment(process_seurat),
                                   qsub_config = qsub_config,
                                   qsub_environment = NULL,
                                   qsub_packages = NULL, path_counts, path_metadata, path_output)


saveRDS(job_liver_sca, "output/job_liver_sca_human.rds")





