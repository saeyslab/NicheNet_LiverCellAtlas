library(Seurat)
library(tidyverse)
library(qsub)

path_seurat_obj = "/group/irc/personal/robinb/Liver_Visium/output/seurat_obj_scrnaseq.rds"
path_zonation = "/group/irc/personal/robinb/Liver_Visium/zonation_metadata.rds"
path_output = "/group/irc/personal/robinb/Liver_Visium/output/"

qsub_config = create_qsub_config(
  remote = "robinb@prism.psb.ugent.be:7777",
  local_tmp_path = "/home/robin/r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/robinb/r2gridengine",
  modules = "R/x86_64/4.0.3",
  memory = "215G",
  wait = FALSE,
  remove_tmp_folder = FALSE,
  name = "Liver-Zonation",
  max_wall_time = "500:00:00",
  stop_on_error = TRUE,
  num_cores = 1
)

celltype_id = "celltype"


include_zonation = function(i, path_seurat_obj, path_zonation, path_output, celltype_id){

  library(dplyr)
  library(tibble)
  library(Seurat)

  print("read in seurat object")
  seurat_obj = readRDS(path_seurat_obj)

  print("read in zonation")
  zonation_metadata  = readRDS(path_zonation)

  print("add zonation to the seuratobj")
  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

  non_zonated_metadata = seurat_obj@meta.data %>% as_tibble() %>% distinct(cell, digest, annot, celltype) %>% filter(! annot %in%  (zonation_metadata$annot %>% unique())) %>% mutate(zonation = NA)
  zonated_metadata = seurat_obj@meta.data %>% as_tibble() %>% distinct(cell, digest, annot) %>% filter(annot %in%  (zonation_metadata$annot %>% unique())) %>% left_join(zonation_metadata)

  new_metadata =  seurat_obj@meta.data %>% as_tibble() %>% select(-celltype) %>% inner_join(bind_rows(non_zonated_metadata, zonated_metadata))
  new_metadata_df = data.frame(new_metadata) %>% magrittr::set_rownames(new_metadata$cell)

  seurat_obj@meta.data = new_metadata_df

  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

  DefaultAssay(seurat_obj) = "SCT"

  seurat_obj %>% saveRDS(paste0(path_output, "seurat_obj_scrnaseq_zonation.rds"))
  p = DimPlot(seurat_obj)
  return(p)
}

set.seed(1)

job_liver_zonation = qsub_lapply(X = 1, FUN = include_zonation,
                            object_envir = environment(include_zonation),
                            qsub_config = qsub_config,
                            qsub_environment = NULL,
                            qsub_packages = NULL, path_seurat_obj, path_zonation, path_output, celltype_id)


saveRDS(job_liver_zonation, "output/job_liver_zonation.rds")
#
#
# #
job_liver_zonation = readRDS("output/job_liver_zonation.rds")
job_liver_zonation = qsub_retrieve(job_liver_zonation) %>% .[[1]]
#
#
saveRDS(job_liver_zonation, "output/umap_seurat_zonation.rds")
#

