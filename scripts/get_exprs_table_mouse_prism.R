library(Seurat)
library(tidyverse)
library(nichenetr)
library(qsub)
source("scripts/plotting_functions.R")
source("scripts/functions_differential_nichenet.R")
expression_pct = 0.10
lfc_cutoff = 0.15
assay_oi = "SCT"

path_seurat_obj = "/group/irc/personal/robinb/Liver_Visium/output/seurat_obj_scrnaseq_zonation.rds"
qsub_config = create_qsub_config(
  remote = "robinb@prism.psb.ugent.be:7777",
  local_tmp_path = "/home/robin/r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/robinb/r2gridengine",
  modules = "R/x86_64/4.0.3",
  memory = "128G",
  wait = FALSE,
  remove_tmp_folder = FALSE,
  name = "KC-ALL-M",
  max_wall_time = "500:00:00",
  stop_on_error = TRUE,
  num_cores = 1
)

get_exprs_table = function(i, path_seurat_obj, features_oi, celltypes_oi){
  
  library(nichenetr)
  library(tidyverse)
  library(Seurat)
  
  print("read in seurat object")
  seurat_obj = readRDS(path_seurat_obj)
  
  # celltypes_to_filter = niches %>% unlist() %>% unique() %>% union( spatial_info$celltype_other_region %>% unique())
  celltype_id = "celltype" # metadata column
  
  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])
  
  # ---------- Get expression information of all genes ---------------------------------------
  
  ## Get exprs information and fractions - combine them with DE information and ligand activities

  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = celltypes_oi), features = features_oi, assay = assay_oi))
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  return(exprs_tbl)
}

set.seed(1)
genes = readRDS("output/genes_oi_exprs_table_mouse.rds")
celltypes = readRDS("output/celltypes_oi_exprs_table_mouse.rds")

celltypes[celltypes == "Stellate cells"] = "Stellate cells_portal"
celltypes[celltypes == "Hepatocytes"] = "Hepatocytes_portal"
celltypes[celltypes == "LSECs"] = "LSECs_portal"


job_liver_de = qsub_lapply(X = 1, FUN = get_exprs_table,
                           object_envir = environment(get_exprs_table),
                           qsub_config = qsub_config,
                           qsub_environment = NULL,
                           qsub_packages = NULL, path_seurat_obj, genes, celltypes)


saveRDS(job_liver_de, "output/job_exprs_table_mouse.rds")

output = readRDS("output/job_exprs_table_mouse.rds")

output = qsub_retrieve(output)

saveRDS(output %>% .[[1]], "output/exprs_table_mouse.rds")



