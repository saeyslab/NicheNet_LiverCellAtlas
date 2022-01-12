library(Seurat)
library(tidyverse)

metadata = readRDS("data/metaData_bigHuman.rds")
metadata = metadata %>% rownames_to_column("cell")

CD45_neg_annotations = read_delim("data/CD45neg_human.csv",delim = ",")
Myeloid_annotations = read_delim("data/Myeloid_human.csv",delim = ",")

CD45_neg_annotations = CD45_neg_annotations %>% rename(cell=barcodes) %>% select(cell, type, annot) %>% distinct()
Myeloid_annotations = Myeloid_annotations  %>% select(cell, type, annot) %>% distinct()

CD45_neg_annotations$type[CD45_neg_annotations$type == "rnaSeq"] = "RnaSeq"
CD45_neg_annotations$type[CD45_neg_annotations$type == "nuqSeq"] = "NucSeq"

new_annotations = Myeloid_annotations %>% bind_rows(CD45_neg_annotations) %>% distinct()

older_annotations = metadata %>% rename(annot_old = annot) %>% as_tibble()
annotations_combined = older_annotations %>% left_join(new_annotations)

annotations_combined$annot %>% table() %>% sort(decreasing = TRUE)
annotations_combined$annot_old %>% table() %>% sort(decreasing = TRUE)

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

annotations_combined = annotations_combined %>% mutate_cond(is.na(annot), annot = annot_old)
annotations_combined$annot %>% table() %>% sort(decreasing = TRUE)

annotations_combined %>% filter(annot != annot_old) %>% select(annot, annot_old) %>% group_by(annot, annot_old) %>% count() %>% arrange(-n) %>% View()

table(annotations_combined$annot, annotations_combined$type)

metadata = annotations_combined %>% data.frame() %>% magrittr::set_rownames(annotations_combined$cell)
metadata %>% saveRDS("data/metaData_bigHuman_Robin.rds") ## use this metadata for the analysis on the PRISM cluster

####### ---- Prepare Seurat object to work with locally -------------------

annotations_combined = annotations_combined %>% filter(annot %in% c(
  "immLAMs","matLAMs","Cholangiocytes","Fibroblasts","Hepatocytes","resKCs","LSECs", "Stellate cells", "Central Vein Endothelial cells", "MoMac1","Macrophages","Portal Vein Endothelial cells", "Pre-moKCs and moKCs"
))

sampling_metadata_frac = annotations_combined %>% distinct(cell, annot, typePatient, type) %>%
  group_by(annot, type, typePatient) %>%
  slice_sample(prop = 0.10) %>% as_tibble()

sampling_metadata_n = annotations_combined %>% distinct(cell, annot, typePatient, type) %>%
  group_by(annot, type, typePatient) %>%
  slice_sample(n = 500) %>% as_tibble()

sampling_metadata = bind_rows(sampling_metadata_frac, sampling_metadata_n) %>% distinct()

counts = readRDS("data/rawCounts_bigHuman.rds")
metadata = readRDS("data/metaData_bigHuman_Robin.rds")

metadata = annotations_combined %>% data.frame() %>% magrittr::set_rownames(annotations_combined$cell)
metadata = metadata[sampling_metadata$cell,]
counts = counts[,sampling_metadata$cell]

###

seurat_obj = CreateSeuratObject(counts = counts, project = "HumanLiverSCA", min.cells = 5, min.features = 500, meta.data = metadata)

# Again: now with integration
seurat_obj@meta.data$type_patient = paste0(seurat_obj@meta.data$type, seurat_obj@meta.data$typePatient)
seurat_obj@meta.data$type_patient %>% table()
seurat_obj@meta.data$typePatient  = factor(seurat_obj@meta.data$typePatient)
seurat_obj@meta.data$type  = factor(seurat_obj@meta.data$type)

liver.list <- SplitObject(seurat_obj, split.by = "type_patient")
# liver.list <- SplitObject(seurat_obj, split.by = "type")

liver.list <- lapply(X = liver.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = liver.list, nfeatures = 5000)
liver.list <- PrepSCTIntegration(object.list = liver.list, anchor.features = features)
liver_all.anchors <- FindIntegrationAnchors(object.list = liver.list, normalization.method = "SCT",
                                            anchor.features = features)
liver_all.combined.sct <- IntegrateData(anchorset = liver_all.anchors, normalization.method = "SCT", dims = 1:50, k.weight = 50)
liver_all.combined.sct <- RunPCA(liver_all.combined.sct, verbose = FALSE)
liver_all.combined.sct <- RunUMAP(liver_all.combined.sct, reduction = "pca", dims = 1:50)
p1 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "type")
p2 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "annot", label = TRUE,
              repel = TRUE)
p3 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "typePatient", label = TRUE,
              repel = TRUE)
p1 + p2 + p3

DefaultAssay(liver_all.combined.sct) = "SCT"
VlnPlot(liver_all.combined.sct, "SDS")
VlnPlot(liver_all.combined.sct, "HAL")

liver_all.combined.sct = liver_all.combined.sct %>% SetIdent(value = "annot")
liver_all.combined.sct@meta.data$celltype = liver_all.combined.sct@meta.data$annot
liver_all.combined.sct@meta.data$celltype[liver_all.combined.sct@meta.data$celltype == "immLAMs"] = "LAMs"
liver_all.combined.sct@meta.data$celltype[liver_all.combined.sct@meta.data$celltype == "matLAMs"] = "LAMs"
liver_all.combined.sct = liver_all.combined.sct %>% SetIdent(value = "celltype")
p1 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "type")
p2 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "celltype", label = TRUE,
              repel = TRUE)
p3 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "typePatient", label = TRUE,
              repel = TRUE)
p1 + p2 + p3
saveRDS(liver_all.combined.sct, "data/seurat_obj_subset_integrated_human.rds")
