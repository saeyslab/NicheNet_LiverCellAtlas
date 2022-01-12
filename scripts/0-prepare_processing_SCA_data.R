library(Seurat)
library(tidyverse)

metadata = readRDS("data/metaData_bigMouse.rds")

CD45_neg_annotations = xlsx::read.xlsx2("data/CD45neg_new.xlsx", sheetIndex = 1,header = TRUE)
Fibroblasts_annotations = xlsx::read.xlsx2("data/Fibroblasts_new.xlsx", sheetIndex = 1,header = TRUE)
Myeloid_annotations = xlsx::read.xlsx2("data/Myeloid_new.xlsx", sheetIndex = 1,header = TRUE)

CD45_neg_annotations = CD45_neg_annotations %>% filter(!cell %in% Fibroblasts_annotations$cell)

new_annotations = Myeloid_annotations %>% bind_rows(CD45_neg_annotations) %>% bind_rows(Fibroblasts_annotations) %>% distinct()

older_annotations = metadata %>% rownames_to_column("cell") %>% rename(annot_old = annot) %>% as_tibble()
annotations_combined = older_annotations %>% left_join(new_annotations %>% select(-sampleName))

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

annotations_combined = annotations_combined %>% mutate_cond(is.na(annot), annot = annot_old)
annotations_combined = annotations_combined %>% mutate_cond(annot == "Hepatocytes NucSeq", annot = "Hepatocytes")
annotations_combined = annotations_combined %>% mutate_cond(annot == "Kupffer cells", annot = "KCs")

table(annotations_combined$annot, annotations_combined$digest)
table(annotations_combined$annot) %>% sort(decreasing = TRUE)

metadata = annotations_combined %>% data.frame() %>% magrittr::set_rownames(annotations_combined$cell)
metadata %>% saveRDS("data/metaData_bigMouse_Robin.rds") ## use this metadata for the analysis on the PRISM cluster


####### ---- Prepare Seurat object to work with locally -------------------

annotations_combined = annotations_combined %>% filter(annot %in% c(
  "Capsule fibroblasts","Central Vein Endothelial cells","Cholangiocytes","Fibroblast 1","Fibroblast 2","Hepatocytes","KCs","LSECs", "Mesothelial cells","MoMac1","MoMac2", "Stellate cells"
))

sampling_metadata_frac = annotations_combined %>% distinct(cell, annot, digest) %>%
  group_by(annot, digest) %>%
  slice_sample(prop = 0.05) %>% as_tibble()

sampling_metadata_n = annotations_combined %>% distinct(cell, annot, digest) %>%
  group_by(annot) %>%
  slice_sample(n = 500) %>% as_tibble()

sampling_metadata = bind_rows(sampling_metadata_frac, sampling_metadata_n) %>% distinct()

counts = readRDS("data/rawCounts_bigMouse.rds")
metadata = readRDS("data/metaData_bigMouse_Robin.rds")

metadata = annotations_combined %>% data.frame() %>% magrittr::set_rownames(annotations_combined$cell)
metadata = metadata[sampling_metadata$cell,]
counts = counts[,sampling_metadata$cell]

###

seurat_obj = CreateSeuratObject(counts = counts, project = "mouseLiverSCA", min.cells = 5, min.features = 500, meta.data = metadata)
rm(counts)
rm(metadata)
rm(annotations_combined)
rm(new_annotations)
rm(older_annotations)
rm(sampling_metadata)
rm(sampling_metadata_frac)
rm(sampling_metadata_n)

# Again: now with integration
liver.list <- SplitObject(seurat_obj, split.by = "digest")
rm(seurat_obj)
liver.list <- lapply(X = liver.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = liver.list, nfeatures = 5000)
liver.list <- PrepSCTIntegration(object.list = liver.list, anchor.features = features)
liver_all.anchors <- FindIntegrationAnchors(object.list = liver.list, normalization.method = "SCT",
                                            anchor.features = features)
liver_all.combined.sct <- IntegrateData(anchorset = liver_all.anchors, normalization.method = "SCT")
liver_all.combined.sct <- RunPCA(liver_all.combined.sct, verbose = FALSE)
liver_all.combined.sct <- RunUMAP(liver_all.combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "digest")
p2 <- DimPlot(liver_all.combined.sct, reduction = "umap", group.by = "annot", label = TRUE,
              repel = TRUE)
p1 + p2

rm(liver_all.anchors)

liver_all.combined.sct = liver_all.combined.sct %>% SetIdent(value = "annot")
liver_all.combined.sct@meta.data$celltype = liver_all.combined.sct@meta.data$annot

saveRDS(liver_all.combined.sct, "data/seurat_obj_subset_integrated.rds")
