library(Seurat)
library(tidyverse)
library(nichenetr)
library(qsub)
source("scripts/plotting_functions.R")
source("scripts/functions_differential_nichenet.R")
expression_pct = 0.10
lfc_cutoff = 0.15
top_n_target = 250
assay_oi = "SCT"
specificity_score_LR_pairs = "min_lfc"
specificity_score_targets = "min_lfc"
specificity_score_spatial = "lfc"

path_seurat_obj = "/group/irc/personal/robinb/Liver_Visium/output/seurat_obj_scrnaseq_zonation.rds"
qsub_config = create_qsub_config(
  remote = "robinb@prism.psb.ugent.be:7777",
  local_tmp_path = "/home/robin/r2gridengine",
  remote_tmp_path = "/scratch/irc/personal/robinb/r2gridengine",
  modules = "R/x86_64/4.0.3",
  memory = "128G",
  wait = FALSE,
  remove_tmp_folder = FALSE,
  name = "MoMac1",
  max_wall_time = "500:00:00",
  stop_on_error = TRUE,
  num_cores = 1
)

differential_nichenet = function(i, path_seurat_obj){

  library(nichenetr)
  library(tidyverse)
  library(Seurat)

  print("read in seurat object")
  seurat_obj = readRDS(path_seurat_obj)

  niches = list(
    "KC_niche" = list(
      "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
      "receiver" = c("KCs")),
    "MoMac1_niche" = list(
      "sender" = c("Capsule fibroblasts","Mesothelial cells"),
      "receiver" = c("MoMac1"))
  )

  # ---------- Define spatial information (if applicable) ------------------------------------------------------------

  spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"), celltype_other_region = c(c("LSECs_central","Hepatocytes_central","Stellate cells_central"))) %>% mutate(niche =  "KC_niche", celltype_type = "sender")

  # ---------- Filter Seurat object to only contain celltypes of interest (optional - but can be useful for requiring less memory) ------------------------------------------------------------

  # celltypes_to_filter = niches %>% unlist() %>% unique() %>% union( spatial_info$celltype_other_region %>% unique())
  celltype_id = "celltype" # metadata column

  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

  # ---------- Define parameters of the Differential NicheNet pipeline ------------------------------------------------------------

  expression_pct = 0.10
  lfc_cutoff = 0.15
  top_n_target = 250
  assay_oi = "SCT"
  specificity_score_LR_pairs = "min_lfc"
  specificity_score_targets = "min_lfc"
  specificity_score_spatial = "lfc"
  # ---------- Read in NicheNet LR network and ligand activities ------------------------------------------------------------

  # change human to mouse symbols if required. Do it yourself, or use the following code.

  # LR network
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
  lr_network = lr_network %>% mutate(ligand = nichenetr::convert_human_to_mouse_symbols(ligand), receptor = nichenetr::convert_human_to_mouse_symbols(receptor)) %>% drop_na()

  # Ligand-target matrix
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

  # ---------- Determine DE between the different niches for both senders and receivers to define the DE L-R pairs ------------------------------------------------------------

  # Calculate DE
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj, niches = niches, type = "sender", assay_oi = assay_oi)
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj, niches = niches, type = "receiver", assay_oi = assay_oi)

  # Process DE
  DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, type = "receiver")

  # Combine sender-receiver DE based on LR pairs
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

  # ---------- Determine target gene set of interest of the receiver cell types & predict ligand activities and ligand-target links ------------------------------------------------------------
  # you can give your geneset of interest in a different way!

  compare_vs_MoMac1 = TRUE

  if(compare_vs_MoMac1 == TRUE){
    niches_targets = list(
      "KC_niche" = list(
        "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
        "receiver" = c("KCs")),
      "MoMac1_niche" = list(
        "sender" = c("Capsule fibroblasts","Mesothelial cells"),
        "receiver" = c("MoMac1")),
      "MoMac2_niche" = list(
        "sender" = c("Cholangiocytes","Fibroblast 2"),
        "receiver" = c("MoMac2"))
    )

    DE_receiver_targets = calculate_niche_de(seurat_obj = seurat_obj, niches = niches_targets, type = "receiver", assay_oi = assay_oi)
    DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches_targets, specificity_score = specificity_score_targets)

  } else {
    DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver, niches = niches, specificity_score = specificity_score_targets)

  }

  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_KC = DE_receiver_processed_targets %>% filter(receiver == "KCs" & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_MoMac1 = DE_receiver_processed_targets %>% filter(receiver == "MoMac1" & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

  # Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
  # If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
  geneset_KC %>% setdiff(rownames(ligand_target_matrix))
  geneset_MoMac1 %>% setdiff(rownames(ligand_target_matrix))

  niche_geneset_list = list(
    "KC_niche" = list(
      "receiver" = "KCs",
      "geneset" = geneset_KC,
      "background" = background),
    "MoMac1_niche" = list(
      "receiver" = "MoMac1",
      "geneset" = geneset_MoMac1 ,
      "background" = background)
  )

  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

  # ---------- Calculate Zonation DE ------------------------------------------------------------

  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj, spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, specificity_score = specificity_score_spatial)

  # add a neutral zonation score for sender celltypes in which the zonation is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)

  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_zonation = scale_quantile_adapted(ligand_score_zonation))

  # add a neutral zonation score for all receiver celltypes (for none of them, zonation is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)

  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_zonation = scale_quantile_adapted(receptor_score_zonation))

  # ---------- Get expression information of all genes ---------------------------------------

  ## Get exprs information and fractions - combine them with DE information and ligand activities
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target)

  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
  exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

  exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  # ---------- Score ligand-receptor interactions based on expression strength of the receptor + whether the interaction is bona fide or not ---------------------------------------
  # give higher scores to the most strongly expressed receptor of a certain ligand, in a certain celltype

  exprs_sender_receiver = lr_network %>%
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>%
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup()

  return(list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
         ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target))

}

set.seed(1)

job_liver_de = qsub_lapply(X = 1, FUN = differential_nichenet,
                           object_envir = environment(differential_nichenet),
                           qsub_config = qsub_config,
                           qsub_environment = NULL,
                           qsub_packages = NULL, path_seurat_obj)


saveRDS(job_liver_de, "output/job_differential_nichenet_mouse_MoMac1.rds")

output = readRDS("output/job_differential_nichenet_mouse_MoMac1.rds")

output = qsub_retrieve(output)

saveRDS(output %>% .[[1]], "output/differential_nichenet_mouse_KC_MoMac1_new.rds")



