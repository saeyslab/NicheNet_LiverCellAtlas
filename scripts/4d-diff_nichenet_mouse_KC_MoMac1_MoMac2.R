library(tidyverse)
library(Seurat)
library(nichenetr)
library(circlize)

source("scripts/plotting_functions.R")
source("scripts/functions_differential_nichenet.R")

lfc_cutoff = 0.15
# Ligand-target matrix
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

expression_pct = 0.10

prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 2,
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "scaled_avg_score_ligand_receptor" = 0,
                         "scaled_ligand_score_zonation" = 2,
                         "scaled_receptor_score_zonation" = 0,
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_activity" = 0,"scaled_activity_normalized" = 1,
                         "ligand_fraction" = 1, "receptor_fraction" = 1,
                         "bona_fide" = 1)


####################################### ####################################### #######################################
####################################### MOUSE MOMAC2 ####################################### ###############################
####################################### ####################################### #######################################

output_nichenet_analysis_MoMac2 = readRDS("output/differential_nichenet_mouse_KC_MoMac2_new.rds")
prioritization_tables_MoMac2 = get_prioritization_tables(output_nichenet_analysis_MoMac2, prioritizing_weights)

prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs")
prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs")

####################################### ####################################### #######################################
####################################### MOUSE MOMAC1 ####################################### ###############################
####################################### ####################################### #######################################

output_nichenet_analysis_MoMac1 = readRDS("output/differential_nichenet_mouse_KC_MoMac1_new.rds")
prioritization_tables_MoMac1 = get_prioritization_tables(output_nichenet_analysis_MoMac1, prioritizing_weights)

prioritization_tables_MoMac1$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs")
prioritization_tables_MoMac1$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs")

####################################### ####################################### #######################################
####################################### MOUSE MOMAC1-CV ####################################### ###############################
####################################### ####################################### #######################################

output_nichenet_analysis_MoMac1_CV = readRDS("output/differential_nichenet_mouse_KC_MoMac1_CV.rds")
prioritization_tables_MoMac1_CV = get_prioritization_tables(output_nichenet_analysis_MoMac1_CV, prioritizing_weights)

prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs")
prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs")


####################################### ####################################### #######################################
####################################### Combine output in one table for visualization ####################################### ###############################
####################################### ####################################### #######################################

####### Ligand-Receptor #######


# KC
combined_ligand_receptor_tbl_KC = list(
  prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs"),
  prioritization_tables_MoMac1$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs"),
  prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs")
) %>% bind_rows() %>% mutate(id = paste0(sender, receiver, ligand_receptor))

combined_ligand_receptor_tbl_KC = combined_ligand_receptor_tbl_KC %>%
  mutate_cond(sender == "Stellate cells_portal", sender = "Stellate cells") %>%
  mutate_cond(sender == "Hepatocytes_portal", sender = "Hepatocytes") %>%
  mutate_cond(sender == "LSECs_portal", sender = "LSECs")


# check now the receptors expressed by KCs in mouse, and the ones higher in KCs than other MFs
prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs") %>% filter(receptor_fraction > 0.05) %>% distinct(receptor, receptor_fraction) %>% arrange(-receptor_fraction) %>% xlsx::write.xlsx2("output/KC_receptors.xlsx",sheetName = "mouse_expressed", append = FALSE)

inner_join(
  prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs") %>% filter(receptor_fraction > 0.05) %>% distinct(receptor, receptor_score) %>% rename(lfc_vs_LAM = receptor_score),
  prioritization_tables_MoMac1$prioritization_tbl_ligand_receptor %>% filter(receiver == "KCs") %>% filter(receptor_fraction > 0.05) %>% distinct(receptor, receptor_score) %>% rename(lfc_vs_CD207MF = receptor_score)
  ) %>% mutate(lfc_avg = (lfc_vs_LAM + lfc_vs_CD207MF)/2) %>% arrange(-lfc_avg) %>% filter(lfc_vs_LAM > 0 | lfc_vs_CD207MF > 0) %>% xlsx::write.xlsx2("output/KC_receptors.xlsx",sheetName = "mouse_DE", append = TRUE)

combined_ligand_receptor_tbl_KC_summarized = combined_ligand_receptor_tbl_KC %>% select(id, ligand_score:prioritization_score) %>% distinct() %>% group_by(id) %>% summarise_all(mean, na.rm = TRUE)
combined_ligand_receptor_tbl_KC_summarized = combined_ligand_receptor_tbl_KC %>% distinct(id, niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide) %>% inner_join(combined_ligand_receptor_tbl_KC_summarized) %>% select(-id) %>% arrange(-prioritization_score)

# MFs

combined_ligand_receptor_tbl_MF = list(
  prioritization_tables_MoMac2$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs"),
  prioritization_tables_MoMac1$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs"),
  prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_receptor %>% filter(receiver != "KCs")
) %>% bind_rows() %>% distinct()

# combine

prioritization_tbl_ligand_receptor = bind_rows(combined_ligand_receptor_tbl_KC_summarized, combined_ligand_receptor_tbl_MF)
# filtering
quantile_cutoff = prioritization_tbl_ligand_receptor %>% distinct(sender, receiver, ligand, receptor, prioritization_score) %>% pull(prioritization_score) %>% quantile(0.95,na.rm = TRUE)
prioritization_tbl_ligand_receptor_filtered = prioritization_tbl_ligand_receptor %>% filter(prioritization_score > quantile_cutoff)

####### Ligand-Target #######

# KC
combined_ligand_target_tbl_KC = list(
  prioritization_tables_MoMac2$prioritization_tbl_ligand_target %>% filter(receiver == "KCs"),
  prioritization_tables_MoMac1$prioritization_tbl_ligand_target %>% filter(receiver == "KCs"),
  prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_target %>% filter(receiver == "KCs")
) %>% bind_rows() %>% mutate(id = paste0(sender, receiver, ligand_receptor, target))

combined_ligand_target_tbl_KC = combined_ligand_target_tbl_KC %>%
  mutate_cond(sender == "Stellate cells_portal", sender = "Stellate cells") %>%
  mutate_cond(sender == "Hepatocytes_portal", sender = "Hepatocytes") %>%
  mutate_cond(sender == "LSECs_portal", sender = "LSECs")

combined_ligand_target_tbl_KC_summarized = combined_ligand_target_tbl_KC %>% select(id, target_score:prioritization_score) %>% distinct() %>% group_by(id) %>% summarise_all(mean, na.rm = TRUE)
combined_ligand_target_tbl_KC_summarized = combined_ligand_target_tbl_KC %>% distinct(id, niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target) %>% inner_join(combined_ligand_target_tbl_KC_summarized) %>% select(-id) %>% arrange(-prioritization_score)

# MFs

combined_ligand_target_tbl_MF = list(
  prioritization_tables_MoMac2$prioritization_tbl_ligand_target %>% filter(receiver != "KCs"),
  prioritization_tables_MoMac1$prioritization_tbl_ligand_target %>% filter(receiver != "KCs"),
  prioritization_tables_MoMac1_CV$prioritization_tbl_ligand_target %>% filter(receiver != "KCs")
) %>% bind_rows() %>% distinct()

# combine

prioritization_tbl_ligand_target = bind_rows(combined_ligand_target_tbl_KC_summarized, combined_ligand_target_tbl_MF)

####################################### ####################################### #######################################
####################################### Visualization ####################################### #########################
####################################### ####################################### #######################################

prioritization_tbl_ligand_receptor = prioritization_tbl_ligand_receptor %>% mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac1_CV_niche","MoMac2_niche"))) %>% mutate(niche = gsub("_niche","",niche))
prioritization_tbl_ligand_receptor_filtered = prioritization_tbl_ligand_receptor_filtered %>% mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac1_CV_niche","MoMac2_niche")))  %>% mutate(niche = gsub("_niche","",niche))
prioritization_tbl_ligand_target = prioritization_tbl_ligand_target %>% mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac1_CV_niche","MoMac2_niche"))) %>% mutate(niche = gsub("_niche","",niche))

### !!!! for filtering: prioritization score is not allowed to be higher in the opposite niche!!
### disadvantage: risk to lose things that could score high for both niches

prioritization_tbl_ligand_receptor_spread = prioritization_tbl_ligand_receptor %>% select(ligand_receptor, niche, prioritization_score) %>% group_by(ligand_receptor, niche) %>% top_n(1, prioritization_score) %>% distinct() %>% mutate() %>% spread(niche, prioritization_score, fill = 0)

lr_KC_niche = prioritization_tbl_ligand_receptor_spread %>% filter(KC > MoMac1 &  KC > MoMac1_CV & KC > MoMac2) %>% pull(ligand_receptor) %>% unique()
lr_MoMac2_niche = prioritization_tbl_ligand_receptor_spread %>% filter( MoMac2 > KC) %>% pull(ligand_receptor) %>% unique()
lr_MoMac1_niche = prioritization_tbl_ligand_receptor_spread %>% filter( MoMac1 > KC) %>% pull(ligand_receptor) %>% unique()
lr_MoMac1_CV_niche = prioritization_tbl_ligand_receptor_spread %>% filter( MoMac1_CV > KC) %>% pull(ligand_receptor) %>% unique()

### Ligand-Receptor pair: only keep top2 receptor per ligand  -----------------------

receiver_oi = "KCs"

prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

filtered_ligands = prioritized_tbl_oi %>% pull(ligand) %>% unique()

prioritized_tbl_lr_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(1, prioritization_score) %>% ungroup()
prioritized_tbl_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup()  %>% filter(prioritization_score > prioritized_tbl_lr_oi$prioritization_score %>% min())
prioritized_tbl_oi = prioritized_tbl_oi %>% filter(ligand_receptor %in% lr_KC_niche)
lfc_plot = make_ligand_receptor_lfc_zonation_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
ggsave("plots/KC_MoMac1_MoMac2/p_LR_pair_lfc_KCs.pdf",width = 15, height = 25)

receiver_oi = "MoMac2"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

filtered_ligands = prioritized_tbl_oi %>% pull(ligand) %>% unique()

prioritized_tbl_lr_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(1, prioritization_score) %>% ungroup()
prioritized_tbl_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup()  %>% filter(prioritization_score > prioritized_tbl_lr_oi$prioritization_score %>% min())
prioritized_tbl_oi = prioritized_tbl_oi %>% filter(ligand_receptor %in% lr_MoMac2_niche)

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
ggsave("plots/KC_MoMac1_MoMac2/p_LR_pair_lfc_MoMac2.pdf",width = 10, height = 22)

niche_oi = "MoMac1"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

filtered_ligands = prioritized_tbl_oi %>% pull(ligand) %>% unique()

prioritized_tbl_lr_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(niche == niche_oi) %>% top_n(1, prioritization_score) %>% ungroup()
prioritized_tbl_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(niche == niche_oi) %>% top_n(2, prioritization_score) %>% ungroup()  %>% filter(prioritization_score > prioritized_tbl_lr_oi$prioritization_score %>% min())
prioritized_tbl_oi = prioritized_tbl_oi %>% filter(ligand_receptor %in% lr_MoMac1_niche)

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
ggsave("plots/KC_MoMac1_MoMac2/p_LR_pair_lfc_MoMac1.pdf",width = 10, height = 22)

niche_oi = "MoMac1_CV"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

filtered_ligands = prioritized_tbl_oi %>% pull(ligand) %>% unique()

prioritized_tbl_lr_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(niche == niche_oi) %>% top_n(1, prioritization_score) %>% ungroup()
prioritized_tbl_oi = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(niche == niche_oi) %>% top_n(2, prioritization_score) %>% ungroup()  %>% filter(prioritization_score > prioritized_tbl_lr_oi$prioritization_score %>% min())
prioritized_tbl_oi = prioritized_tbl_oi %>% filter(ligand_receptor %in% lr_MoMac1_CV_niche)

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
ggsave("plots/KC_MoMac1_MoMac2/p_LR_pair_lfc_MoMac1_CV.pdf",width = 10, height = 22)


####################################### ####################################### #######################################
####################################### Combine Expression tables for visualization ####################################### ###############################
####################################### ####################################### #######################################

### Ligand-Exprs-Activity-Targets -----------------------

### remark: I will probably need to make an exprs table of all genes over all cell types I have here
### but maybe I can already give a shortlist of genes needed!!
### so it wont be all genes...

no_exprs_table = FALSE
if(no_exprs_table == TRUE){
  celltypes = union(prioritization_tbl_ligand_target %>% pull(sender) %>% unique(), prioritization_tbl_ligand_target %>% pull(receiver) %>% unique())
  genes = union(prioritization_tbl_ligand_target %>% pull(target) %>% unique(), prioritization_tbl_ligand_target %>% pull(ligand) %>% unique()) %>% union(prioritization_tbl_ligand_receptor %>% pull(ligand) %>% unique()) %>% union(prioritization_tbl_ligand_receptor %>% pull(receptor) %>% unique())

  celltypes %>% saveRDS("output/celltypes_oi_exprs_table_mouse.rds")
  genes %>% saveRDS("output/genes_oi_exprs_table_mouse.rds")

  source("scripts/get_exprs_table_mouse_prism.R")

} else {
  exprs_table_mouse = readRDS("output/exprs_table_mouse.rds")

}

exprs_table_mouse = exprs_table_mouse %>% mutate(celltype = as.character(celltype)) %>%
  mutate_cond(celltype == "Stellate cells_portal", celltype = "Stellate cells") %>%
  mutate_cond(celltype == "Hepatocytes_portal", celltype = "Hepatocytes") %>%
  mutate_cond(celltype == "LSECs_portal", celltype = "LSECs")

exprs_tbl_ligand = exprs_table_mouse %>% filter(gene %in% prioritization_tbl_ligand_receptor$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
exprs_tbl_receptor = exprs_table_mouse %>% filter(gene %in% prioritization_tbl_ligand_receptor$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_table_mouse %>% filter(gene %in% prioritization_tbl_ligand_target$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

exprs_tbl_ligand = exprs_tbl_ligand %>% left_join(prioritization_tbl_ligand_receptor %>% distinct(niche, sender))
exprs_tbl_receptor = exprs_tbl_receptor %>% left_join(prioritization_tbl_ligand_receptor %>% distinct(niche, receiver))
exprs_tbl_target = exprs_tbl_target %>% left_join(prioritization_tbl_ligand_receptor %>% distinct(niche, receiver))


receiver_oi = "KCs"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

exprs_plot = make_ligand_zonation_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_plot$combined_plot
ggsave("plots/KC_MoMac1_MoMac2/p_ligand_exprs_KCs.pdf",width = 35, height = 16)

receiver_oi = "MoMac2"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

exprs_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_plot$combined_plot
ggsave("plots/KC_MoMac1_MoMac2/p_ligand_exprs_MoMac2.pdf",width = 35, height = 15)

niche_oi = "MoMac1"
receiver_oi = "MoMac1"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

exprs_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_plot$combined_plot
ggsave("plots/KC_MoMac1_MoMac2/p_ligand_exprs_MoMac1.pdf",width = 20, height = 15)

niche_oi = "MoMac1_CV"
receiver_oi = "MoMac1"
prioritized_tbl_oi = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% group_by(sender) %>% top_n(20, prioritization_score)
prioritized_tbl_oi2 = prioritization_tbl_ligand_receptor %>% filter(ligand_score > 0 & ligand_significant >= 0.50) %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(niche == niche_oi) %>% ungroup() %>% top_n(70, prioritization_score)
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

exprs_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_plot$combined_plot
ggsave("plots/KC_MoMac1_MoMac2/p_ligand_exprs_MoMac1_CV.pdf",width = 20, height = 15)


### PRIORITIZIATION SCORE HEATMAP -----------------------

# prioritization score per ligand-receptor pair and then heatmap

prioritization_tbl_ligand_receptor = prioritization_tbl_ligand_receptor %>%
  mutate_cond(niche == "MoMac1", niche = "MoMac1-Cap") %>%
  mutate_cond(niche == "MoMac1_CV", niche = "MoMac1-CV")

prioritization_tbl_ligand_receptor_spread = prioritization_tbl_ligand_receptor %>% select(ligand, receptor, ligand_receptor, niche, prioritization_score) %>% group_by(ligand_receptor, niche) %>% top_n(1, prioritization_score) %>% distinct()  %>% spread(niche, prioritization_score, fill = 0)

prioritization_score_cutoff = prioritization_tbl_ligand_receptor_spread %>% filter(KC > `MoMac1-Cap` &  KC > `MoMac1-CV` & KC > MoMac2) %>% group_by(ligand) %>% top_n(1, KC) %>% ungroup() %>% top_n(50, KC) %>% pull(KC) %>% min()

lr_KC_niche = prioritization_tbl_ligand_receptor_spread %>% filter(KC > `MoMac1-Cap` &  KC > `MoMac1-CV` & KC > MoMac2) %>% group_by(ligand) %>% top_n(1, KC) %>% ungroup() %>% filter(KC > prioritization_score_cutoff) %>% pull(ligand_receptor) %>% unique()
lr_MoMac2_niche =  prioritization_tbl_ligand_receptor_spread %>% filter( MoMac2 > KC) %>% group_by(ligand) %>% top_n(1, MoMac2) %>% ungroup() %>% filter(MoMac2 > prioritization_score_cutoff)  %>% pull(ligand_receptor) %>% unique()
lr_MoMac1_niche =  prioritization_tbl_ligand_receptor_spread %>% filter( `MoMac1-Cap` > KC)  %>% group_by(ligand) %>% top_n(1, `MoMac1-Cap`) %>% ungroup() %>% filter(`MoMac1-Cap` > prioritization_score_cutoff)  %>% pull(ligand_receptor) %>% unique()
lr_MoMac1_CV_niche =  prioritization_tbl_ligand_receptor_spread %>% filter( `MoMac1-CV` > KC) %>% group_by(ligand) %>% top_n(1, `MoMac1-CV`) %>% ungroup() %>% filter(`MoMac1-CV` > prioritization_score_cutoff)  %>% pull(ligand_receptor) %>% unique()

prioritization_tbl_ligand_receptor_spread_filtered = prioritization_tbl_ligand_receptor_spread %>% filter(ligand_receptor %in% c(lr_KC_niche, lr_MoMac2_niche, lr_MoMac1_niche, lr_MoMac1_CV_niche)) %>% select(-ligand, -receptor) %>% distinct()

prioritization_tbl_ligand_receptor_spread_filtered_mat = prioritization_tbl_ligand_receptor_spread_filtered %>% ungroup() %>% select(-ligand_receptor) %>% data.frame() %>% as.matrix() %>% magrittr::set_rownames(prioritization_tbl_ligand_receptor_spread_filtered$ligand_receptor)

dist_lr = dist(1-cor(t(prioritization_tbl_ligand_receptor_spread_filtered_mat)))
hclust_lr = hclust(dist_lr, method = "ward.D2")

dist_receiver = dist(1-cor(prioritization_tbl_ligand_receptor_spread_filtered_mat))
hclust_receiver = hclust(dist_receiver, method = "ward.D2")
pheatmap::pheatmap(prioritization_tbl_ligand_receptor_spread_filtered_mat, border_color = "white", cluster_rows = hclust_lr, cluster_cols = hclust_receiver)

dev.off()
pdf("plots/KC_MoMac1_MoMac2/prioritization_score_heatmap.pdf", height = 17, width = 4)
pheatmap::pheatmap(prioritization_tbl_ligand_receptor_spread_filtered_mat, border_color =  "white", cluster_rows = hclust_lr, cluster_cols = hclust_receiver)
dev.off()


order_ligands_receptor = hclust_lr$labels[hclust_lr$order]
order_receivers = hclust_receiver$labels[hclust_receiver$order]

order_receivers[order_receivers ==  "MoMac1.CV" ] = "MoMac1-CV"
order_receivers[order_receivers ==  "MoMac1.Cap" ] = "MoMac1-Cap"

### do it via ggplot

prioritization_tbl_ligand_receptor_selection_sender_niche = prioritization_tbl_ligand_receptor %>% filter(ligand_receptor %in% prioritization_tbl_ligand_receptor_spread_filtered$ligand_receptor) %>% group_by(ligand_receptor) %>% top_n(1, prioritization_score) %>% select(sender, ligand_receptor, niche) %>% distinct()
prioritization_tbl_ligand_receptor_selection = prioritization_tbl_ligand_receptor %>% select(ligand, receptor, ligand_receptor, niche, prioritization_score, sender) %>% group_by(ligand_receptor, niche) %>% top_n(1, prioritization_score) %>% distinct() %>% filter(ligand_receptor %in% prioritization_tbl_ligand_receptor_spread_filtered$ligand_receptor)
prioritization_tbl_ligand_receptor_selection = prioritization_tbl_ligand_receptor_selection %>% inner_join(prioritization_tbl_ligand_receptor_selection_sender_niche %>% select(-niche) %>% rename(sender_top = sender))

prioritization_tbl_ligand_receptor_selection = prioritization_tbl_ligand_receptor_selection %>%
  mutate_cond(sender_top == "Central Vein Endothelial cells", sender_top = "Central Vein\nEndothelial cells") %>%
  mutate_cond(sender_top == "Capsule fibroblasts", sender_top = "Capsule\nfibroblasts") %>%
  mutate_cond(sender_top == "Mesothelial cells", sender_top = "Mesothelial\ncells") %>%
  mutate_cond(sender_top == "Stellate cells", sender_top = "Stellate\ncells")

levels_senders = c("Stellate\ncells", "LSECs", "Hepatocytes","Fibroblast 1", "Central Vein\nEndothelial cells", "Capsule\nfibroblasts", "Mesothelial\ncells", "Fibroblast 2", "Cholangiocytes")
prioritization_tbl_ligand_receptor_selection = prioritization_tbl_ligand_receptor_selection %>% mutate(niche = factor(niche, levels = order_receivers), ligand_receptor = factor(ligand_receptor, levels = order_ligands_receptor), sender_top = factor(sender_top, levels = levels_senders))

p1 = prioritization_tbl_ligand_receptor_selection %>%
  ggplot(aes(niche , ligand_receptor , fill = prioritization_score )) +
  geom_tile(color = "white") +
  facet_grid(sender_top~., scales = "free", space = "free") +
  scale_x_discrete(position = "top") +
  theme_light() +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(face = "bold.italic", size = 9),
    axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "bold"),
    strip.text.x.top = element_text(angle = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.25, "lines"),
    panel.spacing.y = unit(0.25, "lines"),
    strip.text.x = element_text(size = 7, color = "black", face = "bold"),
    strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = 7, color = "black", face = "bold")
  ) + labs(fill = "Differential NicheNet\nPrioritization Score")
max_score = prioritization_tbl_ligand_receptor_selection$prioritization_score %>% max() + 0.01
min_score = prioritization_tbl_ligand_receptor_selection$prioritization_score %>% min() - 0.01
custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.30, 0.40, 0.5, 0.70, 0.75, 1),  limits = c(min_score, max_score))

p_prioritization = p1 + custom_scale_fill + ylab("Prioritized Ligand-Receptor Pairs") + xlab("")
p_prioritization
ggsave("plots/KC_MoMac1_MoMac2/p_prioritization.pdf",width = 4.5, height = 16)

p1 = prioritization_tbl_ligand_receptor_selection %>%
  ggplot(aes(niche , ligand_receptor , fill = prioritization_score )) +
  geom_tile(color = "white") +
  facet_grid(sender_top~., scales = "free", space = "free") +
  scale_x_discrete(position = "top") +
  theme_light() +
  theme(
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(face = "bold.italic", size = 9),
    axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "bold"),
    strip.text.x.top = element_text(angle = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.25, "lines"),
    panel.spacing.y = unit(0.25, "lines"),
    strip.text.x = element_text(size = 7, color = "black", face = "bold"),
    strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = 8, color = "black", face = "bold", angle = 0)
  ) + labs(fill = "Differential NicheNet\nPrioritization Score")
p_prioritization = p1 + custom_scale_fill + ylab("Prioritized Ligand-Receptor Pairs") + xlab("")
p_prioritization
ggsave("plots/KC_MoMac1_MoMac2/p_prioritization_v2.pdf",width = 5.15, height = 16)
