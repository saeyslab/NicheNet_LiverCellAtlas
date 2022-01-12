library(nichenetr)
library(tidyverse)

# LR network
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# Ligand-target matrix
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

# human signature
humanKC = xlsx::read.xlsx2("data/humanKC_markers.xlsx", sheetName = "Sheet1") %>% pull(gene) %>% unique()
humanKC_signature = humanKC %>% intersect(rownames(ligand_target_matrix))
humanKC_signature %>% setdiff(rownames(ligand_target_matrix))

# mouse signature
mouseKC = xlsx::read.xlsx2("data/mouseKC_markers.xlsx", sheetName = "Sheet1") %>% pull(gene) %>% unique()
mouseKC_signature = mouseKC %>% intersect(rownames(ligand_target_matrix) %>% convert_human_to_mouse_symbols())
mouseKC_signature %>% setdiff(rownames(ligand_target_matrix) %>% convert_human_to_mouse_symbols())

mouseKC_signature_humansymbols = mouseKC_signature %>% convert_mouse_to_human_symbols() %>% .[!is.na(.)]

# human KC specific genes
humanKC_signature %>% setdiff(mouseKC_signature_humansymbols) %>% saveRDS("output/humanKC_specific_genes.rds")

# define ligands per KC niche
lr_network %>% filter(to %in% c(humanKC_signature %>% setdiff(mouseKC_signature_humansymbols))) %>% pull(from) %>% unique() %>% intersect(readRDS("output/KC_niche_ligands.rds"))
lr_network %>% filter(!database %in% c("ppi_prediction","ppi_prediction_go")) %>% filter(to %in% c(humanKC_signature %>% setdiff(mouseKC_signature_humansymbols))) %>% pull(from) %>% unique() %>% intersect(readRDS("output/KC_niche_ligands.rds"))

# ligand activity analysis: -- with all ligands in the database as potential ligand -- mouseKC or humanKC specific genes vs background (genomic or KC only)
# mouse vs bg
background_expressed_genes = rownames(ligand_target_matrix)
potential_ligands = colnames(ligand_target_matrix)
ligand_activities_mouse_vs_bg = predict_ligand_activities(geneset = mouseKC_signature_humansymbols, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_mouse_vs_bg %>% arrange(-pearson)

# human vs bg
ligand_activities_human_vs_bg = predict_ligand_activities(geneset = humanKC_signature, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_human_vs_bg %>% arrange(-pearson)

# mouse-human-DE vs bg
ligand_activities_mousehuman_vs_bg = predict_ligand_activities(geneset = union(humanKC_signature, mouseKC_signature_humansymbols), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_mousehuman_vs_bg %>% arrange(-pearson)

# mouse-human DE mouse-specific vg bg
ligand_activities_mousespecific_vs_bg = predict_ligand_activities(geneset = mouseKC_signature_humansymbols %>% setdiff(humanKC_signature), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_mousespecific_vs_bg %>% arrange(-pearson)

# mouse-human DE human-specific vg bg
ligand_activities_humanspecific_vs_bg = predict_ligand_activities(geneset = humanKC_signature %>% setdiff(mouseKC_signature_humansymbols), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_humanspecific_vs_bg %>% arrange(-pearson)

# mouse-specific vs mouse+human
ligand_activities_mousespecific_vs_human_specific = predict_ligand_activities(geneset = mouseKC_signature_humansymbols %>% setdiff(humanKC_signature), background_expressed_genes = union(humanKC_signature, mouseKC_signature_humansymbols), ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_mousespecific_vs_human_specific %>% arrange(-pearson)

# human-specific vs mouse+human
ligand_activities_humanspecific_vs_mousespecific = predict_ligand_activities(geneset = humanKC_signature %>% setdiff(mouseKC_signature_humansymbols), background_expressed_genes = union(humanKC_signature, mouseKC_signature_humansymbols), ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_humanspecific_vs_mousespecific %>% arrange(-pearson)

##### Combine all ligand activities
ligand_activities_df = list(ligand_activities_mouse_vs_bg %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(mouse_vs_bg = scaling_zscore(pearson)) %>% select(-pearson), #################################
                            ligand_activities_human_vs_bg %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(human_vs_bg = scaling_zscore(pearson)) %>% select(-pearson),
                            ligand_activities_mousehuman_vs_bg %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(mousehuman_vs_bg = scaling_zscore(pearson)) %>% select(-pearson),
                            ligand_activities_mousespecific_vs_bg %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(mouse_specific_vs_bg = scaling_zscore(pearson)) %>% select(-pearson),
                            ligand_activities_humanspecific_vs_bg %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(human_specific_vs_bg = scaling_zscore(pearson)) %>% select(-pearson),
                            ligand_activities_mousespecific_vs_human_specific %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(mouse_specific_vs_mousehumanKC = scaling_zscore(pearson)) %>% select(-pearson) ,
                            ligand_activities_humanspecific_vs_mousespecific %>% drop_na() %>% distinct(test_ligand, pearson) %>% mutate(human_specific_vs_mousehumanKC = scaling_zscore(pearson))%>% select(-pearson)) %>% reduce(inner_join)

ligand_activities_df %>% saveRDS("output/ligand_activities_humanKC_vs_mouseKC.rds")

