library(nichenetr)
library(tidyverse)
library(circlize)

source("scripts/plotting_functions.R")
source("scripts/functions_differential_nichenet.R")
make_circos_lr= function(prioritized_tbl_oi, colors_sender, colors_receiver, cutoff, scale, transparency = NULL, circos_type, border = TRUE){

  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")

  # Link each cell type to a color
  grid_col_ligand = colors_sender
  # names(grid_col_ligand) = prioritized_tbl_oi$sender %>% unique() %>% sort()

  grid_col_receptor = colors_receiver
  # names(grid_col_receptor) = prioritized_tbl_oi$receiver %>% unique() %>% sort()

  grid_col_tbl_ligand = tibble::tibble(sender = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_receptor = tibble::tibble(receiver = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

  # Make the plot for condition of interest - title of the plot
  circos_links_oi = prioritized_tbl_oi

  # deal with duplicated sector names
  # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
  # only do it with duplicated ones!
  circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)

  #circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score) %>% dplyr::mutate(ligand = paste(sender, ligand, sep = "_"), receptor = paste(receptor, receiver, sep = "_"))

  df = circos_links %>% mutate(ligand_receptor_sender_receiver = paste0(sender, receiver, ligand_receptor))

  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$ligand_receptor_sender_receiver %in% df.i.j$ligand_receptor_sender_receiver] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$ligand_receptor_sender_receiver %in% df.i.j$ligand_receptor_sender_receiver] = df.i.j$receptor
    }
  }

  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

  # print(intersecting_ligands_receptors)

  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(receptor, " ", sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }

  circos_links = df

  # Link ligands/Receptors to the colors of senders/receivers
  circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
  ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)

  # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
  # transparency = circos_links %>% mutate(old_weight = weight) %>% mutate_cond(old_weight < cutoff, weight = 0) %>% mutate_cond(old_weight >= cutoff, weight = 2)  %>% mutate_cond(old_weight > cutoff+0.15, weight = 3) %>% mutate_cond(old_weight > cutoff+0.2, weight = 3.5)  %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
  if(is.null(transparency)) {
    # transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    # transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    transparency = circos_links %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

  }

  # Define order of the ligands and receptors and the gaps
  ligand_order_automatic = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    # ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(weight) %>% dplyr::distinct(ligand)

  }) %>% unlist()


  # hepatocyte_ligands = c( "C5", "APOA1","F9", "F2", "APOB", "FGA", "TF","PLG","TTR", "ANGPTL3", "COL18A1", "COL5A3", "SAA1",  "SERPINA1", "SERPINC1", "KNG1","HAMP")  %>% intersect(prioritized_tbl_oi %>% filter(sender =="Hepatocytes") %>% pull(ligand))
  # lsec_ligands = c("DLL4", "ITGA9", "CDH5", "F8","ANGPT2", "CD274","BMP6","BMP2", "HSP90B1", "CSF1") %>% intersect(prioritized_tbl_oi %>% filter(sender =="LSECs") %>% pull(ligand))
  # stellate_ligands = c("CDH2", "ITGA9","IL34","NGF", "SEMA6D", "RELN",  "LAMB1", "COL4A1", "JAM3", "GDF2","BMP10","BMP5", "ITGB1", "HGF")  %>% intersect(prioritized_tbl_oi %>% filter(sender =="Stellate cells") %>% pull(ligand))
  #
  # # hepatocyte_ligands = c("ANGPTL3", "C3", "C5","F9","F2", "FGA","CFH","TF","TTR","HEBP1","Hebp1","APOA1","APOB","SERPINA1","SERPINC1","VTN","Vtn","SAA1", "PLG","KNG1","A2M","COL18A1", "COL5A3","HAMP")  %>% intersect(prioritized_tbl_oi %>% filter(sender =="Hepatocytes") %>% pull(ligand))
  # # lsec_ligands = c("DLL4","BMP2","BMP6","CALM3","Calm2", "CCL14","CCL23","CXCL10", "Cxcl10","CXCL9", "Cxcl9","CD274", "EFNB1", "Efnb1", "EFNB2","Efnb2","F8","TFPI","Tfpi", "HSP90B1", "ADAM23","Adam23","ANGPT2","RELN", "ICAM1","Icam1", "ITGA9","CDH5","CSF1") %>% intersect(prioritized_tbl_oi %>% filter(sender =="LSECs") %>% pull(ligand))
  # # stellate_ligands = c("IL34","IL7","GDF2","BMP10","BMP5","VEGFC", "HGF","FGF14", "TNFRSF11B","Tnfrsf11b", "MAPT","Mapt","NGF","NRXN1", "Nrxn1" ,"NTF3","PGF","PROS1","SEMA5A", "SEMA6D" ,"RELN","ITGA4","ITGA9","ITGB1","JAM3","LAMA1","Lama1","LAMA2","LAMB1","LAMC3","PCDH9","CDH2","COL4A1")  %>% intersect(prioritized_tbl_oi %>% filter(sender =="Stellate cells") %>% pull(ligand))
  #
  # if("RELN" %in% lsec_ligands){
  #   stellate_ligands[stellate_ligands == "RELN"] = "RELN "
  # }
  # if("ITGA9" %in% lsec_ligands){
  #   stellate_ligands[stellate_ligands == "ITGA9"] = "ITGA9 "
  # }
  # "RELN ", "ITGA9 ": both stellate

  # hepatocyte_ligands = hepatocyte_ligands %>% intersect(prioritized_tbl_oi %>% filter(sender =="Hepatocytes") %>% pull(ligand))
  # lsec_ligands = lsec_ligands %>% intersect(prioritized_tbl_oi %>% filter(sender =="LSECs") %>% pull(ligand))
  # stellate_ligands = stellate_ligands %>% intersect(prioritized_tbl_oi %>% filter(sender =="Stellate cells") %>% pull(ligand))


  # ligand_order = c(hepatocyte_ligands, lsec_ligands, stellate_ligands)
  # ligand_order = ligand_order %>% intersect(ligand_order_automatic)
  ligand_order = ligand_order_automatic
  receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    circos_links_n = circos_links %>% dplyr::filter(receiver == receiver_oi) %>% group_by(receptor) %>% count() %>% ungroup()
    receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>% inner_join(circos_links_n) %>% dplyr::arrange(-n, ligand) %>% dplyr::distinct(receptor)
    # receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(weight) %>% dplyr::distinct(receptor)

  }) %>% unlist()

  # receptor_order_last = c("CSF1R","NOTCH2")
  # receptor_order_first = c("SLC40A1","BMPR1A","BMPR2","ACVRL1","LRP6")
  # receptor_order = c(receptor_order_first, setdiff(receptor_order, c(receptor_order_last,receptor_order_first)), receptor_order_last)

  order = c(ligand_order,receptor_order)
  # print(length(order))

  width_same_cell_same_ligand_type = 1
  width_different_cell = 5
  width_ligand_receptor = 15
  width_same_cell_same_receptor_type = 1

  sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]

  receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]

  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)

  # print(length(gaps))
  # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
  if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }



  # links_circle$weight[links_circle$weight == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)

  if(circos_type == "arrow"){
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = FALSE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight >= cutoff,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.25),
                 # grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                 reduce = 0
                 , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
    )
  } else {
    if(border == TRUE) {
      chordDiagram(links_circle,
                   directional = 1,
                   order=order,
                   link.sort = TRUE,
                   link.decreasing = FALSE,
                   grid.col = grid_col,
                   transparency = transparency,
                   diffHeight = 0.0075,
                   direction.type = c("diffHeight", "arrows"),
                   link.visible = links_circle$weight >= cutoff,
                   annotationTrack = "grid",
                   preAllocateTracks = list(track.height = 0.25),
                   grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                   # link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                   reduce = 0
                   , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
      )
    } else {
      chordDiagram(links_circle,
                   directional = 1,
                   order=order,
                   link.sort = TRUE,
                   link.decreasing = FALSE,
                   grid.col = grid_col,
                   transparency = transparency,
                   diffHeight = 0.0075,
                   direction.type = c("diffHeight", "arrows"),
                   link.visible = links_circle$weight >= cutoff,
                   annotationTrack = "grid",
                   preAllocateTracks = list(track.height = 0.25),
                   link.arr.length = 0.05, link.arr.type = "big.arrow",
                   # link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                   reduce = 0
                   , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
      )
    }

  }




  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  }, bg.border = NA) #

  p_circos = recordPlot()

  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(grid_col_receptor, grid_col_ligand)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_receptor[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))

  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_ligand[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))

  p_legend = grDevices::recordPlot()


  return(list(p_circos = p_circos, p_legend = p_legend))

}


####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################
############## DIFFERENTIAL NICHENET RESULTS  ############## ############## ############## ############## ##############
####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################


####################################### ####################################### #######################################
####################################### MOUSE ####################################### ###############################
####################################### ####################################### #######################################

output_nichenet_analysis_mouse = readRDS("output/differential_nichenet_mouse_KC_MF_new.rds")

expression_pct = 0.10

# classic analysis
# ---------- Combine all the information for prioritization------------------------------------------------------------

combined_information = output_nichenet_analysis_mouse$DE_sender_receiver %>%
  inner_join(output_nichenet_analysis_mouse$ligand_scaled_receptor_expression_fraction_df, by = c("receiver", "ligand", "receptor")) %>%
  inner_join(output_nichenet_analysis_mouse$sender_spatial_DE_processed, by = c("niche", "sender", "ligand")) %>%
  inner_join(output_nichenet_analysis_mouse$receiver_spatial_DE_processed, by = c("niche", "receiver", "receptor")) %>%
  inner_join(output_nichenet_analysis_mouse$ligand_activities_targets, by = c("receiver", "ligand")) %>%
  inner_join(output_nichenet_analysis_mouse$DE_receiver_processed_targets, by = c("niche", "receiver", "target")) %>%
  inner_join(output_nichenet_analysis_mouse$exprs_tbl_ligand, by = c("sender", "ligand")) %>%
  inner_join(output_nichenet_analysis_mouse$exprs_tbl_receptor, by = c("receiver", "receptor")) %>%
  inner_join(output_nichenet_analysis_mouse$exprs_tbl_target, by = c("receiver", "target"))

# reorder the columns

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

combined_information = combined_information %>% mutate(ligand_receptor = paste(ligand, receptor, sep = "--"))  %>%  mutate(bonafide_score = 1) %>%  mutate_cond(bonafide == FALSE, bonafide_score = 0.5)

combined_information = combined_information %>% select(
  niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
  ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_zonation,
  receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_zonation,
  ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor, bonafide_score,
  target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
  activity, activity_normalized,
  scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
  scaled_ligand_score_zonation, scaled_receptor_score_zonation,
  scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
  scaled_activity, scaled_activity_normalized
) %>% distinct()


# combine in a prioritization score

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

combined_information_prioritized = combined_information %>%
  dplyr::mutate(prioritization_score =
                  ((prioritizing_weights["scaled_ligand_score"] * scaled_ligand_score) +
                     (prioritizing_weights["scaled_ligand_expression_scaled"] * scaled_ligand_expression_scaled) +
                     (prioritizing_weights["scaled_receptor_score"] * scaled_receptor_score) +
                     (prioritizing_weights["scaled_receptor_expression_scaled"] * scaled_receptor_expression_scaled) +
                     (prioritizing_weights["scaled_avg_score_ligand_receptor"] * scaled_avg_score_ligand_receptor) +
                     (prioritizing_weights["scaled_ligand_score_zonation"] * scaled_ligand_score_zonation) +
                     (prioritizing_weights["scaled_receptor_score_zonation"] * scaled_receptor_score_zonation) +
                     (prioritizing_weights["ligand_scaled_receptor_expression_fraction"] * ligand_scaled_receptor_expression_fraction) +
                     (prioritizing_weights["scaled_activity"] * scaled_activity) +
                     (prioritizing_weights["scaled_activity_normalized"] * scaled_activity_normalized) +
                     (prioritizing_weights["ligand_fraction"] * scaled_ligand_fraction_adapted ) +
                     (prioritizing_weights["receptor_fraction"] * scaled_receptor_fraction_adapted  ) +
                     (prioritizing_weights["bona_fide"] * bonafide_score)
                  )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)

prioritization_tbl_ligand_receptor_mouse = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide,
                                                                                       ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_zonation,
                                                                                       receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_zonation,
                                                                                       ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor,
                                                                                       activity, activity_normalized,
                                                                                       scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
                                                                                       scaled_ligand_score_zonation, scaled_receptor_score_zonation,
                                                                                       scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
                                                                                       scaled_activity, scaled_activity_normalized,
                                                                                       prioritization_score) %>% distinct()
prioritization_tbl_ligand_target_mouse = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
                                                                                     target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
                                                                                     activity, activity_normalized, scaled_activity, scaled_activity_normalized, prioritization_score) %>% distinct()


####################################### ####################################### #######################################
####################################### HUMAN ####################################### ###############################
####################################### ####################################### #######################################

output_nichenet_analysis_human = readRDS("output/differential_nichenet_human_KC_MF_noPV_EC.rds")

expression_pct = 0.10

# classic analysis
# ---------- Combine all the information for prioritization------------------------------------------------------------

combined_information = output_nichenet_analysis_human$DE_sender_receiver %>%
  inner_join(output_nichenet_analysis_human$ligand_scaled_receptor_expression_fraction_df, by = c("receiver", "ligand", "receptor")) %>%
  # inner_join(output_nichenet_analysis_human$sender_spatial_DE_processed, by = c("niche", "sender", "ligand")) %>%
  # inner_join(output_nichenet_analysis_human$receiver_spatial_DE_processed, by = c("niche", "receiver", "receptor")) %>%
  inner_join(output_nichenet_analysis_human$ligand_activities_targets, by = c("receiver", "ligand")) %>%
  inner_join(output_nichenet_analysis_human$DE_receiver_processed_targets, by = c("niche", "receiver", "target")) %>%
  inner_join(output_nichenet_analysis_human$exprs_tbl_ligand, by = c("sender", "ligand")) %>%
  inner_join(output_nichenet_analysis_human$exprs_tbl_receptor, by = c("receiver", "receptor")) %>%
  inner_join(output_nichenet_analysis_human$exprs_tbl_target, by = c("receiver", "target"))

# reorder the columns

combined_information = combined_information %>% mutate(ligand_receptor = paste(ligand, receptor, sep = "--"))  %>%  mutate(bonafide_score = 1) %>%  mutate_cond(bonafide == FALSE, bonafide_score = 0.5)

# # in the future this:
combined_information = combined_information %>% select(
  niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
  ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction,
  receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction,
  ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor, bonafide_score,
  target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
  activity, activity_normalized,
  scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
  scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
  scaled_activity, scaled_activity_normalized
) %>% distinct()

# combine in a prioritization score

prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 2,
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "scaled_avg_score_ligand_receptor" = 0,
                         # "scaled_ligand_score_zonation" = 2,
                         # "scaled_receptor_score_zonation" = 0,
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_activity" = 0,"scaled_activity_normalized" = 1,
                         "ligand_fraction" = 1, "receptor_fraction" = 1,
                         "bona_fide" = 1)

combined_information_prioritized = combined_information %>%
  dplyr::mutate(prioritization_score =
                  ((prioritizing_weights["scaled_ligand_score"] * scaled_ligand_score) +
                     (prioritizing_weights["scaled_ligand_expression_scaled"] * scaled_ligand_expression_scaled) +
                     (prioritizing_weights["scaled_receptor_score"] * scaled_receptor_score) +
                     (prioritizing_weights["scaled_receptor_expression_scaled"] * scaled_receptor_expression_scaled) +
                     (prioritizing_weights["scaled_avg_score_ligand_receptor"] * scaled_avg_score_ligand_receptor) +
                     # (prioritizing_weights["scaled_ligand_score_zonation"] * scaled_ligand_score_zonation) +
                     # (prioritizing_weights["scaled_receptor_score_zonation"] * scaled_receptor_score_zonation) +
                     (prioritizing_weights["ligand_scaled_receptor_expression_fraction"] * ligand_scaled_receptor_expression_fraction) +
                     (prioritizing_weights["scaled_activity"] * scaled_activity) +
                     (prioritizing_weights["scaled_activity_normalized"] * scaled_activity_normalized) +
                     (prioritizing_weights["ligand_fraction"] * scaled_ligand_fraction_adapted ) +
                     (prioritizing_weights["receptor_fraction"] * scaled_receptor_fraction_adapted  ) +
                     (prioritizing_weights["bona_fide"] * bonafide_score)
                  )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)

prioritization_tbl_ligand_receptor_human = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide,
                                                                                       ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction,
                                                                                       receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction,
                                                                                       ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor,
                                                                                       activity, activity_normalized,
                                                                                       scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
                                                                                       scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
                                                                                       scaled_activity, scaled_activity_normalized,
                                                                                       prioritization_score) %>% distinct()
prioritization_tbl_ligand_target_human = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
                                                                                     target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
                                                                                     activity, activity_normalized, scaled_activity, scaled_activity_normalized, prioritization_score) %>% distinct()

####################################### ####################################### #######################################
####################################### Compare human and mouse ####################################### ###############################
####################################### ####################################### #######################################

# We are mainly interested in the conservation of the KC niche between human and mouse
# because not exactly the same sender cell types
# still need to make the heatmaps showing the output of both human and mouse separately

###### RECEPTOR LEVEL ############

prioritization_tbl_ligand_receptor_mouse = prioritization_tbl_ligand_receptor_mouse %>%
  mutate_cond(sender == "Stellate cells_portal", sender = "Stellate cells") %>%
  mutate_cond(sender == "Hepatocytes_portal", sender = "Hepatocytes") %>%
  mutate_cond(sender == "LSECs_portal", sender = "LSECs")

prioritization_tbl_ligand_receptor_human = prioritization_tbl_ligand_receptor_human %>%
  mutate_cond(receiver  == "resKCs", receiver  = "KCs")


prioritization_tbl_ligand_receptor_human = prioritization_tbl_ligand_receptor_human %>% mutate(prioritization_score = scale_quantile_adapted(prioritization_score))
prioritization_tbl_ligand_receptor_mouse = prioritization_tbl_ligand_receptor_mouse %>% mutate(prioritization_score = scale_quantile_adapted(prioritization_score))


prioritization_tbl_ligand_receptor_human_orthologs = prioritization_tbl_ligand_receptor_human  %>% select(sender, receiver, ligand, receptor, prioritization_score) %>% rename(human_score = prioritization_score) %>% mutate(mouse_ligand = nichenetr::convert_human_to_mouse_symbols(ligand), mouse_receptor = nichenetr::convert_human_to_mouse_symbols(receptor)) %>% distinct()
prioritization_tbl_ligand_receptor_mouse_orthologs = prioritization_tbl_ligand_receptor_mouse  %>% select(sender, receiver, ligand, receptor, prioritization_score) %>% rename(mouse_score = prioritization_score) %>% rename(mouse_ligand = ligand, mouse_receptor = receptor) %>% mutate(ligand = nichenetr::convert_mouse_to_human_symbols(mouse_ligand), receptor = nichenetr::convert_mouse_to_human_symbols(mouse_receptor)) %>% distinct()


## common inner join on both ligand and receptor

prioritization_tbl_ligand_receptor_human_mouse = prioritization_tbl_ligand_receptor_human_orthologs %>% rename(human_receptor = receptor, human_ligand = ligand)  %>% left_join(prioritization_tbl_ligand_receptor_mouse_orthologs %>% rename(receptor_mouse = receptor))
prioritization_tbl_ligand_receptor_human_mouse = prioritization_tbl_ligand_receptor_human_mouse %>% select(sender, receiver, human_ligand, human_receptor, mouse_ligand, mouse_receptor, human_score, mouse_score) %>% mutate_cond(is.na(human_score), human_score = 0)  %>% mutate_cond(is.na(mouse_score), mouse_score = 0) %>% mutate(conservation_score = mouse_score + human_score) %>% arrange(-conservation_score)

prioritization_tbl_ligand_receptor_human_mouse %>% filter(receiver == "KCs") %>% group_by(human_ligand) %>% top_n(1, conservation_score) %>% ungroup() %>% top_n(40, conservation_score) %>% View()

###### Extract the top ligands ############

receiver_oi = "KCs"
human_ligands = prioritization_tbl_ligand_receptor_human_mouse %>% filter(receiver == receiver_oi) %>% group_by(human_ligand) %>% top_n(1, human_score)  %>% ungroup() %>% top_n(50, human_score) %>% arrange(-human_score) %>% pull(human_ligand)
mouse_ligands = prioritization_tbl_ligand_receptor_human_mouse %>% filter(receiver == receiver_oi) %>% group_by(mouse_ligand) %>% top_n(1, mouse_score)  %>% ungroup() %>% top_n(50, mouse_score) %>% arrange(-mouse_score) %>% pull(mouse_ligand)

human_ligands_specific = human_ligands %>% setdiff(mouse_ligands %>% convert_mouse_to_human_symbols())

####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################
############## Ligand Activites  ############## ############## ############## ############## ##############
####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################

ligand_activities_df = readRDS("output/ligand_activities_humanKC_vs_mouseKC.rds")
ligand_activities_df_human_specific = ligand_activities_df %>% filter(test_ligand %in% human_ligands_specific) %>% arrange(-human_specific_vs_mousehumanKC) %>% filter(human_specific_vs_mousehumanKC > 0) %>%
  select(test_ligand, human_specific_vs_mousehumanKC) %>% rename(ligand = test_ligand, human_activity = human_specific_vs_mousehumanKC)


####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################
############## Circos plots  ############## ############## ############## ############## ##############
####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################

### show circos plot of the human-specific ligand-receptor pairs
### size of the arrow accords to the ligand activities! (just like transparancy?)

### now: only when keeping the conserved ones

receiver_oi = "KCs"

prioritized_tbl_human_conserved = prioritization_tbl_ligand_receptor_human %>% filter(receiver == receiver_oi) %>% group_by(ligand) %>% top_n(1, prioritization_score)  %>% ungroup() %>% top_n(50, prioritization_score)
prioritized_tbl_human_conserved_extended = prioritized_tbl_human_conserved %>% filter(receiver == receiver_oi) %>% filter(ligand %in% prioritized_tbl_human_conserved$ligand) %>% filter(prioritization_score >= min(prioritized_tbl_human_conserved$prioritization_score)) %>% group_by(ligand) %>% top_n(3, prioritization_score) %>% ungroup()

prioritized_tbl_all = prioritized_tbl_human_conserved_extended %>% distinct()

prioritized_tbl_all = prioritized_tbl_all %>% inner_join(prioritization_tbl_ligand_receptor_human %>% distinct(sender, receiver, niche) )
prioritized_tbl_all = prioritized_tbl_all %>% arrange(sender, prioritization_score)

prioritized_tbl_oi = prioritized_tbl_all %>% distinct() %>% filter(ligand %in% ligand_activities_df_human_specific$ligand)

prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique() %>% sort()

colors_sender = c("#EC67A7", "#FBB05F", "#A31A2A") %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("#5DA6DB")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

prioritized_tbl_oi = prioritized_tbl_oi %>% select(-prioritization_score) %>% inner_join(ligand_activities_df_human_specific) %>% rename(prioritization_score = human_activity)
prioritized_tbl_oi

cutoff = 0

arr.col = prioritized_tbl_oi %>% mutate(color = colors_sender[sender]) %>% select(ligand, receptor, color)

n_lsec_ligands = prioritized_tbl_oi  %>% group_by(sender) %>% count() %>% filter(sender == "LSECs") %>% pull(n)
n_hep_ligands = prioritized_tbl_oi %>% group_by(sender) %>% count() %>% filter(sender == "Hepatocytes") %>% pull(n)
n_stellate_ligands = prioritized_tbl_oi  %>% group_by(sender) %>% count() %>% filter(sender == "Stellate cells") %>% pull(n)

circos_output_KC = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver, cutoff, scale = FALSE, transparency = NULL, circos_type = "normal", border = TRUE)
circos_output_KC$p_circos

dev.off()
pdf("plots/Human_Mouse_Final/HUMAN_KC_conserved_subset_border.pdf", width = 10, height = 10)
circos_output_KC$p_circos
dev.off()

dev.off()
svg("plots/Human_Mouse_Final/HUMAN_KC_conserved_subset_border.svg", width = 10, height = 10)
circos_output_KC$p_circos
dev.off()

### SO PROTOCOL/
### TOP 50 diFF NICHNEET: in top 50 human not in mouse
### keep only ligands with positive ligand activity: more enrichment human KC vs mouse KC target genes
### size of the arrows and transparance dependent on ligand activity

### new plot of which the arrow size and transparancy does not depend on ligand activity

prioritized_tbl_oi = prioritized_tbl_oi %>% mutate_cond(sender == "LSECs", prioritization_score = rnorm(n_lsec_ligands, 0.90, 0.1)) %>% mutate_cond(sender == "Hepatocytes", prioritization_score = rnorm(n_hep_ligands, 0.75, 0.085)) %>% mutate_cond(sender == "Stellate cells", prioritization_score = rnorm(n_stellate_ligands, 0.70, 0.075))
prioritized_tbl_oi = prioritized_tbl_oi %>% group_by(sender, ligand) %>% mutate(prioritization_score = min(prioritization_score)) %>% ungroup()

circos_output_KC = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver, cutoff, scale = FALSE, transparency = NULL, circos_type = "normal", border = TRUE)
circos_output_KC$p_circos

dev.off()
pdf("plots/Human_Mouse_Final/HUMAN_KC_conserved_subset_border_v2.pdf", width = 10, height = 10)
circos_output_KC$p_circos
dev.off()

dev.off()
svg("plots/Human_Mouse_Final/HUMAN_KC_conserved_subset_border_v2.svg", width = 10, height = 10)
circos_output_KC$p_circos
dev.off()


####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################
############## Get the Human KC-specific target genes for all the shown ligands  ############## ############## ############## ############## ##############
####################################### ####################################### #######################################
#######################################  ####################################### ###############################
####################################### ####################################### #######################################
# Ligand-target matrix
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

humanKC_signature = readRDS("output/humanKC_specific_genes.rds")

active_ligand_target_links_df = prioritized_tbl_oi$ligand %>% lapply(get_weighted_ligand_target_links,geneset = humanKC_signature, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links_df = prioritized_tbl_oi %>% select(ligand, receptor) %>% inner_join(active_ligand_target_links_df) %>% rename(ligand_target_potential = weight) %>% inner_join(ligand_activities_df_human_specific) %>% arrange(-human_activity, -ligand_target_potential) %>% rename(ligand_activity_humanKC = human_activity)
active_ligand_target_links_df %>% xlsx::write.xlsx2("output/diff_nichenet_circos_HUMAN_SPECIFIC.xlsx",sheetName = "ligand_target", append = FALSE)
