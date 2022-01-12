make_ligand_activity_target_lfc_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, DE_table_receiver_processed_targets, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
  
  # ligand lfc
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score) %>% distinct()) %>% arrange(sender, ligand_score) 
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)
  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligands) 
  p1 = plot_data %>% 
    ggplot(aes(sender, ligand_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches") 
  max_lfc = abs(plot_data$ligand_score) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_ligand_lfc = p1 + custom_scale_fill + ylab("Prioritized Ligands") + xlab("Ligand LogFC in Sender")

  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()
  
  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score) 
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = DE_table_receiver_processed_targets %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))
  
  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>% 
  p1 = plot_data %>%
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_score)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Target:\nmin LFC vs\nother niches") + xlab("Target LogFC in Receiver")
  max_exprs = abs(plot_data$target_score ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.41, 0.4850, 0.5, 0.5150, 0.6, 0.7, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_targets = p1 + custom_scale_fill 
  
  
  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )
  
  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }
  
  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)
  
  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()
  
  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }
  
  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }
  
  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {
    
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }
  
  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  
  
  # Ligand-Activity-Scaled
  
  order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  
  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))
  
  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))
  # print(limits)
  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill
  
  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill
  
  
  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc )), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))
  
  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, n_groups, n_groups, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "SAaB
               ###C"
    combined_plot = patchwork::wrap_plots(S = p_ligand_lfc  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "SAaB
               L##C"
    
    combined_plot = patchwork::wrap_plots(S = p_ligand_lfc  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
  
}
make_ligand_zonation_activity_target_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
  
  # ligand expression
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation) %>% distinct()) %>% arrange(sender, ligand_score) 
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = exprs_tbl_ligand %>% inner_join(prioritization_tbl_ligand_receptor %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation))   %>% inner_join(ordered_ligands) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))
  
  
  
  p1 = plot_data  %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligands")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands = p1 + custom_scale_fill 
  p_ligands
  
  
  # ligand zonation
  p1 = plot_data  %>% filter(niche == "KC") %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Periportal-vs-Pericentral LFC") + xlab("Zonation LFC") + ylab("Prioritized Ligands")
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands_zonation = p1 + custom_scale_fill 
  p_ligands_zonation
  
  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()
  
  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score) 
  # ordered_targets = ordered_targets %>% select(-niche) %>% distinct() %>% mutate(niche = receiver)
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = exprs_tbl_target %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))
  plot_data = plot_data %>% group_by(target) %>% mutate(target_expression_scaled_myeloid = nichenetr::scaling_zscore(target_expression))
  
  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>% 
    p1 = plot_data %>%
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_expression_scaled_myeloid)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Scaled Exprs Target") + xlab("Target Expression")
  max_exprs = abs(plot_data$target_expression_scaled_myeloid ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_targets = p1 + custom_scale_fill 
  
  
  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )
  
  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }
  
  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)
  
  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()
  
  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }
  
  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }
  
  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {
    
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }
  
  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  # Ligand-Activity-Scaled
  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))
  
  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))

  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill
  
  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill
  
  
  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligands)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligands_zonation)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))
  
  
  
  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, 3, n_groups, n_groups, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "SZAaB
               ####C"
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          Z = p_ligands_zonation + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""), 
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "SZAaB
               L###C"
    
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          Z = p_ligands_zonation + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""), 
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
  
}


make_ligand_activity_target_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
  
  # ligand expression
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score) %>% distinct()) %>% arrange(sender, ligand_score) 
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligands) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))
  
  p1 = plot_data  %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligands")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands = p1 + custom_scale_fill 
  p_ligands
  
  
  
  
  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()
  
  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score) 
  
  # if duplicated: eg MoMac1-vs-MoMac1_CV
  ordered_targets = ordered_targets %>% select(-niche) %>% distinct() %>% mutate(niche = receiver)
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = exprs_tbl_target %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))
  plot_data = plot_data %>% group_by(target) %>% mutate(target_expression_scaled_myeloid = nichenetr::scaling_zscore(target_expression))
  
  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>% 
    p1 = plot_data %>% 
    
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_expression_scaled_myeloid)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Scaled Exprs Target") + xlab("Target Expression")
  max_exprs = abs(plot_data$target_expression_scaled_myeloid ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_targets = p1 + custom_scale_fill 
  
  
  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )
  
  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }
  
  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)
  
  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()
  
  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }
  
  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0) 
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)
    
    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }
  
  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {
    
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }
  
  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  
  order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  
  # Ligand-Activity-Scaled
  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  # print(ligand_pearson_df)
  
  
  
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))
  
  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))
  # print(limits)
  
  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill
  
  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill
  
  
  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligands)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))
  
  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, n_groups, n_groups, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "SAaB
               ###C"
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "SAaB
               L##C"
    
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
  
}

make_ligand_receptor_lfc_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = TRUE, heights = NULL, widths = NULL){
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score) 
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score) 
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)

    plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors) 
  p_lig_lfc = plot_data %>% 
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_lig_lfc = p_lig_lfc + custom_scale_fill 
  p_lig_lfc
  
  p_rec_lfc = plot_data %>% 
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches")  + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_rec_lfc = p_rec_lfc + custom_scale_fill 
  p_rec_lfc
  
  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc, nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}


make_ligand_receptor_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, exprs_tbl_ligand, exprs_tbl_receiver, plot_legend = TRUE, heights = NULL, widths = NULL){
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score) 
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score) 
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))
  
  p1 = plot_data  %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligand-Receptor pairs")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_lig_lfc = p1 + custom_scale_fill 

  plot_data = exprs_tbl_receptor %>% inner_join(ordered_ligand_receptors) %>% filter(receiver %in% (prioritization_tbl_ligand_receptor$receiver %>% unique()))
  plot_data = plot_data %>% group_by(receptor) %>% mutate(receptor_expression_scaled_sender = nichenetr::scaling_zscore(receptor_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver))
  
  p1 = plot_data  %>% 
    # ggplot(aes(receptor_ordered, sender , color = receptor_expression_scaled_myeloid, size = receptor_fraction )) +
    ggplot(aes(receiver,ligand_receptor_ordered , fill = receptor_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs receptor") + xlab("Receptor Expression") 
  max_exprs = abs(plot_data$receptor_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_rec_lfc = p1 + custom_scale_fill 
  p_rec_lfc
  
  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc, nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}
make_ligand_receptor_exprs_zonation_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, exprs_tbl_ligand, exprs_tbl_receiver, plot_legend = TRUE, heights = NULL, widths = NULL){
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score) 
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score) 
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)
  
  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))
  
  ## ligand figure
  p1 = plot_data  %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligand-Receptor pairs")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_lig_lfc = p1 + custom_scale_fill 
  
  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver, ligand, ligand_score_zonation))
  
  
  # ligand zonation
  p1 = plot_data  %>% filter(niche == "KC") %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nPeriportal-vs-Pericentral\nLigand LFC") + xlab("Sender Zonation\nLFC ligand") 
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands_zonation = p1 + custom_scale_fill 
  # p_ligands_zonation
  
  ## receptor fugure
  plot_data = exprs_tbl_receptor %>% inner_join(ordered_ligand_receptors) %>% filter(receiver %in% (prioritization_tbl_ligand_receptor$receiver %>% unique()))
  plot_data = plot_data %>% group_by(receptor) %>% mutate(receptor_expression_scaled_sender = nichenetr::scaling_zscore(receptor_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver))
  
  p1 = plot_data  %>% 
    # ggplot(aes(receptor_ordered, sender , color = receptor_expression_scaled_myeloid, size = receptor_fraction )) +
    ggplot(aes(receiver,ligand_receptor_ordered , fill = receptor_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs receptor") + xlab("Receptor Expression") 
  max_exprs = abs(plot_data$receptor_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_rec_lfc = p1 + custom_scale_fill 
  # p_rec_lfc
  
  design = "AZ#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, Z = p_ligands_zonation + ylab(""), B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 3, 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}


make_ligand_receptor_lfc_zonation_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = TRUE, heights = NULL, widths = NULL){
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score) 
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score) 
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)
  
  
  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors) 
  p_lig_lfc = plot_data %>% 
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in sender")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_lig_lfc = p_lig_lfc + custom_scale_fill 
  # p_lig_lfc
  
  # ligand zonation
  p1 = plot_data   %>% filter(niche == "KC") %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nPeriportal-vs-Pericentral\nLigand LFC") + xlab("Sender Zonation\nLFC ligand") 
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands_zonation = p1 + custom_scale_fill 
  
  
  p_rec_lfc = plot_data %>% 
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches") + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_rec_lfc = p_rec_lfc + custom_scale_fill 
  # p_rec_lfc
  
  design = "AZ#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, Z = p_ligands_zonation + ylab(""), B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 3, 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}
############################## ############################## ############################## ############################## 
############################## smaller functions############################## ############################## 
############################## ############################## ############################## ############################## 

ligand_lfc_plot = function(plot_data, max_lfc){
  p_lig_lfc = plot_data %>% 
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in sender")
  
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  p_lig_lfc = p_lig_lfc + custom_scale_fill 
  
  return(p_lig_lfc)
  
}

ligand_lfc_zonation_plot = function(plot_data, max_lfc){
  # ligand zonation
  p1 = plot_data   %>% filter(niche == "KC") %>% 
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nPeriportal-vs-Pericentral\nLigand LFC") + xlab("Sender Zonation\nLFC ligand") 
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  
  p_ligands_zonation = p1 + custom_scale_fill 
  
  return(p_ligands_zonation)
  
}

receptor_lfc_plot = function(plot_data, max_lfc){
  p_rec_lfc = plot_data %>% 
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    # facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches") + xlab("Receptor LFC\n in Receiver")
  # max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_rec_lfc = p_rec_lfc + custom_scale_fill 
  # p_rec_lfc
  return(p_rec_lfc)
}
normalized_activity_plot = function(plot_data, max_activity, min_activity){
  p_lig_lfc = plot_data %>% 
    ggplot(aes(receiver, ligand_receptor_ordered, fill = activity_normalized)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Normalized Ligand activity\nin Receiver")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Normalized Ligand activity\nin Receiver")
  
  limits = c(-max_normalized_activity, max_activity)
  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.49, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_activity = p_lig_lfc +custom_scale_fill
  
  return(p_activity)
  
}
activity_plot = function(plot_data, max_activity, min_activity){
  p_lig_lfc = plot_data %>% 
    ggplot(aes(receiver, ligand_receptor_ordered, fill = activity)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand activity\nin Receiver")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand activity\nin Receiver")
  
  limits = c(min_activity, max_activity)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = limits) # non-scaled
  p_activity = p_lig_lfc +custom_scale_fill
  
  return(p_activity)
  
}
prioritization_score_plot = function(plot_data){
  p_score = plot_data %>% 
    ggplot(aes(scoretype   , ligand_receptor_ordered, fill = score)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "LR Prioritization Scores")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Prioritization Scores\nLR pair")
  
  limits = c(0, 1.01)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "YlGn")),values = c(0, 0.50, 0.625, 0.75, 0.825, 0.885, 0.945, 1),  limits = limits) # non-scaled
  p_score = p_score +custom_scale_fill
  
  return(p_score)
  
}
