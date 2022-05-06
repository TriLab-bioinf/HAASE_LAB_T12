# Plotting functions
draw_enrichment_barplot <- function(my_enricher_obj, 
                                    my_pathway_counts = 10, 
                                    file_name_prefix = "enricher_barplot", 
                                    my_width = 11, my_height = 8){
  my_enricher_obj@result$p.adjust <- as.numeric(format(my_enricher_obj@result$p.adjust,  digits=3))
  sp <- ggbarplot(my_enricher_obj@result[my_pathway_counts, ], 
                  x = "Description", 
                  y = "Count",
                  fill = "p.adjust",          # change fill color by cyl
                  color = "white",            # Set bar border colors to white
                  sort.val = "desc",          # Sort the value in dscending order
                  sort.by.groups = FALSE,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  rotate = TRUE,
                  ggtheme = theme_minimal(),
                  ylab = c("Gene counts"),
  ) + gradient_fill("RdYlBu")
  
  print(sp)
  
  ggsave2(filename = paste0(file_name_prefix,".pdf"), plot = sp, width = my_width, height = my_height)
  
  return(sp)
}

draw_GSEA_barplot <- function(my_gsea_obj, 
                              my_pathway_counts = 10, 
                              file_name_prefix = "gsea_barplot", 
                              my_width = 11, my_height = 8){
  my_gsea_obj@result$p.adjust <- as.numeric(format(my_gsea_obj@result$p.adjust,  digits=3))
  sp <- ggbarplot(head(my_gsea_obj@result, n = my_pathway_counts), 
                  x = "Description", 
                  y = "NES",
                  fill = "p.adjust",          # change fill color by cyl
                  color = "white",            # Set bar border colors to white
                  sort.val = "desc",          # Sort the value in dscending order
                  sort.by.groups = FALSE,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  rotate = TRUE,
                  ggtheme = theme_minimal(),
                  ylab = c("Normalized Enrichment Score (NES)"),
                  lab.size = 3
  ) + gradient_fill("RdYlBu")
  print(sp)
  ggsave2(filename = paste0(file_name_prefix,".pdf"), plot = sp, width = my_width, height = my_height)
  return(sp)
}

# Function for overrepresentation analysis (ONLY FOR HUMAN GENES SO FAR)
run_overrepresentation_analysis <- function(gene_set, genes, 
                                            analysis_name = "dds.result", 
                                            gs_name = "MSigDB", 
                                            adj_p_cutoff = 0.05, 
                                            type = "reactome"){
  # Inputs: gene_set=Gene set from MSigDB, 
  # genes=vector with ensembl gene IDs.
  # analysis_name=prefix describing DE comparison for output file names
  # gs_name=identifier of the MSigDB used (gene_set)
  
  # Create output directories
  my_path = paste0("./ClusterProfiler/",analysis_name,"/")
  dir.create(path = my_path, recursive = TRUE)
  
  if(length(genes > 0)){
    if(type == "reactome"){
      # Replace ensembl IDs by Entrez IDs for compatibility with gsePathway
      entrez_ids.df <- clusterProfiler::bitr(geneID = genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
      msig_h.oa <- enrichPathway(gene = entrez_ids.df$ENTREZID, qvalueCutoff = adj_p_cutoff, readable = TRUE)
    } else {
      msig_h.oa <- enricher(gene = genes, TERM2GENE=gene_set,  pvalueCutoff = adj_p_cutoff)
    }
    
    if(dim(msig_h.oa@result)[1] > 0){
      # Convert Entrez gene IDs to gene symbols
      msig_h.oa.gs <- setReadable(msig_h.oa, org.Hs.eg.db, keyType = "ENSEMBL")
      
      if(dim(msig_h.oa.gs)[1] > 0){
        
        # Write output to a file (OR stands for Overrepresenation analysis). Change file names for diff conditions.
        write.table(as.data.frame(msig_h.oa.gs), file = paste0(my_path,"OR_msigdb.",gs_name,".",analysis_name,".txt"), sep = "\t")
        
        # Dotplot
        p1 <- dotplot(msig_h.oa.gs, showCategory=36) + 
          ggtitle(paste0(analysis_name," dotplot for ORA")) + 
          theme(axis.text=element_text(size=3))
        
        ggsave2(filename = paste0(my_path,"OA_msigdb.",gs_name,".",analysis_name,"_dotplot.pdf"), plot = p1, width = 11, height = 8)
        
        # Barplot
        p3 <- draw_enrichment_barplot(my_enricher_obj = msig_h.oa.gs, 
                                      my_pathway_counts = 36, 
                                      file_name_prefix = paste0(my_path,
                                                                "OA_msigdb.",
                                                                gs_name,".",
                                                                analysis_name,
                                                                "_barplot"
                                                                ), 
                                      my_width = 11, 
                                      my_height = 11)
        print(p1)
        print(p3)
        
      }
    }  
  }
}

# Function for enrichment analysis (ONLY FOR HUMAN GENES SO FAR)
run_enrichment_analysis <- function(gene_set, geneList, analysis_name = "dds.result", gs_name = "MSigDB", adj_p_cutoff = 0.05, type = "general"){
  # Inputs: gene_set=Gene set from MSigDB, 
  # geneList=a list with ensembl_ids as names and decremented order of Log2FCs as values
  # analysis_name=prefix describing DE comparison for output file names
  # gs_name=identifier of the MSigDB used (gene_set)
  
  # Create output directories
  my_path = paste0("./ClusterProfiler/",analysis_name,"/")
  dir.create(path = my_path, recursive = TRUE)
  
  # Run Gene Set Enrichment Analysis. Reactome type display reactome IDs and Descriptions properly
  if(type == "reactome"){
    # Replace ensembl IDs by Entrez IDs for compatibility with gsePathway
    entrez_ids.df <- clusterProfiler::bitr(geneID = names(geneList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
    geneList.entrez <- geneList[c(entrez_ids.df$ENSEMBL)]
    names(geneList.entrez) <- entrez_ids.df$ENTREZID
    msig_h.gsea <- gsePathway(geneList = geneList.entrez, pvalueCutoff = adj_p_cutoff, verbose = TRUE, eps = 0)
    if(dim(msig_h.gsea@result)[1] > 0){ 
      # Convert Entrez gene IDs to gene symbols
      msig_h.gsea.gs <- setReadable(msig_h.gsea, org.Hs.eg.db, keyType = "ENTREZID")
    } 
  } else {
    msig_h.gsea <- GSEA(geneList, TERM2GENE = gene_set, pvalueCutoff = adj_p_cutoff, verbose = TRUE, eps = 0)
    if(dim(msig_h.gsea@result)[1] > 0){ 
      # Convert Entrez gene IDs to gene symbols
      msig_h.gsea.gs <- setReadable(msig_h.gsea, org.Hs.eg.db, keyType = "ENSEMBL")
    }  
    
  } 
  
  if(dim(msig_h.gsea@result)[1] > 0){ 
    # Write output to a file (OR stands for Overrepresenation analysis). Change file names for diff conditions.
    write.table(as.data.frame(msig_h.gsea.gs), file = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,".txt"), sep = "\t")
    
    # Generate plots summarizing enrichment results (You might need to tweak the height and width parameters of the plot to make it fit within the page limits and look nice, otherwise it might shrink to fit within the page)
    
    # Upset plots
    # filter out signatures with NES < 0 since they are likely due to RNAseL degradation to simplify the plot
    msig_h.gsea_bt_0 <- filter(msig_h.gsea, NES > 0) # filter out signatures with NES < 0
    msig_h.gsea_lt_0 <- filter(msig_h.gsea, NES < 0) # filter out signatures with NES < 0
    
    if(dim(msig_h.gsea_bt_0@result)[1] > 0){
      p_ups_up <- upsetplot(msig_h.gsea_bt_0, n=36) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_UP_upset_plot.pdf"), 
             plot = p_ups_up, width = 22, 
             height = 4 * (max(2,min(36, length(msig_h.gsea_bt_0@result$ID))/15)))
    }
    if(dim(msig_h.gsea_lt_0@result)[1] > 0){
      p_ups_down <- upsetplot(msig_h.gsea_lt_0, n=36) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_DOWN_upset_plot.pdf"), 
             plot = p_ups_down, width = 22, 
             height = 4 * (max(2,min(36, length(msig_h.gsea_lt_0@result$ID))/15)))
    }
    
    # Dotplot
    p2 <- dotplot(msig_h.gsea.gs, showCategory=36, 
                  font.size = 7, 
                  color = "p.adjust", 
                  label_format = 100) + 
      ggtitle(paste0(my_prefix," dotplot for GSEA"))
    
    ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_dotplot.pdf"), 
            plot = p2, width = 15, height = 8)
    
    # Barplot
    p4 <- draw_GSEA_barplot(my_gsea_obj =  msig_h.gsea.gs, my_pathway_counts = 36, 
                            file_name_prefix = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_barplot"), 
                                                my_width = 11, 
                                                my_height =11
                            )
    
    # Network plots and tree plots for all, up and downregulated pathways
    if(dim(msig_h.gsea)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix <- pairwise_termsim(msig_h.gsea, showCategory = 400)
      
      # Network plot
      p6 <- emapplot(distance_matrix, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network.pdf"), 
              plot = p6, width = 11, height = 8)
      
      # Treeplots
      number_of_categories = min(80, length(distance_matrix@result$ID))
      p12 <- treeplot(distance_matrix, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree.pdf"), 
              plot = p12, width = 11, height = 8 * (number_of_categories/40))
    }
    
    if(dim(msig_h.gsea_bt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.up <- pairwise_termsim(msig_h.gsea_bt_0, showCategory = 400)
      
      # Network plot
      p8 <- emapplot(distance_matrix.up, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network_UP.pdf"), 
              plot = p8, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.up = min(80, length(distance_matrix.up@result$ID))
      p14 <- treeplot(distance_matrix.up, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.up), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree_UP.pdf"), 
              plot = p14, width = 11, height = 8 * (number_of_categories.up/40))
      
    }
    
    if(dim(msig_h.gsea_lt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.down <- pairwise_termsim(msig_h.gsea_lt_0, showCategory = 400)
      
      # Network plot
      p10 <- emapplot(distance_matrix.down, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network_DOWN.pdf"), 
              plot = p10, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.down = min(80, length(distance_matrix.down@result$ID))
      p16 <- treeplot(distance_matrix.down, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.down), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree_DOWN.pdf"), 
              plot = p16, width = 11, height = 8 * (number_of_categories.down/40))
      
    }
    
    
    
    #print(p2)
    #print(p4)
    #print(p_ups_up)
    #print(p_ups_down)
    #print(p6)
    #print(p8)
    #print(p10)
    #print(p12)
    #print(p14)
    #print(p16)
  }
}
