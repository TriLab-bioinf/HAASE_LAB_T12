# Define function to generate plots

sort_and_write_res_table <- function(result_table, file_name){
  dir.create(path = "./DE", showWarnings = FALSE)
  # Sort genes by (padj)
  result_table_sorted <- result_table[order(result_table$padj, decreasing = FALSE),]
  # Add gene symbols
  gene_list <- rownames(result_table_sorted)
  result_table_sorted$gene_name <- ensembl_to_symbol$gene_name[match(gene_list, ensembl_to_symbol$Ensembl_ID)]
  #df <-as.data.frame(cbind(result_table_sorted, Gene_name = symbol_list))
  
  # Write sorted table to file
  write.table(result_table_sorted, file = paste0("./DE/",file_name,".txt"), 
              sep = "\t", col.names=NA)
  return(result_table_sorted)
}

plot_gene <- function(my_dds, gene_annotation, gene_name){
  
  # retrieve gene_id from gene_name
  keep <- gene_annotation$gene_name == gene_name
  gene_id     <- gene_annotation[keep, "Ensembl_ID"]
  
  
  d <- plotCounts(my_dds, gene = gene_id, intgroup="group", returnData=TRUE)
  
  ggplot(d, aes(x=group, y=count, colour=group)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(25,100,400)) + ggtitle(paste(gene_id, gene_name))
  
  ggboxplot(d, x = "Diff", 
            y = "count",
            color = "Diff", 
            palette =c("#00AFBB", "#E7B800", "#FC4E07", "#CC00FF", "#32DCFF"),
            add = "jitter", 
            shape = "Diff",
            yscale = "log10",
            title = paste(gene_id , gene_name)
  )
  
  ggviolin(d, x = "Diff", 
           y = "count", 
           fill = "Diff",
           palette = c("#00AFBB", "#E7B800", "#FC4E07", "#CC00FF", "#32DCFF"),
           add = "boxplot", 
           add.params = list(fill = "white"),
           yscale = "log10",
           title = paste(gene_id, gene_name)
  )
}

replace_gene_acc_by_symbols_in_dds <- function(res.tmp){
  ensemble_ids <- rownames(res.tmp)
  symbols <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = c('SYMBOL'), keytype = 'ENSEMBL') #, multiVals = "first"
  symbols <- symbols[!is.na(symbols)]
  to_name <- rownames(res.tmp) %in% names(symbols)
  res.tmp@rownames[to_name] <- as.vector(symbols)
  return(res.tmp)
}

# Replace ensembl_ids by genbank gene_ids.
replace_gene_acc_by_gb_ids <- function(ensemble_ids, return_all = TRUE){
  entrezids <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = c('ENTREZID'), keytype = 'ENSEMBL', multiVals = "first")
  entrezids <- entrezids[!is.na(entrezids)]
  to_name <- ensemble_ids %in% names(entrezids)
  ensemble_ids[to_name] <- as.vector(entrezids)
  if (return_all){
    return(ensemble_ids)
  }
  else {
    return(ensemble_ids[to_name])
  }
}

replace_gene_acc_by_symbol_ids <- function(ensemble_ids, return_all = TRUE){
  entrezids <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = c('SYMBOL'), keytype = 'ENSEMBL', multiVals = "first")
  entrezids <- entrezids[!is.na(entrezids)]
  to_name <- ensemble_ids %in% names(entrezids)
  ensemble_ids[to_name] <- as.vector(entrezids)
  if (return_all){
    return(ensemble_ids)
  }
  else {
    return(ensemble_ids[to_name])
  }
}

replace_ensemble_acc_by_symbols <- function(ensemble_ids){
  symbols <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = c('SYMBOL'), keytype = 'ENSEMBL', multiVals = "first")
  symbols <- symbols[!is.na(symbols)]
  to_name <- ensemble_ids %in% names(symbols)
  ensemble_ids[to_name] <- as.vector(symbols)
  return(ensemble_ids)
}
generate_volcano_plot <- function(res.tmp, my_file_name, log_scale = FALSE){
  res.tmp <- replace_gene_acc_by_symbols_in_dds(res.tmp)
  vp <- EnhancedVolcano(res.tmp, 
                        lab = rownames(res.tmp), 
                        x = 'log2FoldChange', 
                        y = 'pvalue',
                        pointSize = 1,
                        colAlpha = 4/5,
                        labSize = 2,  # Controls labels size
                        title = res.tmp@elementMetadata$description[2],
                        titleLabSize = 14,
                        subtitle = '', # add subtitle here
                        subtitleLabSize = 12,
                        legendPosition = 'right',
                        legendLabSize = 12,
                        legendIconSize = 2.0,
                        axisLabSize = 10
  )
  
  if (log_scale){
    vp <- vp + scale_x_log10()
  }
  ggsave(filename = paste0(my_file_name,".pdf"), plot = vp )
  
  print(vp)
}
generate_volcano_plot_with_ids <- function(res.tmp, my_file_name, log_scale = FALSE, gene_list){
  res.tmp <- replace_gene_acc_by_symbols_in_dds(res.tmp)
  vp <- EnhancedVolcano(res.tmp, 
                        lab = rownames(res.tmp), 
                        x = 'log2FoldChange', 
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        pointSize = 1,
                        colAlpha = 4/5,
                        labSize = 2,  # Controls labels size
                        labCol = "black",
                        title = res.tmp@elementMetadata$description[2],
                        titleLabSize = 10,
                        subtitle = '', # add subtitle here
                        subtitleLabSize = 10,
                        legendPosition = 'right',
                        legendLabSize = 10,
                        legendIconSize = 4.0,
                        axisLabSize = 10,
                        drawConnectors = TRUE,
                        selectLab = gene_list, # vector of gene symbols to label on volcanoplot
                        boxedLabels = FALSE
  )
  
  if (log_scale){
    vp <- vp + scale_x_log10()
  }
  ggsave(filename = paste0(my_file_name,".pdf"), plot = vp )
  
  print(vp)
}
plot_heat_map <- function( my_vstd, gene_list, file_name, variables){
  
  # Replace gene IDs by gene symbols
  gene_list <- replace_ensemble_acc_by_symbols(gene_list)
  rownames(my_vstd) <- replace_ensemble_acc_by_symbols(rownames(my_vstd))
  
  # Plot the heat map
  hmp <- pheatmap(assay(my_vstd)[gene_list,], cluster_rows=T, show_rownames=TRUE,
                  cluster_cols=T, annotation_col = variables, fontsize_col = 5, tile = file_name)
  
  ggsave(filename = paste0(file_name,"_heatmap.pdf"), plot = hmp, width = 8.5, height = 11, units = "in")
  
  print(hmp)
}
### It will take list of gene symbols as gene_list
plot_heat_map_from_gene_symbols <- function( my_vstd, gene_list, file_name, variables){
  
  # Replace gene IDs by gene symbols
  #gene_list <- replace_ensemble_acc_by_symbols(gene_list)
  rownames(my_vstd) <- replace_ensemble_acc_by_symbols(rownames(my_vstd))
  
  # Plot the heat map
  hmp <- pheatmap(assay(my_vstd)[gene_list,], cluster_rows=T, show_rownames=TRUE,
                  cluster_cols=T, annotation_col = variables, fontsize_col = 5, tile = file_name)
  
  ggsave(filename = paste0(file_name,"_heatmap.pdf"), plot = hmp, width = 8.5, height = 11, units = "in")
  
  print(hmp)
}

# Function to remove all-zero rows (by adjusting min_total_count the function can filter out rows based on total counts other than 0)
remove_all_zero_rows <- function(df, min_total_count = 0){
  df <- df[rowSums(df) > min_total_count,]
  return(df)
}


# Function to normalize by TPMs based on transcript length
# Normalize counts to TPMs
# Fetch exon length from BioMart
normalize_by_TPM <- function(counts.df) {
  
  listMarts(host="https://uswest.ensembl.org/")
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  
  my_mart <- getBM(attributes=c("ensembl_gene_id","exon_chrom_start","exon_chrom_end"), 
                   mart = ensembl, 
                   filters = c("ensembl_gene_id"), 
                   values = rownames(counts.df)
  )
  my_mart$exon_length <- (abs(my_mart$exon_chrom_start - my_mart$exon_chrom_end) + 1) / 1000 # length in Kb
  
  # Calculate transcript length
  transcript_lengths <- aggregate(my_mart$exon_length, by=list(Category=my_mart$ensembl_gene_id), FUN=sum)
  
  # Eliminate gene IDs from counts.df without transcript length info in transcript_lengths
  transcript_lengths <- transcript_lengths[transcript_lengths$Category %in% rownames(counts.df),]
  counts.df <- counts.df[transcript_lengths$Category,]
  
  # Sort transcripts_length df by rownames of counts.df
  transcript_lengths <- transcript_lengths[match(rownames(counts.df), transcript_lengths$Category), ]

  # See reference for formula
  # https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
  x.df <- apply(counts.df, MARGIN = 2, FUN = function(x) x/transcript_lengths$x)
  x.df <- apply(x.df, MARGIN = 2, FUN = function(x) x * 1e6 / sum(x)) 
  return(x.df)
}

print_gene_expression <- function(dds, gene_id, symbol = "", title_prefix = ""){
  counts.p <- plotCounts(dds, gene = gene_id, intgroup = c("Inducer", "Genotype"), 
                         transform = FALSE, returnData = TRUE)
  p <- ggerrorplot(counts.p, x = "Inducer", y = "count",
                   desc_stat = "mean_sd", 
                   color = "Inducer",
                   add = "dotplot", 
                   add.params = list(color = "darkgray"), 
                   facet.by = "Genotype",
                   title = paste(title_prefix,"Gene",gene_id, symbol), 
                   ylab = "Normalized counts"
  )
  return(p)
}

go_enrichment <- function(geneList, my_file, golevel = 3, 
                          onthology = "BP", qvalue = 0.05, 
                          fdr = 0.05, log2FC_cutoff = 0){
  # Convert Ensembl IDs to EntrezIDs
  rownames(geneList) <- replace_gene_acc_by_gb_ids(rownames(geneList),
                                                   return_all = TRUE)
  # Discard genes missing Entrez Gene IDs
  keep <- grep('ENSG', rownames(geneList), invert = TRUE)
  geneList <- geneList[keep,]
  
  geneList.filtered <- rownames(geneList)[abs(geneList$log2FoldChange) > log2FC_cutoff 
                                          & geneList$padj < qvalue]
  ego <- enrichGO(gene = geneList.filtered,
                  OrgDb= "org.Hs.eg.db",
                  ont = onthology, # Either MF, BP, CC or ALL
                  pAdjustMethod = "BH",
                  qvalueCutoff  = fdr,
                  keyType = 'ENTREZID',
                  readable      = TRUE,
                  minGSSize = 10,
                  maxGSSize = 500
  )
  write.table(ego, file = paste0(my_file,"_",onthology,"_go_enrich.txt"), sep = "\t", col.names=NA)
  return(ego)
}

go_overrep <- function(geneList, my_file, golevel = 3, onthology = "BP"){
  # Convert Ensembl IDs to EntrezIDs
  rownames(geneList) <- replace_gene_acc_by_gb_ids(rownames(geneList),
                                                   return_all = TRUE)
  # Discard genes missing Entrez Gene IDs
  keep <- grep('ENSG', rownames(geneList), invert = TRUE)
  geneList <- geneList[keep,]
  
  geneList.filtered <- rownames(geneList)[abs(geneList$log2FoldChange) > log2FC_cutoff & geneList$padj < alpha]
  
  print(paste( "Total number of Entrez gene IDs = ",length(geneList.filtered)))
  
  ggo <- groupGO(gene     = geneList.filtered,
                 OrgDb    = "org.Hs.eg.db",
                 ont      = onthology,  # Either MF, BP or CC
                 level    = golevel, # The higher the number, the more specific the GO terms.
                 readable = TRUE, # if readable is TRUE, the gene IDs will mapping to gene symbols.
                 keyType = 'ENTREZID'
  )
  write.table(ggo, file = paste0(my_file,"_",onthology,"_go_over.txt"), sep = "\t",col.names=NA)
  return(ggo)
}


