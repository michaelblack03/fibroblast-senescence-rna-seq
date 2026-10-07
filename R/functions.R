#FUNCTIONS

#function for creating faceted expression-density plots
plot_expression_density = function(em, out_path)
{
  #melt em table 
  em.m = melt(em)
  
  #draw plot
  ggp = ggplot(em.m, aes(x = log10(value)))+
    geom_density(colour="coral2", fill="coral", size = 0.8, alpha=0.3)+
    facet_wrap(~variable,ncol = 3)+
    labs(x="expression (log10)", y="density", title = "")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, linewidth=1),
          strip.text.x = element_text(size=12, family="Arial", face="bold", vjust=1),
          panel.spacing = unit(1, "lines"))
  
  svglite(out_path, height = 6, width = 6)
  print(ggp)
  dev.off()
  
  return(ggp)
}


#plotting PCA function: 
plot_pca = function(pca, pca_coordinates, pcx, pcy, out_path)
{
   #calculate %var for axes
   vars = apply(pca$x, 2, var)
   prop_x = round(vars[pcx] / sum(vars), 4) * 100
   prop_y = round(vars[pcy] / sum(vars), 4) * 100

   #save axes titles as object
   x_axis_label = paste(pcx, " (",prop_x, "%)",sep="")
   y_axis_label = paste(pcy, " (",prop_y, "%)",sep="")
   
   #draw plot
   ggp = ggplot(pca_coordinates, aes_string(x=pcx, y=pcy, colour="SAMPLE_GROUP"))+
     geom_point(size=3)+
     geom_text_repel(aes(label=SAMPLE_GROUP), colour="black", size=2)+
     scale_colour_manual(values=c("bisque2", "cadetblue3", "coral2"),labels = c("Prolif", "Senes", "Senes_MTKO"),name="")+
     labs(x=x_axis_label, y=y_axis_label)+
     theme_bw()+
     theme(panel.grid = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(colour="black", fill=NA, linewidth=1),
           axis.title = element_text(size=14),
           axis.text = element_text(size=10))
  
   svglite(out_path, height = 6, width = 6)
   print(ggp)
   dev.off()
    
  return(ggp)
}


#plotting heatmap function:
plot_heatmap = function(scaled_table, out_path)
{
  #make hm.matrix
  hm.matrix = as.matrix(scaled_table)
  
  #cluster x- and y-axes
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist,method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN = "average")
  y.order = order.dendrogram(y.dd.reorder)
  
  x.dist = Dist(t(hm.matrix),method="spearman")
  x.cluster = hclust(x.dist,method = "average")
  x.dd = as.dendrogram(x.cluster)
  x.dd.reorder = reorder(x.dd,0,FUN = "average")
  x.order = order.dendrogram((x.dd.reorder))
  
  #reorder and melt
  hm.matrix_clustered = hm.matrix[y.order, x.order]
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  #set colour palette
  colours = c("bisque", "coral4", "coral2")
  palette = colorRampPalette(colours)(50)
  
  #plot heatmap
  ggp = ggplot(hm.matrix_clustered, aes(x=Var1, y=Var2, fill = value))+
    geom_tile()+
    scale_fill_gradientn(colours = palette)+ 
    ylab("Genes")+ 
    xlab("")+ 
    theme( axis.ticks=element_blank(), 
          legend.title = element_blank(), 
          legend.spacing.x = unit(0.25, 'cm'),
            axis.text.x = element_text(size = 3, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 2))
  
  svglite(out_path, height = 6, width = 6)
  print(ggp)
  dev.off()
  
  return(ggp)
}


#single-gene boxplot function:
plot_single_boxplot = function(gene_name, em_scaled, sample_groups, out_path)
{
  gene_data = em_scaled[,gene_name]
  gene_data = data.frame(gene_data)
  gene_data$sample_group = sample_groups
  names(gene_data) = c("expression", "sample_group")
  
  ggp = ggplot(gene_data, aes(x = sample_group, y = expression, fill = sample_group))+
    geom_boxplot()+
    scale_colour_manual(values=c("black", "black", "black"))+
    scale_fill_manual(values = c("coral2", "cadetblue3", "bisque2"), labels=c("Prolif", "Senes", "Senes_MTKO"), name="")+
    labs(x = "Sample Group", y = "Expression")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, linewidth=1),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
  
  svglite(out_path, height = 6, width = 6)
  print(ggp)
  dev.off()
  
  return(ggp)
}


#multi-gene boxplot function:
plot_multi_boxplot = function(top6table, em_scaled, sample_groups, out_path)
{
  #get gene data
  candidate_genes = row.names(top6table)
  gene_data = em_scaled[,candidate_genes]
  gene_data = data.frame(gene_data)
  gene_data$sample_group = sample_groups
  
  #melt gene data table
  gene_data.m = melt(gene_data, id.vars = "sample_group")
  
  #create boxplot
  ggp = ggplot(gene_data.m, aes(x = sample_group, y = value, fill = sample_group))+
    geom_boxplot()+
    scale_colour_manual(values=c("black", "black", "black"))+
    scale_fill_manual(values = c("coral2", "cadetblue3", "bisque2"),
                      labels=c("Prolif", "Senes", "Senes_MTKO"), name="")+
    labs(x = "Sample Group", y = "Expression")+
    facet_wrap(~variable, ncol=3)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour="black", fill=NA, linewidth=1),
      strip.text.x = element_text(size=12, face="bold", vjust=1),
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
    )
  
  svglite(out_path, height = 6, width = 6)
  print(ggp)
  dev.off()
  
  return(ggp)
}

#have to make separate multi-gene boxplot function for ORA results as gene input will be a list, not a data frame
#if i had more time i would try and implement if/else code to handle different data formats
#i don't, so right now it is just easier to copy code and remove one line for a new, separate function
#multi-gene boxplot function:
plot_multi_boxplot_ora = function(candidate_genes, em_scaled, sample_groups, out_path)
{
  #get gene data
  gene_data = em_scaled[,candidate_genes]
  gene_data = data.frame(gene_data)
  gene_data$sample_group = sample_groups
  
  #melt gene data table
  gene_data.m = melt(gene_data, id.vars = "sample_group")
  
  #create boxplot
  ggp = ggplot(gene_data.m, aes(x = sample_group, y = value, fill = sample_group))+
    geom_boxplot()+
    scale_colour_manual(values=c("black", "black", "black"))+
    scale_fill_manual(values = c("coral2", "cadetblue3", "bisque2"),
                      labels=c("Prolif", "Senes", "Senes_MTKO"), name="")+
    labs(x = "Sample Group", y = "Expression")+
    facet_wrap(~variable, ncol=3)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour="black", fill=NA, linewidth=1),
      strip.text.x = element_text(size=12, face="bold", vjust=1),
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1)
    )
  
  svglite(out_path, height = 6, width = 6)
  print(ggp)
  dev.off()
  
  return(ggp)
}

#create pathway function
do_pathway = function(gene_input, em_scaled, sample_groups, master, fold_change, out_dir)
{
  #run ORA  
  sig_genes_entrez = bitr(gene_input, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  
  #make plots from ORA results
  ggp1 = barplot(ora_results, showCategory = 10)+
    scale_fill_gradient(high = "cadetblue3", low =  "coral2")

  svglite(file.path(out_dir, "ORA_barplot.svg"), height = 6, width = 6)
  print(ggp1)
  dev.off()
  
  #run GSEA
  gsea_input = fold_change
  names(gsea_input)=row.names(master)
  gsea_input=na.omit(gsea_input)
  gsea_input=sort(gsea_input, decreasing = TRUE)
  
  gse_results = gseGO(geneList = gsea_input,
                      ont = "MF",
                      keyType = "SYMBOL",
                      nPerm = 1000,
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      verbose = TRUE,
                      OrgDb = org.Hs.eg.db,
                      pAdjustMethod = "none")
  
  ggp2 = ridgeplot(gse_results)
  
  svglite(file.path(out_dir, "GSEA_ridgeplot.svg"), height = 6, width = 6)
  print(ggp2)
  dev.off()
  
  #extract genes from top 3 pathways
  gene_sets = ora_results$geneID
  description = ora_results$Description
  p.adj = ora_results$p.adjust
  
  ora_results_table = data.frame(cbind(gene_sets, p.adj))
  row.names(ora_results_table) = description
  
  enriched_gene_set = as.character(ora_results_table[1:3,1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  
  #plot heatmap - function within function?
  scaled_table = data.frame(t(scale(candidate_em)))
  ggp3 = plot_heatmap(scaled_table,  file.path(out_dir, "ORA_heatmap.svg"))
  
  #run STRING
  candidate_genes_table = data.frame(candidate_genes)
  names(candidate_genes_table) = "gene"
  
  string_db = STRINGdb$new(version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE)
  string = string_db$plot_network(string_mapped)
  
  svglite(file.path(out_dir, "STRING_network.svg"), height = 6, width = 6)
  print(string)
  dev.off()
  
  #return all 4 plots
  return(list(ora_barplot = ggp1,
              gsea_plot = ggp2,
              heatmap = ggp3,
              string_plot = string))
}
