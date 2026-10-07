#user note - run Functions.R to load functions into global environment before running this script!

#load packages
library(ggplot2)
library(ggrepel)
library(reshape2)
library(amap)
library(eulerr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(devEMF)
library(svglite)

#load in data tables
annotations = read.table("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/MT_KO/ANNO.csv", header=TRUE, row.names = 1, sep = "\t")
de_mtko_vs_senes = read.table("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/MT_KO/DE_Senes_MTKO_vs_Senes.csv", header=TRUE, row.names = 1, sep = "\t")
de_senes_vs_prolif = read.table("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/MT_KO/DE_Senes_vs_Prolif.csv", header=TRUE, row.names = 1, sep = "\t")
em = read.table("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/MT_KO/EM.csv", header=TRUE, row.names = 1, sep = "\t")
samples = read.table("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/MT_KO/SS.csv", header=TRUE, sep = "\t")

#merge 2x de tables
master_temp1 = merge(de_mtko_vs_senes, de_senes_vs_prolif,by.x = 0, by.y = 0, suffixes=c(".mvs",".svp"))

#merge with expression matrix
master_temp2 = merge(master_temp1, em, by.x = 1, by.y = 0)

#merge with annotations
master = merge(master_temp2, annotations, by.x = 1, by.y = 0)

#change row names and column 1 name
row.names(master)=master[,17]
master = master[,-17]
colnames(master)[1] = "Gene_ID"

#the na.omit() problem - i lose majority of genes with na.omit()
#just do it and justify
master = na.omit(master)

#add new columns - will be useful for subsetting and plotting certain graphs
#mean expression 
master$"Mean_Expression" = rowMeans(master[,8:16])

#-log10p
master$"mlog10p.mvs" = -log10(master$p.adj.mvs)
master$"mlog10p.svp" = -log10(master$p.adj.svp)

#column flagging significance
master$"sig.mvs" = as.factor(master$p.adj.mvs < 0.05 & abs(master$log2fold.mvs) > 1.0)
master$"sig.svp" = as.factor(master$p.adj.svp < 0.05 & abs(master$log2fold.svp) > 1.0)

#make scaled em
#first make em_symbols
em_symbols = master[,8:16]
#then make em_scaled
em_scaled = data.frame(t(scale(em_symbols)))
em_scaled = na.omit(em_scaled)

#make list of significant genes
#table of genes significant in EITHER fibroblast type
master_sig = subset(master, master$sig.mvs==TRUE | master$sig.svp==TRUE)
#take list of names
sig_genes = row.names(master_sig)

#make scaled expression values for sig genes
em_symbols_sig = master_sig[,8:16]
em_scaled_sig = data.frame(t(scale(em_symbols_sig)))
em_scaled_sig = na.omit(em_scaled_sig)

#make two master tables sorted by p-value
sorted_p_mvs = order(master[,"p.adj.mvs"], decreasing = FALSE)
sorted_p_svp = order(master[,"p.adj.svp"], decreasing = FALSE)
master_p_mvs = master[sorted_p_mvs,]
master_p_svp = master[sorted_p_svp,]

#more subsetting for volcano plots
#find most significantly upregulated and downregulated genes
#one volcano plot = one DE table
sig_up_mvs = subset(master_p_mvs, log2fold.mvs > 1.0 & sig.mvs == TRUE)
sig_down_mvs = subset(master_p_mvs, log2fold.mvs < -1.0 & sig.mvs == TRUE)
sig_up_svp = subset(master_p_svp, log2fold.svp > 1.0 & sig.svp == TRUE)
sig_down_svp = subset(master_p_svp, log2fold.svp < -1.0 & sig.svp == TRUE)
#top 6 totally arbitrary - only chose it because it looks nicer in faceted plots than 5 does
#top6mvs_up = sig_up_mvs[1:6,]
#top6mvs_down = sig_down_mvs[1:6,]
#top6svp_up = sig_up_svp[1:6,]
#top6svp_down = sig_down_svp[1:6,]
#^^^ THIS DOES NOT WORK - i think na.omit() on different tables caused loss of different genes in different tables, which had big impact in plotting boxplots later on
#need to find way around this - ensure top 6 genes are DEFINITELY IN em_scaled
mvs_up = sig_up_mvs[rownames(sig_up_mvs) %in% colnames(em_scaled), ]
mvs_down = sig_down_mvs[rownames(sig_down_mvs) %in% colnames(em_scaled), ]
svp_up = sig_up_svp[rownames(sig_up_svp) %in% colnames(em_scaled), ]
svp_down = sig_down_svp[rownames(sig_down_svp) %in% colnames(em_scaled), ]

top6mvs_up = mvs_up[1:6,]
top6mvs_down = mvs_down[1:6,]
top6svp_up = svp_up[1:6,]
top6svp_down = svp_down[1:6,]


#volcano plot with top 6 most sig genes in MvS annotated
ggp = ggplot(master, aes(x=log2fold.mvs, y=mlog10p.mvs))+
  geom_point(aes(colour = "a"), size = 1)+
  geom_point(data = sig_up_mvs, aes(colour = "b"), size = 1)+
  geom_point(data = sig_down_mvs, aes(colour = "c"), size = 1)+
  labs(title = "MTKO vs Senescent", x="Log2 Fold Change", y="-Log10p")+
  geom_vline(xintercept = -1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey", size = 0.5)+
  geom_text_repel(data=top6mvs_up, aes(label=row.names(top6mvs_up)), size = 2)+
  geom_text_repel(data=top6mvs_down, aes(label=row.names(top6mvs_down)), size = 2)+
  xlim(-12, 12)+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    panel.background = element_blank(),
    panel_border = element_rect(colour="black", fill=NA, size=1.5),
    panel.spacing=unit(1,"lines")
  )+
  scale_colour_manual(values = c("black", "coral2", "cadetblue3"), labels=c("No Change", "Upregulated", "Downregulated"), name="")
ggp

#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/volcano_mvs.svg", height = 6, width = 6)
print(ggp)
dev.off()

#same volcano plot, but with top 6 most sig SvP annotated
ggp = ggplot(master, aes(x=log2fold.svp, y=mlog10p.svp))+
  geom_point(aes(colour = "a"), size = 1)+
  geom_point(data = sig_up_svp, aes(colour = "b"), size = 1)+
  geom_point(data = sig_down_svp, aes(colour = "c"), size = 1)+
  labs(title = "Senescent vs Proliferative", x="Log2 Fold Change", y="-Log10p")+
  geom_vline(xintercept = -1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey", size = 0.5)+
  geom_text_repel(data=top6svp_up, aes(label=row.names(top6svp_up)), size = 2)+
  geom_text_repel(data=top6svp_down, aes(label=row.names(top6svp_down)), size = 2)+
  xlim(-12, 12)+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    panel.background = element_blank(),
    panel_border = element_rect(colour="black", fill=NA, size=1.5),
    panel.spacing=unit(1,"lines")
  )+
  scale_colour_manual(values = c("black", "coral2", "cadetblue3"), labels=c("No Change", "Upregulated", "Downregulated"), name="")
ggp

#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/volcano_svp.svg", height = 6, width = 6)
print(ggp)
dev.off()

#MA plot for MvS
ggp = ggplot(master, aes(x=log10(Mean_Expression), y=log2fold.mvs))+
  geom_point(aes(colour = "a"), size = 1)+
  geom_point(data = sig_up_mvs, aes(colour = "b"), size = 1)+
  geom_point(data = sig_down_mvs, aes(colour = "c"), size = 1)+
  labs(title = "MTKO vs Senescent", x="Log10(Mean Expression)", y="Log2(Fold Change)")+
  geom_hline(yintercept = -1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5)+
  ylim(-10, 10)+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    panel.background = element_blank(),
    panel_border = element_rect(colour="black", fill=NA, size=1.5),
    panel.spacing=unit(1,"lines")
  )+
  scale_colour_manual(values = c("black", "coral2", "cadetblue3"), labels=c("No Change", "Upregulated", "Downregulated"), name="")
ggp
#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/MA_mvs.svg", height = 6, width = 6)
print(ggp)
dev.off()


#MA plot for SvP
ggp = ggplot(master, aes(x=log10(Mean_Expression), y=log2fold.svp))+
  geom_point(aes(colour = "a"), size = 1)+
  geom_point(data = sig_up_svp, aes(colour = "b"), size = 1)+
  geom_point(data = sig_down_svp, aes(colour = "c"), size = 1)+
  labs(title = "Senescent vs Proliferative", x="Log10(Mean Expression)", y="Log2(Fold Change)")+
  geom_hline(yintercept = -1, linetype = "dashed", colour = "grey", size = 0.5)+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey", size = 0.5)+
  ylim(-10, 10)+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    panel.background = element_blank(),
    panel_border = element_rect(colour="black", fill=NA, size=1.5),
    panel.spacing=unit(1,"lines")
  )+
  scale_colour_manual(values = c("black", "coral2", "cadetblue3"), labels=c("No Change", "Upregulated", "Downregulated"), name="")
ggp

#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/MA_svp.svg", height = 6, width = 6)
print(ggp)
dev.off()

#use created function to plot faceted density plots
plot_expression_density(em, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/expression_density.svg")


#quality of data has now been assessed with MA and expression density plots
#PCA next
#do fibroblast sample cluster or separate by condition?
#is MTKO acting like senescent or proliferative?
pca = prcomp(as.matrix(sapply(em_scaled, as.numeric)))
pca_coordinates = data.frame(pca$x)

#could i make a function for easy PCA plotting?
#would save me hard coding 
#add sample group to pca_coordinates to help function run properly
pca_coordinates$SAMPLE_GROUP = samples$SAMPLE_GROUP

#function was run in control panel mainly - hard coding PCA plot with highest clustering and most notable separation
plot_pca(pca, pca_coordinates, "PC1", "PC2", "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/pca.svg")
#PC1 is ~88%... 88% of variability explained by PC1

#overlap analysis of genes in MvS and SvP
#collect list of significant genes for each DE
sig_mvs = row.names(subset(master, sig.mvs == TRUE))
sig_svp = row.names(subset(master, sig.svp == TRUE))

#vectors in list must be named - save to object Venn_data
venn_data = list("MTKO v Senes"=sig_mvs,"Senes v Prolif"=sig_svp)

#subset master to find genes that are in only one DE
sig_mvs_only = row.names(subset(master,master$sig.mvs==TRUE & master$sig.svp==FALSE))
sig_svp_only = row.names(subset(master,sig.mvs==FALSE & sig.svp==TRUE))
sig_both = row.names(subset(master,sig.mvs==TRUE & sig.svp==TRUE))

#hypergeometric test to find overlap p-value
group1 = nrow(subset(master, master$sig.mvs==TRUE))
group2 = nrow(subset(master, master$sig.svp==TRUE))
overlap = nrow(subset(master, master$sig.mvs==TRUE & master$sig.svp ==TRUE))
total = nrow(master)
overlap.p = phyper(overlap-1, group2, total-group2, group1, lower.tail=FALSE)

#plot venn diagram
#find out how to change colours of V-D
fit = euler(venn_data)
ggp = plot(fit, fills=c("coral2", "bisque2"), shape="ellipse", edges=FALSE, quantities = TRUE)
#overlap of 188 genes is significant (overlap.p = 1.18e-64) - include in figure legend when presenting in report
#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/overlap.svg", height=6, width=6)
print(ggp)
dev.off()

#fold vs fold plot
ggp = ggplot(master, aes(x = log2fold.mvs, y = log2fold.svp))+
  geom_point()+
  labs(x="Log2(Fold Change in MTKO vs Senes)", y="Log2(Fold Change in Senes vs Prolif)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA, linewidth=1),
        axis.title = element_text(size=14),
        axis.text = element_text(size=10))
ggp

#export
svglite("C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/fold_vs_fold.svg", height=6, width=6)
print(ggp)
dev.off()

#now have clear overlap - 188 genes significant in both MvS and SvP
#heatmaps next
#use em_scaled_sig- too many genes
#can't interpret anything from just a colour block
#heatmap of overlapping genes? 
#i have list of overlap gene names
#can obtain gene info from master table
candidate_genes = master[sig_both,]
candidate_em = candidate_genes[,8:16]
candidate_em_scaled = data.frame(t(scale(candidate_em)))

#plot heatmap with created function
plot_heatmap(candidate_em_scaled, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/overlapping_genes_heatmap.svg")

#replot with top 50 most significant genes in both MvS and SvP
mvs_top_50 = data.frame(t(scale(master_p_mvs[1:50,8:16])))
plot_heatmap(mvs_top_50, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/mvs_top_50_heatmap.svg")

svp_top_50 = data.frame(t(scale(master_p_svp[1:50,8:16])))
plot_heatmap(svp_top_50, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/svp_top_50_heatmap.svg")

#identifiable gene of interest from heatmap for boxplot - FTH1
gene_data = em_symbols["FTH1",]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = samples$SAMPLE_GROUP
names(gene_data) = c("expression", "sample_group")

#FTH1 boxplot
plot_single_boxplot("FTH1", em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/FTH1_boxplot.svg")
#boxplot shows low expression in healthy prolifilative cells, but greatly upregulated in senescent and senes_MTKO
#likely makes it a cell-cycle arrest gene
#research shows FTH1 involved in iron-induced cellular damage and ferroptosis

#another one is GAPDH
plot_single_boxplot("GAPDH", em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/GADPH_boxplot.svg")
#highly upregulated in MTKO compared to Senes... however prolif expression is only slightly upregulated
#this gene probably not as important for restorative function
#research anyway to be sure

rownames(top6mvs_down)
setdiff(rownames(top6mvs_down), colnames(em_scaled))

#multi boxplots
#look at top 5 most up/downregulated in both DE tables
#already identified and labelled from volcano plot
#modify single boxplot function to make faceted multi-gene boxplot function
#make list of gene names
#actually, i can include the list making process in function for less hard-coding
plot_multi_boxplot(top6mvs_up, em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/mvs_top_6_up_boxplot.svg")
plot_multi_boxplot(top6mvs_down, em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/mvs_top_6_down_boxplot.svg")
plot_multi_boxplot(top6svp_up, em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/svp_top_6_up_boxplot.svg")
plot_multi_boxplot(top6svp_down, em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/svp_top_6_down_boxplot.svg")

#ORA
#convert IDs
sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#15.49% of genes failed to map - big data loss

#run ORA
#look at ALL first - then would probably look at BP to define what biological processes exactly are implicated in both cell states
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#set colour palette
colours = c("cadetblue3", "coral2")
palette = colorRampPalette(colours)(50)

#barplot
ggp = barplot(ora_results, showCategory = 10)+
  scale_fill_gradient(high = "cadetblue3", low =  "coral2")
ggp

#dotplot
ggp = dotplot(ora_results, showCategory = 10)+
  scale_fill_gradient(high = "cadetblue3", low =  "coral2")
ggp

#GO plot
#ggp = goplot(ora_results, showCategory = 10)+
  #scale_fill_gradient(high = "cadetblue3", low =  "coral2")
#ggp

#cnet plot
#ggp = cnetplot(ora_results, showCategory = 10)+
  #scale_fill_gradient(high = "cadetblue3", low =  "coral2")
#ggp

#didn't work, couldn't connect to GOALL.sqlite.gz so commented them out
#pathway of overlapping genes
#give insight into biology of genes significant in both cell DE 
#16.5% of genes lost
sig_genes_entrez = bitr(sig_both, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
ggp = barplot(ora_results, showCategory = 10)+
  scale_fill_gradient(high = "cadetblue3", low =  "coral2")
ggp

#pathway of genes significant in MvS and SvP
#significantly upregulated and downregulated
#use function
#limitation of code - i can't figure out how to code to any user's cwd. I only know how to hard code to my own folders. following lines will not work properly on another person's computer :/
gene_input = row.names(sig_up_mvs)
do_pathway(gene_input, em_scaled, samples$SAMPLE_GROUP, master, master$log2fold.mvs, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/Sig_Up_MvS_Path")

gene_input = row.names(sig_down_mvs)
do_pathway(gene_input, em_scaled, samples$SAMPLE_GROUP, master, master$log2fold.mvs, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/Sig_Down_MvS_Path")

gene_input = row.names(sig_up_svp)
do_pathway(gene_input, em_scaled, samples$SAMPLE_GROUP, master, master$log2fold.svp, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/Sig_Up_SvP_Path")

gene_input = row.names(sig_down_svp)
do_pathway(gene_input, em_scaled, samples$SAMPLE_GROUP, master, master$log2fold.svp, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/Sig_Down_SvP_Path")

do_pathway(sig_both, em_scaled, samples$SAMPLE_GROUP, master, master$log2fold.svp, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/Sig_Both_Path")

#found a better gene to analyse individually
plot_single_boxplot("TAGLN2", em_scaled, samples$SAMPLE_GROUP, "C:/UofG/Year 5/Data Exploration and Interpretation for Bioinformatics/Assessments/Figures/TAGLN2.svg")

write.csv(master, file="master.csv", quote = FALSE)

