# library(monocle); ; library(ggplot2); library(cowplot); library(data.table);
# library(GSA); library(GSVA); library(ggpubr); library(ggthemes); library(car); library(qvalue);
# library(data.table); library(limma); library(sva); library(randomForest)
# library(FactoMineR);library(factoextra); library(survival); library(Seurat); library(biomaRt);
# library(survival); library(survminer); library(maftools);
# library(gplots); library(igraph); library(irlba); library(Rtsne); library(densityClust);
# library(plyr); library(dplyr); library(caret); library(broom); library(impute);
# library(genefu)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

library(cellrangerRkit)
genome <- "GRCh38"
prad = load_cellranger_matrix("duke_aggr/",genome=genome)
dim(exprs(prad))

library(plyr); library(car)
pradTMP = exprs(prad)
ind = laply(1:ncol(pradTMP),function(i) strsplit(colnames(pradTMP)[i],split="-")[[1]][2],.progress="time")
ind = as.numeric(ind)
pradTMP1 = data.frame("Barcodes"=colnames(pradTMP),"Index"=ind)
pradTMP1$PatientID = ifelse(pradTMP1$Index==1|pradTMP1$Index==2,"Patient_1","Patient_2")
pradTMP1$Cell_Type = as.factor(pradTMP1$Index)
pradTMP1$Cell_Type = car::recode(pradTMP1$Cell_Type,"'1'='Benign_Basal';'2'='Cancer_Luminal';'3'='Benign_Basal';'4'='Benign_Luminal';'5'='Benign_NE';'6'='Cancer_Luminal';'7'='Cancer_NE';")
rownames(pradTMP1) = pradTMP1$Barcodes
pradTMP2 = pradTMP1[pradTMP1$PatientID=="Patient_2",]
dim(pradTMP2)
pradEXP = prad[,match(pradTMP2$Barcodes,colnames(prad))]

library(data.table)
probes = fread("~/Research/pathways/Homo_sapiens.GRCh38.79_table.txt",header=F,data.table=F)
colnames(probes)=c("ENSEMBL","Genes","Coordinates","Strand")
probes$Chr = gsub('\\:.*', '', probes$Coordinates)
probes=probes[c(which(probes$Chr %in% c(1:22)),grep("^X$",probes$Chr),grep("^Y$",probes$Chr)),]
probes = probes[-grep("^MIR",probes$Genes),]
probes = probes[-grep("^RPS",probes$Genes),]
probes = probes[-grep("^RPL",probes$Genes),]

keep_probes = intersect(rownames(pradEXP),probes$ENSEMBL)
seu_exprs = pradEXP[match(keep_probes,rownames(pradEXP)),]
probes_exprs = probes[match(keep_probes,probes$ENSEMBL),]

seu_exprs = as.data.frame(as.matrix(exprs(seu_exprs)))
seu_exprs$Genes = as.factor(as.character(probes_exprs$Genes))
dat = setDT(seu_exprs)[, lapply(.SD, median), by = Genes]
dat = as.data.frame(dat)
rownames(dat) = dat[,1]
dat = dat[,-1]
fwrite(dat,"/Users/ca31/Research/PRAD/single_cell/PRAD_scRNAseq_gene_level.txt",sep="\t",col.names=T,row.names=T,quote=F)

rm(pradTMP, pradTMP1, pradEXP, seu_exprs, probes_exprs)
gc()

###############
library(Seurat)
seed = 123
set.seed(seed)
seu = CreateSeuratObject(raw.data = dat, min.cells=3, min.genes=200,is.expr = 0.5, project = "PRAD", meta.data = pradTMP2)

pdf("SEURAT_PRAD_plots.pdf",width=11,height=6.71)
VlnPlot(object=seu,features.plot=c("nGene", "nUMI"),nCol=2)
GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")
seu = FilterCells(object = seu,subset.names = c("nUMI"), high.thresholds = c(2e7))
seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu = FindVariableGenes(seu, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
seu = ScaleData(object = seu, vars.to.regress = c("nUMI"))
seu = RunPCA(seu, pcs.print = 0,pc.genes = seu@var.genes)
PCAPlot(object = seu, dim.1 = 1, dim.2 = 2)
PCHeatmap(object = seu, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
# seu = JackStraw(object = seu, num.replicate = 100)
# JackStrawPlot(object = seu, PCs = 1:20)
PCElbowPlot(seu)
seu = FindClusters(object = seu, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = seu)
table(seu@ident)
set.seed(seed)
seu = RunTSNE(seu, dims.use = 1:10,do.fast=T)
TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
TSNEPlot(seu, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")

fwrite(data.frame(as.matrix(seu@data),check.names=F),"SEURAT_logNormalized_data.txt",sep="\t",col.names=T,row.names=T,quote=F)

ne_gene_list = c("CGA","CXCR2","SYN1","SYN2","SYN3","SYP","NCAM1")
lum_gene_list = c("AR","KLK3")
basal_gene_list = c("TP63","KRT5","KRT18")
new_list = c("FOXA1", "TMPRSS2","HOXB13","KRT8","NKX3-1","SOX2")
new_list1 = c("ABCF1","EYA3","DENND4C","MSTO1","FUBP3")

FeaturePlot(seu, features.plot = ne_gene_list, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = ne_gene_list,group.by="Cell_Type")
DotPlot(object = seu, genes.plot = ne_gene_list, plot.legend = TRUE)

FeaturePlot(seu, features.plot = lum_gene_list, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = lum_gene_list,group.by="Cell_Type")
DotPlot(object = seu, genes.plot = lum_gene_list, plot.legend = TRUE)

FeaturePlot(seu, features.plot = basal_gene_list, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = basal_gene_list,group.by="Cell_Type")
DotPlot(object = seu, genes.plot = basal_gene_list, plot.legend = TRUE)

FeaturePlot(seu, features.plot = new_list, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = new_list,group.by="Cell_Type")
DotPlot(object = seu, genes.plot = new_list, plot.legend = TRUE)

FeaturePlot(seu, features.plot = new_list1, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = new_list1,group.by="Cell_Type")
DotPlot(object = seu, genes.plot = new_list1, plot.legend = TRUE)

FeaturePlot(seu, features.plot = c("CXCL8"), min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(seu,features.plot = c("CXCL8"),group.by="Cell_Type")
DotPlot(object = seu, genes.plot = c("CXCL8"), plot.legend = TRUE)

# markers <- FindAllMarkers(object = seu,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
# top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# DoHeatmap(object = seu, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5)
clust_ave_cell = AverageExpression(object = seu)
clust_ave = AverageExpression(object = seu, return.seurat = TRUE, show.progress = FALSE)
DoHeatmap(object = clust_ave, genes.use = PCTopGenes(object = seu, pc.use = 1,do.balanced = TRUE), group.label.rot = TRUE, group.cex = 0)
dev.off()

new_list3 = c("AR","FOXA1","TP63","KLK3","TMPRSS2","KRT5","HOXB13","KRT8","KRT18","CXCL8")
FeaturePlot(seu, features.plot = new_list3, min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
ggsave("feature_plot_cheng.pdf")

## NE = CGA (ENSG00000135346), CXCR2 (ENSG00000180871), SYN1/SYN2/SYN3
## LUM = AR (ENSG00000169083), KLK3 (ENSG00000142515)
## BASAL = TP63 (ENSG00000073282)

pdf("SEURAT_PRAD_new_clusters_plots.pdf",width=11,height=6.71)
new_clusters = seu@ident
new_clusters = car::recode(new_clusters,"1=0;4=0")
new_clusters = car::recode(new_clusters,"6=2; 7=2;10=2")
new_clusters = car::recode(new_clusters,"11=9;12=9;13=9")
new_clusters = car::recode(new_clusters,"5=3;8=3")
new_cell_type=ifelse(new_clusters=="0", "Basal",ifelse(new_clusters=="2","LUM_1",ifelse(new_clusters=="3","LUM_2","NE")))
new_clusters = car::recode(new_clusters,"2=1;3=2;9=3")
seu@meta.data$new_clusters = new_clusters
seu@meta.data$new_cell_type = new_cell_type
table(new_clusters)

TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
TSNEPlot(seu, do.label = FALSE, pt.size = 0.5,group.by="new_cell_type")
clust_ave = AverageExpression(object = seu, return.seurat = TRUE, show.progress = FALSE)
DoHeatmap(object = clust_ave, genes.use = PCTopGenes(object = seu, pc.use = 1,do.balanced = TRUE), group.label.rot = TRUE, group.cex = 0)
markers <- FindAllMarkers(object = seu,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seu, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5)
dev.off()

### p53 geneset
library(GSVA)
p53 = scan("~/Research/PRAD/single_cell/p53_geneset.txt",what="")
p53_gsva = gsva(as.matrix(seu@data),gset.idx.list=list(p53),method="ssgsea",kcdf="Poisson",min.sz=1,max.sz=1000)
seu@meta.data$p53 = as.numeric(scale(as.numeric(p53_gsva)))
p53_agg = aggregate(p53 ~ res.1,data=seu@meta.data,FUN="mean")
p53_agg$res.1 = as.numeric(p53_agg$res.1)
p53_agg = p53_agg[order(p53_agg$p53),]
colnames(p53_agg)[1] = c("Cluster")
p53_agg$Cluster = as.factor(p53_agg$Cluster)
p53_agg$Cluster = factor(p53_agg$Cluster,levels=p53_agg$Cluster)
p53_agg$Type = ifelse(p53_agg$p53 < 0, "down_reg", "up_reg")
ggplot(p53_agg, aes(x=Cluster, y=p53, label=Type)) + geom_bar(stat='identity', aes(fill=Type), width=.5)  + scale_fill_manual(name="p53 pathway",
                    labels = c("Under expression", "Over expression"),
                    values = c("up_reg"="#00ba38", "down_reg"="#f8766d")) + coord_flip()
ggsave("p53_pathway.pdf")

#############################
### Extract specific clusters
## NE: 9, 11, 12 and 13
## Basal: 0, 1, and 4
##############################

seu1 = SubsetData(seu,ident.use = c("0","1","4","9","11","12","13"))
seu1@raw.data = seu1@raw.data[,match(colnames(seu1@data),colnames(seu1@raw.data))]
TSNEPlot(seu1, do.label = FALSE, pt.size = 0.5,group.by="new_cell_type")
TSNEPlot(seu1, do.label = FALSE, pt.size = 0.5,group.by="new_clusters")
seu1@raw.data = seu1@raw.data[,match(colnames(seu1@data),colnames(seu1@raw.data))]

# seu3 = SubsetData(seu,ident.use = c("0","1","4"))
# TSNEPlot(seu3, do.label = TRUE, pt.size = 0.5)
# seu3@raw.data = seu3@raw.data[,match(colnames(seu3@data),colnames(seu3@raw.data))]

NE_seu1 = SubsetData(seu1,cells.use= rownames(subset(seu1@meta.data,Cell_Type=="Cancer_NE")))
NE_seu1@raw.data = NE_seu1@raw.data[,match(colnames(NE_seu1@data),colnames(NE_seu1@raw.data))]
TSNEPlot(NE_seu1, do.label = TRUE, pt.size = 0.5)

FeaturePlot(NE_seu1, features.plot = c("CXCL8","CXCR2","TP53"), min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

###########
## Monocle
######################
library(monocle)

###
## Only NE cells
##
dataSet = NE_seu1

pdf("NE_cell_trajectory_maps.pdf",width=11,height=6.7)

gene_df = data.frame(gene_short_name=rownames(dataSet@data))
rownames(gene_df) = gene_df$gene_short_name
fd = new("AnnotatedDataFrame", data = gene_df)

mycdsNE = newCellDataSet(as.matrix(dataSet@raw.data),phenoData = new("AnnotatedDataFrame", data = dataSet@meta.data), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycdsNE = estimateSizeFactors(mycdsNE)
mycdsNE = estimateDispersions(mycdsNE)
dim(mycdsNE)
mycdsNE = detectGenes(mycdsNE, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycdsNE),num_cells_expressed >= 10))
mycdsNE = mycdsNE[expressed_genes,]

L <- log(exprs(mycdsNE))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")

set.seed(123)
diff_test_res <- differentialGeneTest(mycdsNE[expressed_genes,],fullModelFormulaStr = "~res.1")
diff_test_res = diff_test_res[order(diff_test_res$pval),]
diff_test_sig = subset(diff_test_res,diff_test_res$qval<=0.01)
write.table(diff_test_sig,"NE_diff_exp_cluster.txt",col.names=NA,quote=F,sep="\t")
ordering_genes <- row.names(diff_test_sig)
mycdsNE <- setOrderingFilter(mycdsNE, ordering_genes)
mycdsNE <- reduceDimension(mycdsNE, max_components = 2,method = 'DDRTree')
mycdsNE <- orderCells(mycdsNE)
TSNEPlot(dataSet, do.label = TRUE, pt.size = 0.5)
plot_cell_trajectory(mycdsNE, color_by = "res.1",cell_size = 2.5)
plot_cell_trajectory(mycdsNE, color_by = "Pseudotime",cell_size = 2.5)
plot_cell_trajectory(mycdsNE, color_by = "State",cell_size = 2.5)
write.table(pData(mycdsNE),"NE_sampleInfo.txt",col.names=NA,quote=F,sep="\t")

diff_test_res_pseudo <- differentialGeneTest(mycdsNE,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res_pseudo,"NE_diff_exp_pseudotime.txt",col.names=NA,quote=F,sep="\t")
sig_gene_names = subset(diff_test_res_pseudo, qval < 0.1)
sig_gene_names = rownames(sig_gene_names[order(sig_gene_names$pval),])[1:10]

plot.new()
plot_pseudotime_heatmap(mycdsNE[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
plot.new()
plot_pseudotime_heatmap(mycdsNE[ne_gene_list[ne_gene_list %in% rownames(exprs(mycdsNE))],],num_clusters = 2,cores = 1,show_rownames = T)
plot_genes_in_pseudotime(mycdsNE[ne_gene_list[ne_gene_list %in% rownames(exprs(mycdsNE))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycdsNE[basal_gene_list[basal_gene_list %in% rownames(exprs(mycdsNE))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycdsNE[lum_gene_list[lum_gene_list %in% rownames(exprs(mycdsNE))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycdsNE[new_list[new_list %in% rownames(exprs(mycdsNE))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycdsNE[new_list1[new_list1 %in% rownames(exprs(mycdsNE))],],color_by = "new_cell_type", ncol = 1)
dev.off()

####
# Basal + NE
####

pdf("NE+Basal_cell_trajectory_maps.pdf",width=11,height=6.7)

seu1 = SubsetData(seu,ident.use = c("0","1","4","9","11","12","13"))
seu1@raw.data = seu1@raw.data[,match(colnames(seu1@data),colnames(seu1@raw.data))]
plot1 = TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
plot1a = plot1 + theme(legend.position='none')
TSNEPlot(seu1, do.label = TRUE, pt.size = 0.5)
plot2 = TSNEPlot(seu1, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")
plot2a = plot2 + theme(legend.position='none')
plot_grid(plot1a, plot2a)

dataSet = seu1

gene_df = data.frame(gene_short_name=rownames(dataSet@data))
rownames(gene_df) = gene_df$gene_short_name
fd = new("AnnotatedDataFrame", data = gene_df)

mycds = newCellDataSet(as.matrix(dataSet@raw.data),phenoData = new("AnnotatedDataFrame", data = dataSet@meta.data), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycds = estimateSizeFactors(mycds)
mycds = estimateDispersions(mycds)
dim(mycds)
mycds = detectGenes(mycds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycds),num_cells_expressed >= 10))
mycds = mycds[expressed_genes,]

L <- log(exprs(mycds))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")

set.seed(123)
diff_test_res <- differentialGeneTest(mycds[expressed_genes,],fullModelFormulaStr = "~res.1")
diff_test_res = diff_test_res[order(diff_test_res$pval),]
diff_test_sig = subset(diff_test_res,diff_test_res$qval<=0.01)
write.table(diff_test_sig,"NE+Basal_diff_exp_cluster.txt",col.names=NA,quote=F,sep="\t")
ordering_genes <- row.names(diff_test_sig)[1:2000]
mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2,method = 'DDRTree')
mycds <- orderCells(mycds)
TSNEPlot(dataSet, do.label = TRUE, pt.size = 0.5)
plot_cell_trajectory(mycds, color_by = "res.1",cell_size = 2.5)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 2.5)
plot_cell_trajectory(mycds, color_by = "State",cell_size = 2.5)
write.table(pData(mycds),"NE+Basal_sampleInfo.txt",col.names=NA,quote=F,sep="\t")

diff_test_res_pseudo <- differentialGeneTest(mycds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res_pseudo,"NE+Basal_diff_exp_pseudotime.txt",col.names=NA,quote=F,sep="\t")
sig_gene_names = subset(diff_test_res_pseudo, qval < 0.1)
sig_gene_names = rownames(sig_gene_names[order(sig_gene_names$pval),])[1:10]

plot.new()
plot_pseudotime_heatmap(mycds[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
plot.new()
plot_pseudotime_heatmap(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],num_clusters = 2,cores = 1,show_rownames = T)
plot_genes_in_pseudotime(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[basal_gene_list[basal_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[lum_gene_list[lum_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list[new_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list1[new_list1 %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)

dev.off()

####
# Luminal A and B
################

pdf("luminal_trajectory_maps.pdf",width=11,height=6.7)

seu2 = SubsetData(seu,ident.use = c("2","3","5","6","7","8","10"),subset.raw=T)
TSNEPlot(seu2, do.label = TRUE, pt.size = 0.5)
TSNEPlot(seu2, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")
TSNEPlot(seu2, do.label = FALSE, pt.size = 0.5,group.by="new_cell_type")
TSNEPlot(seu2, do.label = FALSE, pt.size = 0.5,group.by="new_clusters")

plot1 = TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
plot1a = plot1 + theme(legend.position='none')
TSNEPlot(seu2, do.label = TRUE, pt.size = 0.5)
plot2 = TSNEPlot(seu2, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")
plot2a = plot2 + theme(legend.position='none')
plot_grid(plot1a, plot2a)

dataSet = seu2

gene_df = data.frame(gene_short_name=rownames(dataSet@data))
rownames(gene_df) = gene_df$gene_short_name
fd = new("AnnotatedDataFrame", data = gene_df)

mycds = newCellDataSet(as.matrix(dataSet@raw.data),phenoData = new("AnnotatedDataFrame", data = dataSet@meta.data), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycds = estimateSizeFactors(mycds)
mycds = estimateDispersions(mycds)
dim(mycds)
mycds = detectGenes(mycds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycds),num_cells_expressed >= 10))
mycds = mycds[expressed_genes,]

L <- log(exprs(mycds))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")

set.seed(123)
diff_test_res <- differentialGeneTest(mycds[expressed_genes,],fullModelFormulaStr = "~res.1")
diff_test_res = diff_test_res[order(diff_test_res$pval),]
diff_test_sig = subset(diff_test_res,diff_test_res$qval<=0.01)
write.table(diff_test_sig,"Luminal_diff_exp_cluster.txt",col.names=NA,quote=F,sep="\t")
ordering_genes <- row.names(diff_test_sig)[1:2000]
mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2,method = 'DDRTree')
mycds <- orderCells(mycds)

TSNEPlot(dataSet, do.label = TRUE, pt.size = 0.5)
plot_cell_trajectory(mycds, color_by = "res.1",cell_size = 1.5)
mycds = orderCells(mycds,root_state = 3)
plot_cell_trajectory(mycds, color_by = "res.1",cell_size = 1.5)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 2.5)
# plot_cell_trajectory(mycds, color_by = "State",cell_size = 2.5)
write.table(pData(mycds),"Luminal_sampleInfo.txt",col.names=NA,quote=F,sep="\t")

diff_test_res_pseudo <- differentialGeneTest(mycds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res_pseudo,"Luminal_diff_exp_pseudotime.txt",col.names=NA,quote=F,sep="\t")
sig_gene_names = subset(diff_test_res_pseudo, qval < 0.1)
sig_gene_names = rownames(sig_gene_names[order(sig_gene_names$pval),])[1:10]

plot.new()
plot_pseudotime_heatmap(mycds[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
plot.new()
plot_pseudotime_heatmap(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],num_clusters = 2,cores = 1,show_rownames = T)
plot_genes_in_pseudotime(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[basal_gene_list[basal_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[lum_gene_list[lum_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list[new_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list1[new_list1 %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)

dev.off()

####
# NE + Luminal A and B
################

pdf("NE+lum_trajectory_maps.pdf",width=11,height=6.7)

seu3 = SubsetData(seu,ident.use = c("2","3","5","6","7","8","9","10","11","12","13"),subset.raw=T)
TSNEPlot(seu3, do.label = TRUE, pt.size = 0.5)
TSNEPlot(seu3, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")
TSNEPlot(seu3, do.label = FALSE, pt.size = 0.5,group.by="new_cell_type")
TSNEPlot(seu3, do.label = FALSE, pt.size = 0.5,group.by="new_clusters")

plot1 = TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
plot1a = plot1 + theme(legend.position='none')
TSNEPlot(seu3, do.label = TRUE, pt.size = 0.5)
plot2 = TSNEPlot(seu3, do.label = FALSE, pt.size = 0.5,group.by="Cell_Type")
plot2a = plot2 + theme(legend.position='none')
plot_grid(plot1a, plot2a)

dataSet = seu3

gene_df = data.frame(gene_short_name=rownames(dataSet@data))
rownames(gene_df) = gene_df$gene_short_name
fd = new("AnnotatedDataFrame", data = gene_df)

mycds = newCellDataSet(as.matrix(dataSet@raw.data),phenoData = new("AnnotatedDataFrame", data = dataSet@meta.data), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycds = estimateSizeFactors(mycds)
mycds = estimateDispersions(mycds)
dim(mycds)
mycds = detectGenes(mycds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycds),num_cells_expressed >= 10))
mycds = mycds[expressed_genes,]

L <- log(exprs(mycds))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")

set.seed(123)
diff_test_res <- differentialGeneTest(mycds[expressed_genes,],fullModelFormulaStr = "~res.1")
diff_test_res = diff_test_res[order(diff_test_res$pval),]
diff_test_sig = subset(diff_test_res,diff_test_res$qval<=0.01)
write.table(diff_test_sig,"NE+lum_diff_exp_cluster.txt",col.names=NA,quote=F,sep="\t")
ordering_genes <- row.names(diff_test_sig)[1:5000]
mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2,method = 'DDRTree')
mycds <- orderCells(mycds)

TSNEPlot(dataSet, do.label = TRUE, pt.size = 0.5)
plot_cell_trajectory(mycds, color_by = "res.1",cell_size = 1.5)
mycds = orderCells(mycds,root_state = 3)
plot_cell_trajectory(mycds, color_by = "res.1",cell_size = 1.5)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 2.5)
# plot_cell_trajectory(mycds, color_by = "State",cell_size = 2.5)
write.table(pData(mycds),"NE+lum_sampleInfo.txt",col.names=NA,quote=F,sep="\t")

diff_test_res_pseudo <- differentialGeneTest(mycds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res_pseudo,"NE+lum_diff_exp_pseudotime.txt",col.names=NA,quote=F,sep="\t")
sig_gene_names = subset(diff_test_res_pseudo, qval < 0.1)
sig_gene_names = rownames(sig_gene_names[order(sig_gene_names$pval),])[1:10]

plot.new()
plot_pseudotime_heatmap(mycds[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
plot.new()
plot_pseudotime_heatmap(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],num_clusters = 2,cores = 1,show_rownames = T)
plot_genes_in_pseudotime(mycds[ne_gene_list[ne_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[basal_gene_list[basal_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[lum_gene_list[lum_gene_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list[new_list %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)
plot_genes_in_pseudotime(mycds[new_list1[new_list1 %in% rownames(exprs(mycds))],],color_by = "new_cell_type", ncol = 1)

dev.off()
