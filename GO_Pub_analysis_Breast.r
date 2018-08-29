#Pipeline for the Controls vs Breast DE analysis

defpath <- "/results/CvsB" #set the default path
usel <- "Breast" #Make Controls vs Breast analysis

##START OF ANALYSIS  --  ALL VS Controls ANALYSIS
dir.create(path=defpath, showWarnings = TRUE, recursive = FALSE, mode = "0777")

con <- file(paste(defpath, "Analysis_log.txt", sep="/"), "w")

cat(paste("Analysis started at:", Sys.time(),"\n", sep=" "), file=con) 
cat(paste("Default path of analysis:", defpath,"\n", sep =" "), file=con) 
cat(paste("Comparison selected:", "Controls vs", usel,"\n", sep =" "), file=con) 

print("--START OF THE ANALYSIS---")
cat("--START OF THE ANALYSIS--- \n", file=con) 
set.seed(12345)
##Loading of the libraries required
library(limma)
library(Rsubread)
library(edgeR)
library(DESeq2)
library("rentrez")
library(gplots)
library(topGO)
library(GO.db)
library(ggplot2)
cat("Loading required libraries: OK \n", file=con) 
print("Loading required libraries: OK", file=con) 

###Loading of the count tables
fc <- readRDS("../../output/featurecounts/counts.rds")
fc2 <- readRDS("../../output/featurecounts/counts.rds")
fc_17_05_2018 <- readRDS("../../output/featurecounts/counts_new_run_17_05_2018.rds")


###Getting the targets table
targets2606 <- read.csv("../../input/experiment_data/targets26-06.txt", sep="\t")
rownames(targets2606) <- targets2606$unique_id


###Old and new runs are being merged
conc_temp <- cbind(fc$counts, fc_17_05_2018$counts)
conc_temp <- cbind(conc_temp[,225:229], conc_temp[,1:224])
colnames(conc_temp) <- targets2606$unique_id


###Count matrix building
fc$counts <- conc_temp
fc$targets <- targets2606


###Get the Mapping coverage for each top-up from the featurecounts output file
mapped_coverage <- read.csv("../../output/featurecounts/detailed_files/dio_test.csv", sep=" ", header = FALSE, row.names=1)
mapped_coverage$transcriptome_bp <- rep_len(2161600, length.out=length(rownames(mapped_coverage)))
mapped_coverage$reads_mean_length <- rep_len(103.1454545, length.out=length(rownames(mapped_coverage)))
mapped_coverage$mapped_reads_coverage <- (mapped_coverage[,1]*mapped_coverage[,3])/mapped_coverage[,2]
mapped_coverage[mapped_coverage[,4]<25,]

##Get the IDs of Samples with less than 22X mapped reads coverage
excluded_ids <- rownames(mapped_coverage[mapped_coverage[,4]<25,])

##Get the IDs of the outliers observed in PCA plot: sample WN00067 - CEPH
#excCEPH <- rownames(targets2606[targets2606$sample =="CEPH",])
#excWN67 <- rownames(targets2606[targets2606$sample =="WN00067",])

##Get the IDs of cancers which are not Breast/Cervix/Skin
excnocancer <- rownames(targets2606[targets2606$Tissue %in% c("Exclude"),])

##Exclude the above top-ups
excluded_ids_all <- unique(c(excluded_ids, excnocancer))


###Create the final matrix and the final targets file
fc_filt <- fc$counts[, !colnames(fc$counts) %in% excluded_ids_all]
targets_filt <- targets2606[!rownames(targets2606) %in% excluded_ids_all,]
#targets_filt$phenotype2 <- ifelse(targets_filt$Tissue=="Control", "Control", "Case")

##Create the factor for the final comparison
Treat <- factor(targets_filt$Tissue, levels=c("Control","Brain", "Oesophagus", "Vulva", "Breast","Cervix","Skin"))
group <- Treat
group <- factor(group)
print("Importing the Data: OK")
cat("Importing the Data: OK \n", file=con) 

###Create the DGE object for EdgeR from the count matrix and the targets
y <- DGEList(fc_filt, group=group)
y <- y[, rownames(targets_filt)]
y$genes <-fc$annotation[,"Length", drop=FALSE] #Annotation
y$genes$Symbol <- fc$annotation$GeneID #Genes Symbol

##Sum the top-ups counts
y <- sumTechReps(y, ID=targets_filt$sample) #Collapse top-ups
Treat <- factor(y$samples$group, levels=c("Control","Brain", "Oesophagus", 
				"Vulva", "Breast","Cervix","Skin")) #Create the phenotype factor
y_cpm2 <- cpm(y)/y$genes$Length #Calculate CPM' values
print("Import the data to EdgeR: OK")
cat("Import the data to EdgeR: OK \n", file=con) 
colors_ins <- c("#4885ed", "#db3236", "cyan", "purple", "#f4c20d", "#3cba54", "pink") #Color palette
points_ins <- c(0,2,2,2,2,2,2,2,2,2,2,2,2) #Points palette

#PCA plot for CPM' values before excluding outliers
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_with_outliers_normal", save=TRUE)
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_with_outliers_small", save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=2, ltyl=0, lwdl=3)
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_with_outliers_small2", save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=1.5, ltyl=0, lwdl=3)
print("PCA plot on CPM values before filtering: OK")
cat("PCA plot on CPM values before filtering: OK \n", file=con) 
dir.create(path=paste(defpath, "/R_objects/", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

###Exclude the outliers WN00067 and CEPH
y <- y[,colnames(y)[!colnames(y) %in%c("WN00067","CEPH")], keep.lib.sizes=FALSE] #Exclude outliers

saveRDS(y, paste(defpath, "/R_objects/y.rds", sep="")) #Save the y object

Treat <- factor(y$samples$group, levels=c("Control","Brain", "Oesophagus", 
				"Vulva", "Breast","Cervix","Skin")) #Create the phenotype factor
y_cpm2 <- cpm(y)/y$genes$Length #Calculate CPM' values

#PCA plot for CPM' values after excluding outliers
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_no_outliers_normal", save=TRUE)
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_no_outliers_small", save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=2, ltyl=0, lwdl=3)
pca_plot(y_cpm2, phenotype=Treat, name="PCA_all_groups_no_outliers_small2", save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=1.5, ltyl=0, lwdl=3)
print("PCA plot on CPM values after excluding outliers: OK")
cat("PCA plot on CPM values after excluding outliers: OK \n", file=con) 
#Select comparison
y_all <- specify_cancer(usel=usel)
print(paste("Comparison Selected: Controls vs.", usel, sep=" ")) 

##Find which ids are cases or controls
y <- y_all
print("Samples phenotype QC:")
print(table(y$samples$group))
cat("Samples phenotype QC: \n", file=con) 
cat(table(y$samples$group), file=con) 
cat("\n", file=con) 
control_ids <- rownames(y$samples[y$samples$group==c("Control"),])
case_ids <- rownames(y$samples[y$samples$group==c("Case"),])

##Get a matrix for controls and a matrix for cases
ysum_controls <- y$counts[ , -which(colnames(y$counts) %in% case_ids)]
ysum_case <- y$counts[ , -which(colnames(y$counts) %in% control_ids)]

##Names of all the genes used in the analysis
all_genes_detected <- rownames(y)
saveRDS(all_genes_detected, paste(defpath, "/R_objects/all_genes_detected.rds", sep=""))


##Filter out Genes that had either controls or cases with <20 median coverage 
keep20 <- rowMedians(ysum_controls) > 20 | rowMedians(ysum_case) > 20
y <- y[keep20, , keep.lib.sizes=FALSE]
Treat <- factor(y$samples$group, levels=c("Control","Case")) #Create the phenotype factor
Treat_old <- factor(y$samples$group_old, levels=unique(y$samples$group_old)) #Create the old phenotype factor
print(paste(table(keep20)[1], "genes were filtered out.", sep=" "))
print(paste(table(keep20)[2], "genes remaining in the analysis.", sep=" "))
print("Filtering out genes: OK")
cat("Filtering out genes: OK \n", file=con) 
cat(paste(table(keep20)[1], "genes were filtered out. \n", sep=" "), file=con) 
cat(paste(table(keep20)[2], "genes remaining in the analysis. \n", sep=" "), file=con) 
cat("Filtering out genes: OK \n", file=con) 


y_cpm <- cpm(y) #Calculate CPM' values
y_cpm2 <- cpm(y)/y$genes$Length #Calculate CPM' values

#PCA plot for CPM' values after filtering out genes
pca_plot(y_cpm2, phenotype=Treat, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_normal", sep=""), save=TRUE)
pca_plot(y_cpm2, phenotype=Treat, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_small", sep=""), save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=2, ltyl=0, lwdl=3)
pca_plot(y_cpm2, phenotype=Treat, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_small2", sep=""), save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=1.4, ltyl=0, lwdl=3)
pca_plot(y_cpm2, phenotype=Treat_old, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_all_phen_normal", sep=""), save=TRUE)
pca_plot(y_cpm2, phenotype=Treat_old, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_all_phen_small", sep=""), save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=2, ltyl=0, lwdl=3)
pca_plot(y_cpm2, phenotype=Treat_old, name=paste("PCA_Controls_vs_", usel, "_filtered_genes_all_phen_small2", sep=""), save=TRUE, cmains=3, caxis=2, lwdp=3, cleg=1.4, ltyl=0, lwdl=3)
print("PCA plot after filtering genes out: OK")
cat("PCA plot after filtering genes out: OK \n", file=con) 
###Create a finaly filter back-up of y
finaly_filt <- y

###Get a matrix for controls and a matrix for cases
control_ids <- rownames(finaly_filt$samples[finaly_filt$samples$group==c("Control"),])
case_ids <- rownames(finaly_filt$samples[finaly_filt$samples$group==c("Case"),])

y_ctrl <- finaly_filt[, control_ids, keep.lib.sizes=FALSE]
y_cases <- finaly_filt[, case_ids, keep.lib.sizes=FALSE]
y_cpms <- cpm(y)

###Design Matrix
finaly_filt$samples$group <- factor(finaly_filt$samples$group, levels=c("Control", "Case"))
Treat <- factor(finaly_filt$samples$group, levels=c("Control","Case")) #Create the phenotype factor
design=model.matrix(~finaly_filt$samples$group)
saveRDS(finaly_filt, paste(defpath, "/R_objects/finaly_filt.rds", sep=""))

###Deseq2 differential expression analysis
dds <- DESeqDataSetFromMatrix(countData=finaly_filt$counts, colData=finaly_filt$samples, design= ~group)
featureData <- data.frame(gene=rownames(finaly_filt$counts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
print("Importing data to DESeq2: OK")
cat("Importing data to DESeq2: OK \n", file=con) 
##Wald & LRT tests

dds_wald = DESeq(dds)
dds_lrt = DESeq(dds, test="LRT", reduced=~1)
saveRDS(dds_wald, paste(defpath, "/R_objects/dds_wald.rds", sep=""))
saveRDS(dds_lrt, paste(defpath, "/R_objects/dds_lrt.rds", sep=""))

res_wald <- results(dds_wald)
res_lrt <- results(dds_lrt)
saveRDS(res_wald, paste(defpath, "/R_objects/res_wald.rds", sep=""))
saveRDS(res_lrt, paste(defpath, "/R_objects/res_lrt.rds", sep=""))


###Get the final normalised gene expression values
deseq_norm <- counts(dds_wald, normalized=TRUE)
saveRDS(deseq_norm, paste(defpath, "/R_objects/deseq_norm.rds", sep=""))


##Create the objects of the final result
res_dewald <- do_test_deseq(res_wald, dds_wald)
res_delrt <- do_test_deseq(res_lrt, dds_lrt)

##Save the ocjects and write them to CSV files
saveRDS(res_dewald, paste(defpath, "/R_objects/res_dewald.rds", sep=""))
saveRDS(res_delrt, paste(defpath, "/R_objects/res_delrt.rds", sep=""))
write.csv(res_dewald, paste(defpath,"deseqres_wald.csv",sep="/"))
write.csv(res_delrt, paste(defpath,"deseqres_lrt.csv",sep="/"))

deseq_plots(defpath, usel)
print("DESeq2 DE analysis: OK")
cat("DESeq2 DE analysis: OK \n", file=con) 


###Extract important genes from the DE analysis
#Genes with pvalue<0.05 and abs(logFC)>0.5 in the wald test 

#imp_genes_fc15_p05_wald <- decide_important(res_dewald, abs_logFC_cutoff = 1.5, pval_cutoff=0.05)
#saveRDS(imp_genes_fc15_p05_wald,  paste(defpath, "/R_objects/imp_genes_fc15_p05_wald.rds", sep=""))

#imp_genes_fc17_wald <- decide_important(res_dewald, abs_logFC_cutoff = 1.78)
#saveRDS(imp_genes_fc17_wald,  paste(defpath, "/R_objects/imp_genes_fc17_wald.rds", sep=""))

#imp_genes_fc3_p05_wald <- decide_important(res_dewald, abs_logFC_cutoff = 1.278, pval_cutoff=0.05)
#saveRDS(imp_genes_fc3_p05_wald,  paste(defpath, "/R_objects/imp_genes_fc3_p05_wald.rds", sep=""))

imp_genes_fc10_p05_wald <- decide_important(res_dewald, abs_logFC_cutoff = 1, pval_cutoff=0.05)
saveRDS(imp_genes_fc10_p05_wald,  paste(defpath, "/R_objects/imp_genes_fc10_p05_wald.rds", sep=""))

imp_genes_fc12_p05_lrt <- decide_important(res_delrt, abs_logFC_cutoff = 1.2, pval_cutoff=0.05)
saveRDS(imp_genes_fc10_p05_wald,  paste(defpath, "/R_objects/imp_genes_fc10_p05_wald.rds", sep=""))

imp_genes_fc13_p05_wald <- decide_important(res_dewald, abs_logFC_cutoff = 1.3, pval_cutoff=0.05)
saveRDS(imp_genes_fc13_p05_wald,  paste(defpath, "/R_objects/imp_genes_fc13_p05_wald.rds", sep=""))
#Genes with pvalue<0.05 and abs(logFC)>0.5 in the lrt test 
					   
imp_genes_merge <- merge_gene_sets(a=list(imp_genes_fc10_p05_wald, imp_genes_fc12_p05_lrt,imp_genes_fc13_p05_wald)) #merging wald important_genes objects
saveRDS(imp_genes_merge,  paste(defpath, "/R_objects/imp_genes_merge.rds", sep=""))
print("Extracting significant gene-sets: OK")
cat("Extracting significant gene-sets: OK \n", file=con) 

##Sort the Deseq2 Normalised values - first Controls then Cases samples
values_sorted <- deseq_norm[,c(control_ids, case_ids)]

#Inspect the differences in the final deseq2 normalised values between cases and controls by creating
#barplot and a whisker's plot
dir.create(path=paste(defpath, "raw_genes",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
for(i in imp_genes_merge[[1]]){

	bargene(gene=i, path=paste(defpath, "raw_genes",sep="/"))

}
print("Saving gene-wise bar/box plots: OK")
cat("Saving gene-wise bar/box plots: OK \n", file=con) 

#Generate heat-maps for all different gene-sets
all_imp_genes <- ls(pattern = "imp_genes_")
hmpath <- paste(defpath, "/Heatmaps", sep="")
dir.create(path=paste(defpath, "/Heatmaps/", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

for(i in c(1:length(all_imp_genes))){
	
	#Heatmaps using euclidean and manhattan distances
	make_heatmap(get(all_imp_genes[i]), Treat_old, name=paste(all_imp_genes[i],"_manhattan_all",sep=""), 
	path=hmpath, save=TRUE, wi=25, h=15, dmethod="manhattan")
	
	make_heatmap(get(all_imp_genes[i]), Treat_old, name=paste(all_imp_genes[i],"_euclidean_all",sep=""),
	path=hmpath, save=TRUE, wi=25, h=15, dmethod="euclidean")
	
	make_heatmap(get(all_imp_genes[i]), Treat, name=paste(all_imp_genes[i],"_manhattan_CvC",sep=""), 
	path=hmpath, save=TRUE, wi=25, h=15, dmethod="manhattan")
	
	make_heatmap(get(all_imp_genes[i]), Treat, name=paste(all_imp_genes[i],"_euclidean_CvC",sep=""),
	path=hmpath, save=TRUE, wi=25, h=15, dmethod="euclidean")
	
	}
print("Creating heatmaps: OK")
cat("Creating heatmaps: OK \n", file=con) 

#GO pathway analysis
tgpath <- paste(defpath, "/GO_analysis/", sep="")

for(i in c(1:length(all_imp_genes))){

	assign(paste("topgo", all_imp_genes[i], sep="_"), go_analysis(get(all_imp_genes[i]), path=tgpath, name=all_imp_genes[i]))
  saveRDS(get(paste("topgo", all_imp_genes[i], sep="_")), file=paste(defpath, "/R_objects/", "topgo", "_", all_imp_genes[i], ".rds", sep=""))
}
print("TopGO analysis: OK")
cat("TopGO analysis: OK \n", file=con) 

#Interpretation of the GO pathway analysis - best table and figure	
all_topgo <- c(ls(pattern = "topgo_imp_genes"))

for(i in c(1:length(all_topgo))){
	
	go_barplot(get(all_topgo[i]), name=all_topgo[i], save=TRUE, path=paste(defpath, "GO_analysis/", sep="/"), test="weight01Fisher", 
		       mars=c(5.1,48,4.1,1.1), he=15, cen=1.5, cela=1.5, cele=2, legpos="right")
	a <- term_wise(get(all_topgo[i]), save=TRUE, path=paste(defpath, "GO_analysis/", sep="/"), name=all_topgo[i], ont="BP")
	b <- term_wise(get(all_topgo[i]), save=TRUE, path=paste(defpath, "GO_analysis/", sep="/"), name=all_topgo[i], ont="CC")
	c <- term_wise(get(all_topgo[i]), save=TRUE, path=paste(defpath, "GO_analysis/", sep="/"), name=all_topgo[i], ont="MF")
	
	}
 
go_barplot(get(all_topgo[1]), name=all_topgo[1], save=TRUE, path=paste(defpath, "GO_analysis/", sep="/"), test="weight01Fisher", 
	       mars=c(5.1,48,4.1,1.1), he=8, cen=1.5, cela=1.5, cele=2, legpos="bottomright")
print("GO Barplots and GO Term-wise: OK")
cat("GO Barplots and GO Term-wise: OK \n", file=con) 

pmpath <- paste(defpath, "Pubmed_search", sep="/")
dir.create(path=pmpath, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#NCBI automated literature search
for(i in c(1:length(all_imp_genes))){
	namez <- unlist(strsplit(all_imp_genes[i], "_"))
	namez <- paste(namez[!namez %in% c("imp", "genes")], collapse="_")
	assign(paste("pms", all_imp_genes[i], sep="_"), pubmed_search(get(all_imp_genes[i]), coterm=co_term, "Homo Sapiens", path=pmpath, name=namez))
	pubmed_barplot(pubmed=get(paste("pms", all_imp_genes[i], sep="_")), path=pmpath, save=TRUE, name=namez, wid=6.5)
}
print("Pubmed search and plots: OK")
cat("Pubmed search and plots: OK \n", file=con) 
close(con)