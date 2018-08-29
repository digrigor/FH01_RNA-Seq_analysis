load_variables <- function(path){
#Function that loads the stored variables from the differential 
#expression analysis given a specific path.

	finaly_filt <<- readRDS(paste(path, "/finaly_filt.rds", sep="")) #To get the phenotype data
	deseq_norm <<- readRDS(paste(path, "/deseq_norm.rds", sep="")) #Matrix of the final normalised values
	res_delrt <<- readRDS(paste(path, "/res_delrt.rds", sep="")) #LRT test outcome
	res_dewald <<- readRDS(paste(path, "/res_dewald.rds", sep="")) #Wald test outcome
	all_genes_detected <<- readRDS(paste(path, "/all_genes_detected.rds", sep="")) #all genes detected
}

###PCA
pca_plot <- function(ymatrix, phenotype, name, save=TRUE, legpos="bottomleft", hi, wi, cmains=1.5, caxis=1.2, lwdp=2, cleg=1.2, ltyl=0, lwdl=2, legtext){
	#Function that takes as its argument a gene expression count matrix with Samples (columns) 
	#and Genes (rows) and a phenotype vector with their phenotypes (in the same order as the samples),
	#performs a PCA and returns a PCA plot. 
	#Arguments:
		#ymatrix: gene expression count matrix. Columns: Samples, Rows: Genes.
		#phenotype: Factor. Phenotypes of the samples. The samples should be in the same
		#			order as in the matrix
		#save: Boolean. Whether or not the PCA plot should be saved. Default value is TRUE.
		#name: name of the output PNG image (if save=TRUE)
		#	   legpos: Keyword to be used to position the legend. Accepted keywords: 
		#      "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" 
		#      and "center"
		
	dir.create(path=paste(defpath, "/PCA/", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

	Xpca <- prcomp(t(ymatrix))
	Xscores <- Xpca$x
	s <- summary(Xpca)

	#Colors and points of the plot
	colors <- rep(colors_ins[1:length(levels(phenotype))], 2)
	points <- rep(points_ins[1:length(levels(phenotype))], 2)
	condition <- phenotype
  
  if(missing(legtext)){
		legtext <- c(levels(phenotype))
	}
  
	#Create the plot
	if(save==TRUE){
		png(paste(defpath, "/PCA/", name, ".png", sep=""), units="in", width=11.5, height=10, res=600)
	}

	par(mar=c(5.1,6,4.1,2.1))
	plot(Xscores, xlab=paste("Component 1 (", round(s$importance[2,1]*100, 2),"%)",sep=""),
		ylab=paste("Component 2 (",round(s$importance[2,2]*100, 2),"%)",sep=""), col=colors[condition], 
		pch=points[condition], cex=cmains, cex.lab=cmains, cex.main=cmains, cex.axis=caxis, lwd=lwdp)
	legend(legpos, legend=legtext, col=unique(colors), pch=points_ins[1:length(levels(phenotype))], cex=cleg, lty=ltyl, lwd=lwdl, box.lwd=0.9)

	if(save==TRUE){
		dev.off()
	}

}

specify_cancer <- function(usel){
	#Function that takes a keyword as an input and returns a new DGE object.
	#Arguments:
		#usel: Keyword specifying the type of comparison which will be analysed. 
		#	   The Keyword could be one of:
		#	       "all": All controls versus all cancers
		#		   "BSC": All controls versus all Skin, Breast and Cervical Cancers
		#		   "Breast": All controls versus breast cancers
		#		   "Skin": All controls versus skin cancers
		#		   "Cervix": All controls versus cervical cancers
		
	if(usel %in% c("Breast", "Skin", "Cervix")){
		selids <- rownames(y$samples[y$samples$group %in% c("Control", usel),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Control", "Control", "Case")
		}
  else if(usel=="all_no_skin"){
		selids <- rownames(y$samples[y$samples$group %in% c("Control","Brain", "Oesophagus", "Vulva", "Breast","Cervix"),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Control", "Control", "Case")
	  }
	else if(usel=="all"){
		newy <- y
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Control", "Control", "Case")
		} 
	else if(usel=="BSC"){
		selids <- rownames(y$samples[y$samples$group %in% c("Control", "Breast", "Skin", "Cervix"),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Control", "Control", "Case")
		}
	else if(usel=="BvC"){
		selids <- rownames(y$samples[y$samples$group %in% c("Breast", "Cervix"),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Breast", "Control", "Case")
		}
	else if(usel=="BvS"){
		selids <- rownames(y$samples[y$samples$group %in% c("Breast", "Skin"),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Breast", "Control", "Case")
		}
 	else if(usel=="CvS"){
		selids <- rownames(y$samples[y$samples$group %in% c("Cervix", "Skin"),])
		newy <- y[,selids, keep.lib.sizes=FALSE]
		newy$samples$group_old <- newy$samples$group
		newy$samples$group <- ifelse(newy$samples$group=="Cervix", "Control", "Case")
		}
	
	return(newy) #Returns the new DGE object
	}
 
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#### DIFFERENTIAL EXPRESSION ANALYSIS ####

do_test_deseq <- function(res, dds){
#Function that takes the result from DESeq2 and the initial Deseq2 object as its arguments
#and returns a detailed table for the given comparison with metrics for each reported gene.
#Arguments:
	#res: Result object of the deseq2 analysis. This object is the output of the 'results' 
	#	  function of the deseq2 package.
	#dds: Initial object used from the deseq2 package for the analysis. This object is the 
	#output of the 'DESeqDataSetFromMatrix' function of the deseq2 package.

	#Calculate how many cases are different (up or down) to controls based on their deseq2
	#final normalised values. If a gene is being reported as up-regulated, its Case samples
	#are reported different if they have at least two fold greater values than the median 
	#value of the Controls. If a gene is being reported as down-regulated, its Case samples
	#are reported different if they have at least two fold lower values than the median value
	#of the Controls.
	howc <- c()
	howcup <- c()
	howcdown <- c()
	for(i in rownames(deseq_norm)){
	howc[i] <- sum(abs(log2(deseq_norm[i, case_ids]/(median(deseq_norm[i,control_ids])+0.0001)+0.0001))>1) 
	howcup[i] <- sum(log2(deseq_norm[i, case_ids]/(median(deseq_norm[i,control_ids])+0.0001)+0.0001) > 1)
	howcdown[i] <- sum(log2(deseq_norm[i, case_ids]/(median(deseq_norm[i,control_ids])+0.0001)+0.0001)< -1)
	}
	res$cases_diff_controls <- howc
	res$cases_up_controls <- howcup
	res$cases_down_controls <- howcdown

	res$Dispersion <- dispersions(dds) #Report Dispersion
	res$logCPM <- log2(rowMeans(cpm(y))) #Report CPM
	res$logCPM2 <- log2(rowMeans(cpm(y)/finaly_filt$genes$Length)) #Report CPM'
	res$ctrl_mean_cpm <- cpmByGroup(y)[,1] #Report Controls CPM
	res$cas_mean_cpm <- cpmByGroup(y)[,2] #Report Controls CPM'
	res$ctrl_median_cpm <- rowMedians(y_cpms[,control_ids]) #Report Controls CPM median
	res$cas_median_cpm <- rowMedians(y_cpms[,case_ids]) #Report Cases CPM median
	res$ctrl_mean_raw <- rowMeans(y_ctrl$counts) #Report Controls CPM mean
	res$cas_mean_raw <- rowMeans(y_cases$counts) #Report Cases CPM mean
	res$ctrl_median_raw <- rowMedians(y_ctrl$counts) #Report Controls raw mean
	res$cas_median_raw <- rowMedians(y_cases$counts )#Report Cases raw mean
	res$ctrl_median_deseq <- rowMedians(deseq_norm[,control_ids]) #Report Controls deseq2 median normalised value
	res$cas_median_deseq <- rowMedians(deseq_norm[,case_ids]) #Report Cases deseq2 median normalised value
	res$ctrl_mean_deseq <- rowMeans(deseq_norm[,control_ids]) #Report Controls deseq2 mean normalised value
	res$cas_mean_deseq <- rowMeans(deseq_norm[,case_ids]) #Report Cases deseq2 mean normalised value
	res$LogFC_median <- log2(res$cas_median_cpm/res$ctrl_median_cpm) #Report Log2FC of mean CPM Controls vs. Cases
	res$color_fc <- ifelse(abs(res$log2FoldChange) > 0.5,"red", "black") #Color for plot. Not important
	res$color_cpm <- ifelse(res$logCPM2 > -5,"red", "black") #Color for plot. Not important

	return(res)
}
####DESeq2 plots

deseq_plots <- function(path, X){
#Function that takes a path and a prefix string as an input and creates
#Volcano plots, Pvalue to CPM plots and Dispersion plot for the wald DESeq2 outcome
#Arguments:
	#Path: Path that the plots will be saved in.
	#X: String that will be used as a part of the plot names.

	#Dispersion plot

	outind <- which(mcols(dds_wald)[,"dispOutlier"]=="TRUE")
	incind <- which(mcols(dds_wald)[,"dispOutlier"]=="FALSE")

	poch <- NULL
	poch[outind] = 1
	poch[incind] = 19

	cexp <- NULL
	cexp[outind] <- 0.8
	cexp[incind] <- 0.3


	png(filename=paste(path, "/FH01_dispersion_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispGeneEst"]), ylim=c(0,2), pch=19, cex=0.3, cex.lab=1.5, xlab="Log2 Ave CPM'", ylab="square-root Dispersion")
	points(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispersion"]), col="#1e90ff", pch=poch, cex=cexp)
	lines(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispFit"]), col="#FF0000")
	legend("bottomright", legend=c("gene-est", "fitted", "final", "outlier"), col=c("black", "#FF0000", "#1e90ff", "#1e90ff"), pch=c(19,19,19,1), cex=1.2)
	dev.off()
 
 	png(filename=paste(path, "/FH01_dispersion_small_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispGeneEst"]), ylim=c(0,2), pch=19, cex=0.5, cex.lab=3, cex.axis=2, xlab="Log2 Ave CPM'", ylab="square-root Dispersion")
	points(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispersion"]), col="#1e90ff", pch=poch, cex=cexp)
	#lines(res_dewald$logCPM2, sqrt(mcols(dds_wald)[,"dispFit"]), col="#FF0000")
	legend("bottomright", legend=c("gene-est", "final", "outlier"), col=c("black", "#1e90ff", "#1e90ff"), pch=c(19,19,1), cex=2)
	dev.off()

	#Colored volcano
	png(filename=paste(path, "/FH01_volcano_2_dewald_small_", X, "a.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$log2FoldChange, -log10(res_dewald$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, xaxt = "n", col=res_dewald$color_cpm)
	axis(1, at=seq((-ceiling(abs(min(res_dewald$log2FoldChange)))), ceiling(max(res_dewald$log2FoldChange)), by=0.5), cex.axis=2)
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=2)
	dev.off()
 
  png(filename=paste(path, "/FH01_volcano_2_dewald_", X, "a.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$log2FoldChange, -log10(res_dewald$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, xaxt = "n", col=res_dewald$color_cpm)
	axis(1, at=seq((-ceiling(abs(min(res_dewald$log2FoldChange)))), ceiling(max(res_dewald$log2FoldChange)), by=0.5))
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()
 

	#Volcano
	png(filename=paste(path, "/FH01_volcano_2_dewald_small_", X, "b.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$log2FoldChange, -log10(res_dewald$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, xaxt = "n")
	axis(1, at=seq((-ceiling(abs(min(res_dewald$log2FoldChange)))), ceiling(max(res_dewald$log2FoldChange)), by=0.5), cex.axis=2)
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	#legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()
 
 	#Volcano
	png(filename=paste(path, "/FH01_volcano_2_dewald_", X, "b.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$log2FoldChange, -log10(res_dewald$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, xaxt = "n")
	axis(1, at=seq((-ceiling(abs(min(res_dewald$log2FoldChange)))), ceiling(max(res_dewald$log2FoldChange)), by=0.5))
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	#legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()

	#Colored volcano
	png(filename=paste(path, "/FH01_volcano_2_delrt_", X, "a.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$log2FoldChange, -log10(res_delrt$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, xaxt = "n", col=res_delrt$color_cpm)
	axis(1, at=seq((-ceiling(abs(min(res_delrt$log2FoldChange)))), ceiling(max(res_delrt$log2FoldChange)), by=0.5))
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()
 
 #Colored volcano
	png(filename=paste(path, "/FH01_volcano_2_delrt_small_", X, "a.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$log2FoldChange, -log10(res_delrt$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, xaxt = "n", col=res_delrt$color_cpm)
	axis(1, at=seq((-ceiling(abs(min(res_delrt$log2FoldChange)))), ceiling(max(res_delrt$log2FoldChange)), by=0.5), cex.axis=2)
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()

	#Volcano
	png(filename=paste(path, "/FH01_volcano_2_delrt_", X, "b.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$log2FoldChange, -log10(res_delrt$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, xaxt = "n")
	axis(1, at=seq((-ceiling(abs(min(res_delrt$log2FoldChange)))), ceiling(max(res_delrt$log2FoldChange)), by=0.5), cex.axis=2)
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	#legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()
 
 #Volcano
	png(filename=paste(path, "/FH01_volcano_2_delrt_small_", X, "b.png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$log2FoldChange, -log10(res_delrt$pvalue), xlab = "log2 Fold-change", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, xaxt = "n")
	axis(1, at=seq((-ceiling(abs(min(res_delrt$log2FoldChange)))), ceiling(max(res_delrt$log2FoldChange)), by=0.5), cex.axis=2)
	abline(v=-0.5, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=0.5, col="black", lty=2)
	abline(v=1, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	#legend("bottomright", legend=c("logCPM2 > -5","logCPM2 < -5"), col = c("red", "black"), pch=19, cex=1.2)
	dev.off()

	#P-value to CPM' plot
	png(filename=paste(path, "/FH01_P_to_CPM2_dewald_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$logCPM2, -log10(res_dewald$pvalue), xlab = "log CPM'", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, col=res_dewald$color_fc)
	abline(v=-5, col="black", lty=2)
	abline(v=-4, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	axis(1, at=-4, cex.lab=1.5)
	legend("topright", legend=c("Absolute LogFC > 0.5","Absolute LogFC < 0.5"), col = c("red", "black"), pch=19)
	dev.off()
 
 #P-value to CPM' plot
	png(filename=paste(path, "/FH01_P_to_CPM2_dewald_small_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_dewald$logCPM2, -log10(res_dewald$pvalue), xlab = "log CPM'", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, col=res_dewald$color_fc)
	abline(v=-5, col="black", lty=2)
	abline(v=-4, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	axis(1, at=-4, cex.lab=1.5)
	legend("topright", legend=c("Absolute LogFC > 0.5","Absolute LogFC < 0.5"), col = c("red", "black"), pch=19, cex=2)
	dev.off()

	#P-value to CPM' plot
	png(filename=paste(path, "/FH01_P_to_CPM2_delrt_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$logCPM2, -log10(res_delrt$pvalue), xlab = "log CPM'", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=1.5, col=res_delrt$color_fc)
	abline(v=-5, col="black", lty=2)
	abline(v=-4, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	axis(1, at=-4, cex.lab=1.5)
	legend("topright", legend=c("Absolute LogFC > 0.5","Absolute LogFC < 0.5"), col = c("red", "black"), pch=19)
	dev.off()
 
 	#P-value to CPM' plot
	png(filename=paste(path, "/FH01_P_to_CPM2_delrt_small_", X, ".png", sep=""), units="in", width=11.5, height=8, res=600)
	par(mar=c(5.1,6,4.1,2.1))
	plot(res_delrt$logCPM2, -log10(res_delrt$pvalue), xlab = "log CPM'", ylab = "-log10 Pvalue ", pch=19, cex = 0.5, cex.lab=3, cex.axis=2, col=res_delrt$color_fc)
	abline(v=-5, col="black", lty=2)
	abline(v=-4, col="black", lty=2)
	abline(a=-log10(0.05), b=0, col="blue")
	axis(1, at=-4, cex.lab=2)
	legend("topright", legend=c("Absolute LogFC > 0.5","Absolute LogFC < 0.5"), col = c("red", "black"), pch=19, cex=2)
	dev.off()
}

###Further analysis of the results

decide_important <- function(table, pval_cutoff, padj_cutoff, abs_logFC_cutoff, mean_med_cutoff){
	#Function that takes a test outcome and pval cutoff, logFC and how many cases are different to 
	#mean/median of controls and returns an object with four items: 
	#first: list of gene names after the filtering (important genes)
	#second: list of 0s or 1s with the names of all the genes in the test. 
			 #0 if not important, 1 if important.
	#third: result of the test for the important genes
	#four: list of 0s or 1s with the names of all the genes in the analysis. 
			 #0 if not important, 1 if important.
		 
	#Setting default values	 
	if(missing(pval_cutoff)){
		pval_cutoff=1
	}
	if(missing(padj_cutoff)){
		padj_cutoff=1
	}
	if(missing(abs_logFC_cutoff)){
		abs_logFC_cutoff=0
	}
	if(missing(mean_med_cutoff)){
		mean_med_cutoff=0
	}
	
	#1st element
	imp_genes <- rownames(subset(table, pvalue < pval_cutoff
		& padj < padj_cutoff
		& log2FoldChange > 0
		& abs(log2FoldChange) > abs_logFC_cutoff 
		& cases_up_controls >= mean_med_cutoff))
	
	imp_genes <- c(imp_genes, rownames(subset(table, pvalue < pval_cutoff
		& padj < padj_cutoff
		& log2FoldChange < 0
		& abs(log2FoldChange) > abs_logFC_cutoff 
		& cases_down_controls >= mean_med_cutoff)))
	
	#2rd element:
	imp_genes_inte <- ifelse(rownames(table) %in% imp_genes, 1, 0)
	names(imp_genes_inte) <- rownames(table)
	
	#3rd element:
	imp_genes_table <- table[imp_genes,]
	
	#4th element:
	all_genes_inte <- ifelse(all_genes_detected %in% imp_genes, 1, 0)
	names(all_genes_inte) <- all_genes_detected
	
	return(list(imp_genes, imp_genes_inte, imp_genes_table, all_genes_inte))
}

merge_gene_sets <- function(a){
	#Function that takes list of important_genes objects as an input and returns 
	#a merged important_genes object by merging the input objects.
	#Example:
	#merge_gene_sets(a=list(imp_genes_fc5_p05_wald, imp_genes_fc5_p05_lrt, imp_genes_25perc_p05_wald, imp_genes_25perc_p05_lrt))

	if(length(a) == 2){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]]))
		}
	if(length(a) == 3){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]]))
		}
	if(length(a) == 4){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]]))
		}
	if(length(a) == 5){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]], a[[5]][[1]]))
		}
	if(length(a) == 6){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]], a[[5]][[1]], a[[6]][[1]]))
		}
	
	imp_genes_merge_inte <- ifelse(rownames(res_dewald) %in% imp_genes_merge_names, 1, 0)
	names(imp_genes_merge_inte) <- rownames(res_dewald)
	imp_genes_merge_table <- res_dewald[imp_genes_merge_names,]
	all_genes_merge_inte <- ifelse(all_genes_detected %in% imp_genes_merge_names, 1, 0)
	names(all_genes_merge_inte) <- all_genes_detected

	merge <- list(imp_genes_merge_names, imp_genes_merge_inte, imp_genes_merge_table, all_genes_merge_inte)
	return(merge)
}

bargene <- function(gene, path){
	#Function that takes a gene name and a path as its arguments and stores a barplot and a whiskers plot  
	#of the Deseq2 normalised values of this gene between Cases and Controls in the given path.
	#Arguments: 
		#gene: String. Should be a name of a genes included in the analysis.
		#path: Path that the plots will be saved in. 

	fname <- paste(path, "/", gene, "_barplot.png", sep="")
	png(fname, units="in", width=15.5, height=7, res=600)
	par(mar=c(6.1,5.1,4.1,2.1), mfrow=c(1,2))
	barplot(values_sorted[gene,], col=c(rep("blue", length(control_ids)), rep("red", length(case_ids))), las=2, ylab="CPM", main="gene", cex.names=0.5, xaxt="n")
	legend("topleft", col=c("blue","red"), pch=15, legend=c("Control", "Case"), cex=1.5)


	boxplot(values_sorted[gene,1:length(control_ids)], values_sorted[gene,c((length(control_ids)+1):ncol(finaly_filt$counts))], col=c("blue", "red"), 
	ylab="CPM", cex.lab=1.5, main=gene, outline=FALSE)
	axis(side=1, at=c(1,2), labels=c("Controls","Cases"), cex.axis=1.5, lty=0)
	#legend("bottom", legend=leg, bty="n")
	dev.off()
}

###Create a heatmap

make_heatmap <- function(genes, phenotype, path, name, save, width, height, dmethod){
	#Function that takes the output of the 'decide_important' or 'merge_gene_sets' function as its argument and 
	#returns a heatmap for the deseq2 normalised values of Cases vs Controls. The function produces a ggplot
	#heatmap.2. 
	#Arguments:
		#genes: Output object of 'decide_important' or 'merge_gene_sets' function.
		#phenotype: Factor. Phenotypes of the samples. The samples should be in the same
				   #order as in the 'deseq_norm' object.
		#save: Boolean. Whether to save or not the generated heatmap.
		#path: Path that the plots will be saved in. (if save=TRUE)
		#width: Set the width of the png output.
		#height: Set the height of the png output.
	
	
	heat_values <- deseq_norm[genes[[1]],]
	heat_genes <- rownames(genes[[3]][order(genes[[3]]$pvalue),])
	heat_values <- heat_values[heat_genes,]
	
	fin_heat_values <- t(scale(t(heat_values)))
	
	colors <- rep(colors_ins[1:length(levels(phenotype))], 2)
	points <- rep(points_ins[1:length(levels(phenotype))], 2)
	condition <- phenotype
		
	col.pan <- colorpanel(100, "blue", "white", "red")
	
	
	if(save==TRUE){
	png(paste(path, "/heatmap_", name, ".png" ,sep=""), units="in", width=width, height=height, res=600)
	}
	
	heatmap.2(fin_heat_values, distfun = function(x) dist(x, method = dmethod), col=col.pan, Rowv=TRUE, 
	scale="none", trace="none", dendrogram="both", cexRow=1.2, cexCol=1.2, density.info="none", margin=c(10,9), 
	lhei=c(1,5), lwid=c(0.4,3), colCol = colors[condition])
	#legend("legpos", legend=c(levels(phenotype)), col=unique(colors), lty=1, cex=1.2, lty=0, lwd=2, box.lwd=0.9)
	if(save==TRUE){
	dev.off()
	}
	
}

get_the_OR <- function(gores, sigs, all){
	#Function that takes the result of the go enrichment analysis as its argument
	#and returns the Fisher's Exact Odds Ratio for each reported GO term.
	#Arguments:
		#gores: Output of the 'GenTable' function of the 'topgo' package.
		#sigs: Vector of gene names which are DE and annotated to GO terms.
		#all: Vector of the all the gene names included in the analysis.

	s1 <- gores$Significant
	s2 <- sigs - s1
	b1 <- gores$Annotated
	b2 <- all - b1
	OR <- (s1*b2)/(s2*b1)

return(OR)
}

go_analysis <- function(allgenes, path, name){
	#Function that takes the output of the 'decide_important' or 'merge_gene_sets' function 
	#as its argument and returns a list of 9 objects, 3 for each specific ontology (BP=Biological Process,
	#CC=Cellular Component, MF=Molecular Function). Four different combinations of the topgo data object
	#were used for each ontology: 'elim' algorithm with KS test, 'elim' algorithm with exact fisher's test,
	#'weight' algorithm with exact fisher's test and 'weight01algorithm' with exact fisher's test. The result 
	#matrices are being saved as CSV files, one for each ontology.
	#First three elements: Three topgodata objects built for the analysis (see topgo vignette). BP, CC, MF respectively.
	#Second three elements: The result matrices of the 'GenTable' function (see topgo vignette.) BP, CC, MF respectively. 
						   #Each result matrix contains the results from all different combinations of algorithms and 
						   #statistical tests. These matrices also include columns reporting the names of significant genes of 
						   #each GO term (a column for the up-regulated and a column for the down-regulated genes) and a column
						   #reporting the Fisher's test Odds Ratio calculated for each term.
	#Third three elemets: Lists of DE genes assigned to each GO term. BP, CC, MF respectively. BP, CC, MF respectively.
	#Arguments:
		#allgenes: output of the 'decide_important' or 'merge_gene_sets' function.		
		#path: path to store the output CSV files. 

	dir.create(path=path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
		
	type <- deparse(substitute(type))
	nameit <- name

	#Annotate genes to Biological Process (BP) GO terms
	topgo_bp  <- new("topGOdata", ontology = "BP", allGenes = allgenes[[4]],
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Biological Process (BP) GO terms
	sigs_bp <- as.numeric(table(topgo_bp@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_bp <- length(topgo_bp@feasible[topgo_bp@feasible==TRUE])

	#Annotate genes to Cellular Component (CC) GO terms
	topgo_cc  <- new("topGOdata", ontology = "CC", allGenes = allgenes[[4]],
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Cellular Component (CC) GO terms
	sigs_cc <- as.numeric(table(topgo_cc@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_cc <- length(topgo_cc@feasible[topgo_bp@feasible==TRUE])

	#Annotate genes to Molecular Function (MF) GO terms
	topgo_mf  <- new("topGOdata", ontology = "MF", allGenes = allgenes[[4]],
	geneSel = function(k) k > 0, annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	#Genes annotated to Molecular Function (MF) GO terms
	sigs_mf <- as.numeric(table(topgo_mf@feasible[which(topgo_bp@allScores == 1)])["TRUE"])
	all_mf <- length(topgo_mf@feasible[topgo_bp@feasible==TRUE])

	#Gemes reported up- or down- regulated
	ups <- rownames(subset(allgenes[[3]], log2FoldChange > 0))
	downs <- rownames(subset(allgenes[[3]], log2FoldChange < 0))

	#Run the tests for BP, CC and MF	
	topgo_t1_bp <- runTest(topgo_bp, algorithm = "elim", statistic = "ks")
	topgo_t1_cc <- runTest(topgo_cc, algorithm = "elim", statistic = "ks")
	topgo_t1_mf <- runTest(topgo_mf, algorithm = "elim", statistic = "ks")
	topgo_t2_bp <- runTest(topgo_bp, algorithm = "elim", statistic = "fisher")
	topgo_t2_cc <- runTest(topgo_cc, algorithm = "elim", statistic = "fisher")
	topgo_t2_mf <- runTest(topgo_mf, algorithm = "elim", statistic = "fisher")
	topgo_t3_bp <- runTest(topgo_bp, algorithm = "weight", statistic = "fisher")
	topgo_t3_cc <- runTest(topgo_cc, algorithm = "weight", statistic = "fisher")
	topgo_t3_mf <- runTest(topgo_mf, algorithm = "weight", statistic = "fisher")
	topgo_t4_bp <- runTest(topgo_bp, algorithm = "weight01", statistic = "fisher")
	topgo_t4_cc <- runTest(topgo_cc, algorithm = "weight01", statistic = "fisher")
	topgo_t4_mf <- runTest(topgo_mf, algorithm = "weight01", statistic = "fisher")

	#Get the genes annotated per term
	all_GO_bp <- genesInTerm(topgo_bp)
	all_GO_cc <- genesInTerm(topgo_cc)
	all_GO_mf <- genesInTerm(topgo_mf)

	#Get the significant genes annotated per term
	SAM_bp <- lapply(all_GO_bp,function(x) x[x %in% names(allgenes[[4]][allgenes[[4]]==1])] )
	SAM_cc <- lapply(all_GO_cc,function(x) x[x %in% names(allgenes[[4]][allgenes[[4]]==1])] )
	SAM_mf <- lapply(all_GO_mf,function(x) x[x %in% names(allgenes[[4]][allgenes[[4]]==1])] )

	#Summarise the BP results
	topgo_res_bp <- GenTable(topgo_bp, elimKS = topgo_t1_bp, elimFisher = topgo_t2_bp, 
		weightFisher = topgo_t3_bp, weight01Fisher = topgo_t4_bp, orderBy = "weight01Fisher", 
		ranksOf = "weight01Fisher", topNodes = length(topgo_bp@graph@nodes))
	topgo_res_bp$OR <- get_the_OR(topgo_res_bp, sigs_bp, all_bp) #add the OR as well

	#Find which genes are up- or down- regulated for each reported BP term
	upgenes <- c()
	downgenes <- c()

	for(i in topgo_res_bp[,1]){

		upgenes[i]	<- paste(c(SAM_bp[[i]][SAM_bp[[i]] %in% ups]), collapse=", ")
		downgenes[i] <- paste(c(SAM_bp[[i]][SAM_bp[[i]] %in% downs]), collapse=", ")
	}

	topgo_res_bp$upgenes <- upgenes
	topgo_res_bp$downgenes <- downgenes

	#Summarise the CC results
	topgo_res_cc <- GenTable(topgo_cc, elimKS = topgo_t1_cc, elimFisher = topgo_t2_cc, 
		weightFisher = topgo_t3_cc, weight01Fisher = topgo_t4_cc, orderBy = "weight01Fisher", 
		ranksOf = "weight01Fisher", topNodes = length(topgo_cc@graph@nodes))
	topgo_res_cc$OR <- get_the_OR(topgo_res_cc, sigs_cc, all_cc) #add the OR as well

	#Find which genes are up- or down- regulated for each reported CC term
	upgenes <- c()
	downgenes <- c()

	for(i in topgo_res_cc[,1]){
		upgenes[i] <- paste(c(SAM_cc[[i]][SAM_cc[[i]] %in% ups]), collapse=", ")
		downgenes[i] <- paste(c(SAM_cc[[i]][SAM_cc[[i]] %in% downs]), collapse=", ")

	}

	topgo_res_cc$upgenes <- upgenes
	topgo_res_cc$downgenes <- downgenes

	#Summarise the MF results
	topgo_res_mf <- GenTable(topgo_mf, elimKS = topgo_t1_mf, elimFisher = topgo_t2_mf, 
		weightFisher = topgo_t3_mf, weight01Fisher = topgo_t4_mf, orderBy = "weight01Fisher", 
		ranksOf = "weight01Fisher", topNodes = length(topgo_mf@graph@nodes))
	topgo_res_mf$OR <- get_the_OR(topgo_res_mf, sigs_mf, all_mf) #add the OR as well

	#Find which genes are up- or down- regulated for each reported MF term
	upgenes <- c()
	downgenes <- c()
	z <- 1
	for(i in topgo_res_mf[,1]){
		upgenes[i] <- paste(c(SAM_mf[[i]][SAM_mf[[i]] %in% ups]), collapse=", ")
		downgenes[i] <- paste(c(SAM_mf[[i]][SAM_mf[[i]] %in% downs]), collapse=", ")
		z=z+1
	}

	topgo_res_mf$upgenes <- upgenes
	topgo_res_mf$downgenes <- downgenes

	#Write the output CSV files
	write.csv(data.frame(lapply(topgo_res_bp, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_bp.csv",sep=""))
	write.csv(data.frame(lapply(topgo_res_mf, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_mf.csv",sep=""))
	write.csv(data.frame(lapply(topgo_res_cc, as.character), stringsAsFactors=FALSE), file=paste(path,nameit,"_cc.csv",sep=""))
		
	return(list(topgo_bp, topgo_cc, topgo_mf, topgo_res_bp, topgo_res_cc, topgo_res_mf, SAM_bp, SAM_cc, SAM_mf))
	}

is_diff <- function(i){
	#To be used inside 'term_wise' function only!
	#Function that takes a gene as its argument and returns a list of case ids which upregulate this gene (if the gene 
	#is being reported from DESeq2 as upregulated) or downregulate this gene (if the gene is being reported from DESeq2
	#as downregulated). A case is considered to up/down regulate a gene if its deseq2 normalised value for that gene is
	#2-fold greater/less than the median deseq2 normalised value of the Controls.
		#Arguments:
			#i: gene name

	a <- names(which((res_dewald[i,]$log2FoldChange > 0 & log2(deseq_norm[i, case_ids]/median(deseq_norm[i,control_ids])) > 1) | 
	(res_dewald[i,]$log2FoldChange < 0 & log2(deseq_norm[i, case_ids]/median(deseq_norm[i,control_ids])) < -1)==TRUE ))
	
	return(a)
}

find_max <- function(red){
	#To be used inside 'term_wise' function only!
	#Function that takes several gene-sets as its input and calculates how many cases upregulate or downregulate each gene of each 
	#gene-set and returns the gene-sets that are up/down- regulated for the maximum number of cases and that number. A case is considered
	#to upregulate a gene if this gene was reported as upregulated from deseq2 and the deseq2 normalised value of that case for that gene 
	#is 2-fold greater than the median deseq2 normalised value of the Controls. A case is considered to downregulate a gene if this gene 
	#was reported as downregulated from deseq2 and the deseq2 normalised value of that case for that gene is 2-fold less than the median 
	#deseq2 normalised value of the Controls.
	#Arguments:
		#red: Matrix containing different combination of genes (like output of the 'combine' function). 
	redm <- lapply(seq_len(ncol(red)), function(i) red[,i])
	perlist <- c()
	for(i in redm){
	perlist <- c(perlist, length(Reduce(intersect, lapply(i, is_diff))))
	}
	pergenes <- paste(unique(unlist(redm[which(perlist == max(perlist))])), collapse = ", ")
	maxno <- max(perlist)
	return(list(unique(pergenes), maxno))
	}

term_wise <- function(go_res, save, path, name, ont){
	#Function that takes the output of the 'go_analysis' function as its argument and returns an enriched 
	#version of the 'GenTable' matrix (see 'topgo' vignette) with further sample-wise information for the
	#top reported GO terms (top GO terms) of the selected ontology; criteria: Weight01Fisher pvalue < 0.05,
	#OR>10, Significant genes >= 2). The output matrix contains the same columns as the go_analysis result 
	#matrices plus: 
		#A column reporting how many samples differentially express all of the sign. genes (100%) of each GO term
		#(see 'find_max' and 'is_diff' functions for sample-wise DE criteria). A column reporting which are these 
		#samples. Columns reporting how many samples differentially express other proportions of the sign. genes 
		#(one column for 75% of the genes & one column for 50% of the genes) of each top GO term (see 'find_max' 
		#and 'is_diff' functions for sample-wise DE criteria. And two columns reporting which are those genes 
		#(one column for 75% of the genes & one column for 50% of the genes). The function uses the 'is_diff' and
		#'find_max' functions to get this information. 
	#The function returns the resulted matrix.  
	#Arguments:
		#go_res: Output of the 'go_analysis' function
		#save: Boolean. Whether to save or not the output matrix as a CSV file.
		#path: Path to save the CSV file.
		#ont: Ontology. Can be "BP", "CC" or "MF". 
		

	if(ont %in% c("bp", "BP")){j=4}
	if(ont %in% c("cc", "CC")){j=5}
	if(ont %in% c("mf", "MF")){j=6}
		
	tops <- subset(go_res[[j]], as.numeric(Significant) >=2 & as.numeric(weight01Fisher) < 0.05 & as.numeric(OR) >= 10)
	
	sam_all <- c()
	no_all <- c()
	res75 <- c()
	res50 <- c()
	res25 <- c()
	z <- 1
	
	for(i in tops[[1]]){		
		
		genes <- go_res[[j+3]][[i]]
		
		#get the 100% genes number of samples
		glist <- lapply(genes, is_diff)
		dl <- Reduce(intersect, glist)
		sam_all[z] <- paste(dl,collapse=", ")
		no_all[z] <- length(dl)
		
		#get the reduced no. of genes. 
		#Divide the total number of genes into all the possible
		#gene-sets containing 75% of genes and apply the 'find_max' function.
		res75[z] <- list(unlist(list(find_max(red=combn(genes, 0.75*length(genes))))))
		
		#Same for 50%.
		res50[z] <- list(unlist(list(find_max(red=combn(genes, 0.5*length(genes))))))		
		z <- z+1
	}
	#Create the matrix
	tops$no_cases100_above_controls <- no_all
	tops$cases100_above_controls <- sam_all
	
	tops$no_cases75_above_controls <- unlist(lapply(res75, `[[`, 2))
	tops$gen_cases75_above_controls <- unlist(lapply(res75, `[[`, 1))
	
	tops$no_cases50_above_controls <- unlist(lapply(res50, `[[`, 2))
	tops$gen_cases50_above_controls <- unlist(lapply(res50, `[[`, 1))
	
	#Save the matrxi
	if(save==TRUE){
		write.csv(tops, file=paste(path, name, "_best_",ont, ".csv",sep=""))
	}
	
	return(tops)		
}

go_barplot <- function(topgo, path, name, test, mars, wi, he, cen, cela, cele, legpos, save=TRUE){
	#Function that takes the output of the 'go_analysis' function as its argument and returns a plot
	#visualising the number of each of the significant genes of the top GO terms for all three ontologies
	#(BP, CC, MF). top GO terms are the GO terms fullifilling these criteria: test pvalue < 0.05, 
	#OR>10, Significant genes >= 2).
	#Arguments:
		#go_res: Output of the 'go_analysis' function
		#save: Boolean. Whether to save or not the output matrix as a CSV file.
		#path: Path to save the CSV file.
		#test: For which test should the pvalue of a term be less than 0.05 in order
			  #for it to be included in the plot. Could be 'Weight01Fisher', 'Weight' or 'elim'.
		#mars: margins of the plot. (to be used in par(mar=c())).
		#wi: Weight parameter of the PNG 
		#he: Height parameter of the PNG
		#cen: Font size of yaxis names and inside-the-bars text. To be passed to cex graphical parameter
		#cela: Font size of xaxis label. To be passed to cex graphical parameter
		#cele: Font size of the legend. To be passed to cex graphical parameter
		#legpos:Keyword to be used to position the legend. Accepted keywords: 
		       #"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" 
		       #and "center" 	

	#Name and prefix to be used for naming the output PNG 
	
	#Top terms of the BP ontology analysis
	bp <- subset(topgo[[4]], as.numeric(Significant) >=2 & as.numeric(topgo[[4]][,test]) < 0.05 & as.numeric(OR) >= 10)
	bp <- bp[order(-as.numeric(bp[,test])),]

	#Top terms of the CC ontology analysis
	cc <- subset(topgo[[5]], as.numeric(Significant) >=2 & as.numeric(topgo[[5]][,test]) < 0.05 & as.numeric(OR) >= 10)
	cc <- cc[order(-as.numeric(cc[,test])),]

	#Top terms of the MF ontology analysis
	mf <- subset(topgo[[6]], as.numeric(Significant) >=2 & as.numeric(topgo[[6]][,test]) < 0.05 & as.numeric(OR) >= 10)
	mf <- mf[order(-as.numeric(mf[,test])),]

	bps <- nrow(bp)
	ccs <- nrow(cc)
	mfs <- nrow(mf)

	
	goterms <- Term(GOTERM)
	
	#Extract information from 'go_analysis' output object.
	#Create the empty lists for the loops
	gos_bp <- c()
	gos_val_bp <- c()
	gos_ann_bp <- c()
	gos_or_bp <- c()
	gos_p_bp <- c()

	gos_mf <- c()
	gos_val_mf <- c()
	gos_ann_mf <- c()
	gos_or_mf <- c()
	gos_p_mf <- c()

	gos_cc <- c()
	gos_val_cc <- c()
	gos_ann_cc <- c()
	gos_or_cc <- c()
	gos_p_cc <- c()
	
	#BP information
	for(i in 1:bps){
	gos_bp[i] <- goterms[[bp[i,1]]]
	gos_val_bp[i] <- as.numeric(bp[i,4])
	gos_ann_bp[i] <- as.numeric(bp[i,3])
	gos_or_bp[i] <- round(as.numeric(bp[i,10]), digits=0)
	gos_p_bp[i] <- format(as.numeric(bp[i,test]), scientific=TRUE)
	}

	#MF information
	for(i in 1:mfs){
	gos_mf[i] <- goterms[[mf[i,1]]]
	gos_val_mf[i] <- as.numeric(mf[i,4])
	gos_ann_mf[i] <- as.numeric(mf[i,3])
	gos_or_mf[i] <- round(as.numeric(mf[i,10]), digits=0)
	gos_p_mf[i] <- format(as.numeric(mf[i,test]), scientific=TRUE)
	}
	
	#CC information
  if(nrow(cc)>0){
  	for(i in 1:ccs){
     	gos_cc[i] <- goterms[[cc[i,1]]]
    	gos_val_cc[i] <- as.numeric(cc[i,4])
    	gos_ann_cc[i] <- as.numeric(cc[i,3])
    	gos_or_cc[i] <- round(as.numeric(cc[i,10]), digits=0)
    	gos_p_cc[i] <- format(as.numeric(cc[i,test]), scientific=TRUE)
    }
  } else {
    gos_cc <- c()
  	gos_val_cc <- c()
  	gos_ann_cc <- c()
  	gos_or_cc <- c()
  	gos_p_cc <- c()
	}
 

	#Combine all the information
	gos <- c(gos_cc, gos_mf, gos_bp)
	gos_val <- c(gos_val_cc, gos_val_mf, gos_val_bp)
	gos_ann <- c(gos_ann_cc, gos_ann_mf, gos_ann_bp)
	gos_or <- c(gos_or_cc, gos_or_mf, gos_or_bp)
	gos_p <- c(gos_p_cc, gos_p_mf, gos_p_bp)
  
  max_x <- max(gos_val) 
  
  if(missing(wi)){
    if(max_x %in% c(1:5)){
      wi <- max_x+24
    } else if(max_x %in% c(6:9)){
      wi <- max_x+11
    } else if(max_x %in% c(10:14)){
      wi <- max_x+12
    } else {
      wi <- (max_x*25)/11
    }
	}

	out_of <- paste(gos_val, gos_ann, sep="/")

  
	#Create the plot for less than 25 BP topgo terms

	if(bps < 25){
		if(save==TRUE){
			png(filename=paste(path, name, ".png", sep=""), units="in", width=wi, height=he, res=600)
		}

		par(mar=mars)
		y <- barplot(gos_val, horiz=TRUE, names.arg=gos, col=c(rep("#4885ed", ccs),rep("#3cba54", mfs), rep("orange", bps)),
		cex.names=cen, xaxt="n", las=2, xlab="Number of differentially expressed genes", cex.lab=cela )
		x <- 0.5*gos_val
		stats <- paste(out_of, gos_or, gos_p, sep=", ")
		text(x, y, stats, cex=cen)
		axis(1, at=seq(0, (max_x+1), by=1), labels=c(0:(max_x+1)), xlim=c(0,(max_x+1)), cex.axis=1.7)
		legend(legpos, pch=15, col=c("orange", "#3cba54", "#4885ed"), legend=c("Biological Process", "Molecular Function", "Cellular Component"), cex=cele, bty="n")
		abline(v=0)
		
		if(save==TRUE){
			dev.off()
		}
	}

	
	if(bps > 25){
    if(save==TRUE){
			png(filename=paste(path, name, ".png", sep=""), units="in", width=wi, height=he, res=600)
		}
		
    par(mar=mars)
		y <- barplot(log2(gos_val), horiz=TRUE, names.arg=gos, col=c(rep("#4885ed", ccs),rep("#3cba54", mfs), rep("orange", bps)),
		cex.names=cen, xaxt="n", las=2, xlab="Number of differentially expressed genes", cex.lab=cela )
		x <- 0.5*log2(gos_val)
		stats <- paste(out_of, gos_or, gos_p, sep=", ")
		text(x, y, stats, cex=cen)
		axis(1)
		legend(legpos, pch=15, col=c("orange", "#3cba54", "#4885ed"), legend=c("Biological Process", "Molecular Function", "Cellular Component"), cex=cele, bty="n")
		abline(v=0)
		
		if(save==TRUE){
			dev.off()
		}
	}

	} #End of function

#Pubmed Search

	pubmed_search <- function(genes, coterm, organism, path, name){
	#Function that takes gene or genes as its argument, a term (String)and the organism
	#that these genes belong to. The function will perform ncbi searches in the pubmed, 
	#pmc and gene databases and it will return a table reporting the condition of each gene
	#according to the data (up/down regulated), its gene ID, description, summary and the 
	#number of pubmed and PMC articles comentioning each gene with the co-term. 
	#Arguments:
		#genes: output of 'decide_important' or 'merge_gene_sets' function.
		#coterm: word or phrase (string) that will be searched against pubmed and pmc databases
			    #for co-mentioning with the name of each gene.
		#organism: Organism that the genes belong to. (e.g. Homo Sapiens)
		#path: The path that the output will be saved to.
		#name: Name of the output file in the drive.
	
	#Generating empty lists for the gene loop	
	genez <- c()
	pus <- c()
	pms <- c()
	idz <- c()
	desc <- c()
	summ <- c()
	updn <- c()
	gene_spec <- c()
	z <- 1
	
	#loop over each gene
	for(i in genes[[1]]){
		
		s1 <- entrez_search("pubmed", paste(i,' AND ', '"', coterm, '"', sep=" ")) #Pubmed Search
		s2 <- entrez_search("pmc", paste(i,' AND ', '"', coterm, '"', sep=" ")) #PMC Search
		s3 <- entrez_search("gene", paste(i,'[Gene/Protein Name] AND ', organism, '[Organism]', sep="")) #Gene Search
		
		#Creating the vectors with the information
		genez[z] <- i
		pus[z] <- as.numeric(s1$count)
		pms[z] <- as.numeric(s2$count)

		if(s3$count > 0){ #If the gene search returned results:
			s4 <- entrez_summary(db="gene", id=s3$ids[1])
			idz[z] <- s4$uid
			desc[z] <- s4$description
			summ[z] <- s4$summary
			}

		if(s3$count == 0){ #If the gene search returned no hits:
			idz[z]= "-"
			desc[z]= "-"
			summ[z]= "-"
			}
		gene_spec[z] = ifelse(genes[[3]][i,]$log2FoldChange < 0, "Down", "Up")
		z <- z+1
	}
	
	#Final dataframe
	f_df <- data.frame(gene_spec, pus, pms, idz, desc, summ)
	rownames(f_df) <- genez
	colnames(f_df) <- c("Condition", "Pubmed_articles", "Pmc_articles", "ID", "Description", "Summary")
	
	write.csv(f_df, paste(path, "/", "pms_",name, ".csv", sep=""))

	return(f_df) #Return the dataframe

	}

#Example: pubmed_search(imp_genes_merge_wald, "Cancer", "Homo Sapiens", "~", "dio")
	
pubmed_barplot <- function(pubmed, save, path, name, wid, types="pubmed"){
	#Function that takes the output of the 'pubmed_search'function as its input 
	#and plots the numbers of pubmed articles. The function log2-transforms the 
	#column of pubmed articles of the pubmed_search object and plots them using 
	#ggplot2.
	#Arguments:
		#pubmed: output of the 'pubmed_search'function.
		#path: where the output image will be stored.
		#save: Boolean: Whether to store the image or not.
		#name: Name of the output image in the drive.
		#wid: width of the image i.e. 11.5.
	
	scf <- (nrow(pubmed)/25)*9
	
	
	
	#Log2 transform and sort the data
	df <- pubmed
  
  if(types=="pubmed"){
  	df$lpubmed <- log2(df[,2]+1)
    df <- df[order(df$lpubmed),]
    plabel <- df[,2]
    pylab <- "Log2 Number of Pubmed Articles"
  } else if(types=="pmc"){
    df$lpubmed <- log2(df[,3]+1)
    df <- df[order(df$lpubmed),]
    plabel <- df[,3]
    pylab <- "Log2 Number of PMC Articles"
  }
	
	df$Gene <- factor(rownames(df), levels=rownames(df))
	
	#Generate the plot
	theme_set(theme_bw())
	p <- ggplot(df, aes(x=df$Gene, y=df$lpubmed, label=plabel)) + geom_point(aes(color=df$lpubmed), size=11) + 
	scale_colour_gradient(low = "#80d0c7", high="#13547a") +
	geom_text(color="white", size=3.5, fontface = "bold") + xlab("Genes") + ylab(pylab) + 
	theme(axis.title.x = element_text(size=14), axis.text.x = element_blank(), 
	axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position="none") + coord_flip()
	
	#Save the image
	if(save==TRUE){
		ggsave(filename=paste("pms_ ",types, "_", name, ".png", sep=""), plot = p, device = NULL, path = path,
		scale = 1, width = wid, height = scf, dpi=800)
	}

	
}	

###R script that builds and tests a lasso regression model to predict binomial or multinomial outcome###

set.seed(12345)

#Loading of all the necessary packages
library(limma)
library(Rsubread)
library(edgeR)
library(DESeq2)
library("rentrez")
library(MASS)
library(glmnet)
library(fitdistrplus)
library(DescTools)
library(ROCR)
library(foreach)
library(doMC)
library(Publish)

load_variables <- function(path){
#Function that loads the stored variables from the differential 
#expression analysis given a specific path.

	finaly_filt <<- readRDS(paste(path, "/finaly_filt.rds", sep="")) #To get the phenotype data
	deseq_norm <<- readRDS(paste(path, "/deseq_norm.rds", sep="")) #Matrix of the final normalised values
	res_delrt <<- readRDS(paste(path, "/res_delrt.rds", sep="")) #LRT test outcome
	res_dewald <<- readRDS(paste(path, "/res_dewald.rds", sep="")) #Wald test outcome
	all_genes_detected <<- readRDS(paste(path, "/all_genes_detected.rds", sep="")) #all genes detected
}


decide_important <- function(table, pval_cutoff, padj_cutoff, abs_logFC_cutoff, mean_med_cutoff){
#Function that takes a test outcome and pval cutoff, logFC and how many cases are different to 
#mean/median of controls and returns an object with four items: 
#first: list of gene names after the filtering (important genes)
#second: list of 0s or 1s with the names of all the genes in the test. 
		 #0 if not important, 1 if important.
#third: result of the test for the important genes
#four: list of 0s or 1s with the names of all the genes in the analysis. 
		 #0 if not important, 1 if important.
		 
	#Setting default values	 
	if(missing(pval_cutoff)){
		pval_cutoff=1
	}
	if(missing(padj_cutoff)){
		padj_cutoff=1
	}
	if(missing(abs_logFC_cutoff)){
		abs_logFC_cutoff=0
	}
	if(missing(mean_med_cutoff)){
		mean_med_cutoff=0
	}
	
	#1st element
	imp_genes <- rownames(subset(table, pvalue < pval_cutoff
		& padj < padj_cutoff
		& log2FoldChange > 0
		& abs(log2FoldChange) > abs_logFC_cutoff 
		& cases_up_controls > mean_med_cutoff))
	
	imp_genes <- c(imp_genes, rownames(subset(table, pvalue < pval_cutoff
		& padj < padj_cutoff
		& log2FoldChange < 0
		& abs(log2FoldChange) > abs_logFC_cutoff 
		& cases_down_controls > mean_med_cutoff)))
	
	#2rd element:
	imp_genes_inte <- ifelse(rownames(table) %in% imp_genes, 1, 0)
	names(imp_genes_inte) <- rownames(table)
	
	#3rd element:
	imp_genes_table <- table[imp_genes,]
	
	#4th element:
	all_genes_inte <- ifelse(all_genes_detected %in% imp_genes, 1, 0)
	names(all_genes_inte) <- all_genes_detected
	
	return(list(imp_genes, imp_genes_inte, imp_genes_table, all_genes_inte))
}


merge_gene_sets <- function(a){
	#Function that takes list of important_genes objects as an input and returns 
	#a merged important_genes object by merging the input objects.
	#Example:
	#merge_gene_sets(a=list(imp_genes_fc5_p05_wald, imp_genes_fc5_p05_lrt, imp_genes_25perc_p05_wald, imp_genes_25perc_p05_lrt))

	if(length(a) == 2){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]]))
		}
	if(length(a) == 3){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]]))
		}
	if(length(a) == 4){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]]))
		}
	if(length(a) == 5){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]], 
		                                  a[[5]][[1]]))
		}
	if(length(a) == 6){
		imp_genes_merge_names <- unique(c(a[[1]][[1]], a[[2]][[1]], a[[3]][[1]], a[[4]][[1]], 
		                                  a[[5]][[1]], a[[6]][[1]]))
		}
	
	imp_genes_merge_inte <- ifelse(rownames(res_dewald) %in% imp_genes_merge_names, 1, 0)
	names(imp_genes_merge_inte) <- rownames(res_dewald)
	imp_genes_merge_table <- res_dewald[imp_genes_merge_names,]
	all_genes_merge_inte <- ifelse(all_genes_detected %in% imp_genes_merge_names, 1, 0)
	names(all_genes_merge_inte) <- all_genes_detected

	merge <- list(imp_genes_merge_names, imp_genes_merge_inte, imp_genes_merge_table, all_genes_merge_inte)
	return(merge)
}


check_dist <- function(imp_genes, ids, plot, save, height, width, path, prefix){
	#Function that takes an important_genes object, desired sample ids and some parameters as an 
	#input and returns a graph showing the fit of lognormal, gamma, weibull and normal distribution
	#to the values of each gene and the result of KS test comparing each distr. with the distr. 
	#of each gene. Then returns the best lognormal parameters for each gene that can be used
	#for simulating samples.
	#Arguments:
		#plot: Boolean: do you want to create the fit plot?
		#save: Boolean: do you want to save the plot and the 
						  #table with the metrics of the fittings?
		#height/width: png parameters
		#path: path to save image and table
		#prefix: add a prefix to output files

	z <- 1 #z parameter for loop

	#Define some lists for the loop
	gamma_para <- c()
	wei_para <- c()
	logn_para <- c()
	nor_para <- c()
	genes_all <- c()
	all_ks_p <- c()
	
	#Prepare names for saving files
	image1 <- paste(substr(x=deparse(substitute(imp_genes)), 11, 
	          nchar(deparse(substitute(imp_genes)))-5), ".png", sep="")
	prefix <- paste(prefix, "_", sep="")
	image2 <- paste("Distributions_", prefix, image1, sep="")
	image3 <- paste(path, image2, sep="/")

	#Save the fitting image
	if(save==TRUE){
		
		png(image3, units="in", width=width, height=height, res=600)
		par(mfrow=c(ceiling(length(imp_genes)/5),5))
	}
		
	for(i in c(1:length(imp_genes))){

		x <- deseq_norm[imp_genes, ids][i,]
		#Printing for troubleshooting
		#print("------------------------------------------")
		#print(i)
		#print(ids)
		
		#getting the density of the gene values for the given samples
		h <- hist(x, breaks=30, plot = FALSE)

		#fitting the gamma distribution
		gam <- fitdist(x, distr = "gamma", method = "mge", lower=c(0,0))
		#print("gam")
		gamma_para[z] <- gam$estimate[1]
		gamma_para[z+1] <- gam$estimate[2]

		#fitting the weibull distribution
		wei <- fitdist(x+0.00000000000000000000000000001, distr = "weibull", lower=c(0,0))
		#print("wei")
		wei_para[z] <- wei$estimate[1]
		wei_para[z+1] <- wei$estimate[2]

		#fitting the lognormal distribution
		nlo <- fitdist(x+0.1, distr = "lnorm", method = "mge")
		#print("nlo")
		logn_para[z] <- nlo$estimate[1]
		logn_para[z+1] <- nlo$estimate[2]
	
		#fitting the normal distribution
		nor_para[z] <- mean(x)
		nor_para[z+1] <- sd(x)

		genes_all[i] <- rownames(deseq_norm[imp_genes, ids])[i]


	if(plot==TRUE){	
	
		#Creating the plot: Each gene's distribution and the four fits on top
		
		#Getting the KS tests for each gene vs each fit
		ksn <- ks.test(x, "pnorm", mean(x), sd(x))
		ksg <- ks.test(x, "pgamma", as.numeric(gam$estimate[1]), as.numeric(gam$estimate[2]))
		kwe <- ks.test(x, "pweibull", as.numeric(wei$estimate[1]), as.numeric(wei$estimate[2]))
		knl <- ks.test(x, "plnorm", as.numeric(nlo$estimate[1]), as.numeric(nlo$estimate[2]))
		
		#Plot the original distribution of the gene's values
		xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
		yhist <- do.call(dnorm, c(list(x = xhist), as.list(c(mean(x), sd(x)))))
		ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density))
		
		#Plot the normal distribution fit of the original gene-wise values
		hist(x, freq = FALSE, xlab = "Data", ylim = c(0, ymax),
		main = rownames(deseq_norm[imp_genes, ids])[i],
		breaks = h$breaks)
		lines(xhist, yhist, lty = 1, col = "red")
		norm_AUC <- AUC(xhist, yhist, method="trapezoid", na.rm=FALSE)
		#print(paste("Norm AUC:",norm_AUC, sep=" "))

		#Plot the gamma distribution fit of the original gene-wise values
		yhist <- do.call(dgamma, c(list(x = xhist), as.list(c(gam$estimate[1], gam$estimate[2]))))
		ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density))
		lines(xhist, yhist, lty = 1, col = "blue")
		gamma_AUC <- AUC(xhist, yhist, method="trapezoid", na.rm=FALSE)
		#print(paste("Gamma AUC:",gamma_AUC, sep=" "))
		
		#Plot the weibull distribution fit of the original gene-wise values
		yhist <- do.call(dweibull, c(list(x = xhist), as.list(c(wei$estimate[1], wei$estimate[2]))))
		ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density))
		lines(xhist, yhist, lty = 1, col = "green")
		weibull_AUC <- AUC(xhist, yhist, method="trapezoid", na.rm=FALSE)
		#print(paste("Weibull AUC:",weibull_AUC, sep=" "))

		#Plot the lognormal distribution fit of the original gene-wise values
		yhist <- do.call(dlnorm, c(list(x = xhist), as.list(c(nlo$estimate[1], nlo$estimate[2]))))
		ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density))
		lines(xhist, yhist, lty = 1, col = "purple")
		lnorm_AUC <- AUC(xhist, yhist, method="trapezoid", na.rm=FALSE)
		#print(paste("Lnorm AUC:",lnorm_AUC, sep=" "))

		#Format the KS reported pvalues for adding them to the plot and table
		norm_p <- paste("normal p", formatC(ksn$p.value, format = "f", digits = 3), sep="=")
		gamma_p <- paste("gamma p", formatC(ksg$p.value, format = "f", digits = 3), sep="=")
		wei_p <- paste("weibull p", formatC(kwe$p.value, format = "f", digits = 3), sep="=")
		nlo_p <- paste("lognorm p", formatC(knl$p.value, format = "f", digits = 3), sep="=")
		legend("topright", inset=c(0.08,0), legend=c(norm_p, gamma_p, wei_p, nlo_p), 
		       col=c("red", "blue", "green", "purple"), pch=0, cex=1)

		all_ks_p <- c(all_ks_p, ksn$p.value, norm_AUC, ksg$p.value, gamma_AUC, kwe$p.value, weibull_AUC, knl$p.value, lnorm_AUC)
	}
	
	z <- z+2	
	
	}

	if(plot==TRUE){
		
		all_ks_p_matrix <- matrix(all_ks_p, byrow=TRUE, ncol=8)
		colnames(all_ks_p_matrix) <- c("normal ks p", "normal AUC", "gamma ks p", "gamma AUC", 
		        "weibull ks p", "weibull AUC", "log_normal ks p", "log_normal AUC")
		
		rep1 <- paste(substr(x=deparse(substitute(imp_genes)), 11, nchar(deparse(substitute(imp_genes)))-5), ".csv", sep="")
		prefix <- paste(prefix, "_", sep="")
		rep2 <- paste("Distributions_comparison_", prefix, rep1, sep="")
		rep3 <- paste(path, rep2, sep="/")
		
		
		write.csv(all_ks_p_matrix, file=rep3) #Save the sample-wise fitting metrics

	}


	if(save==TRUE){
			dev.off() #Save the gene-wise fitting plot
	}


	gamma_param <- matrix(c(gamma_para), byrow=TRUE, ncol=2, nrow=length(imp_genes))
	rownames(gamma_param) <- genes_all

	wei_param <- matrix(c(wei_para), byrow=TRUE, ncol=2, nrow=length(imp_genes))
	rownames(wei_param) <- genes_all

	logn_param <- matrix(c(logn_para), byrow=TRUE, ncol=2, nrow=length(imp_genes))
	rownames(logn_param) <- genes_all

	nor_param <- matrix(c(nor_para), byrow=TRUE, ncol=2, nrow=length(imp_genes))
	rownames(nor_param) <- genes_all

	return(list(1, 2, logn_param, nor_param)) #Return the parameters of the several fittings
											  #on the original gene-wise values
}

simulate_cases <- function(dists, howmany, prefix){
#Function that takes the output of the check_dist function (dists), 
#a real number (howmany) and a string (prefix) as its arguments and returns a simulated 
#gene-wise gene expression value matrix of N=howmany samples. The 
#simulations are based on the lognormal distribution parameters of the dists object.  
#Each simulated sample has the prefix with an increased number as name.

z <- 1

genes_for_sim <- c()
simulated_cases <- c()

for(i in c(1 : nrow(dists))){
	genes_for_sim[i] <- rownames(dists)[i] 
	sim <- rlnorm(howmany, dists[i,][1], dists[i,][2]) #Simulation based on the lognormal parameters
	                                                   #of the dists object	
	simulated_cases <- cbind(simulated_cases, sim)
	}

colnames(simulated_cases) <- genes_for_sim
rownames(simulated_cases) <- paste(prefix, "_sim_sample_", c(1:howmany), sep="") 
return(simulated_cases) #return the matrix
}

conf_matrix <- function(all_test, pred){
#Function that takes two factors of binary values (0,1) as an input:
#Arguments:
	#all_test: Binary factor. Represents the actual phenotype of the 
	#		   samples of a test set.
	#          0 for Controls, 1 for Non-Controls
	#pred:     Binary factor. Represents the predicted phenotype of the 
	#		   samples of a test set. 
	#          0 for Controls, 1 for Non-Controls
#The two factors should have the same order of samples. The function
#compares between the two factors and returns a confusion matrix of type:
#
#    				| True Control | True Case
#___________________|______________|__________
#Predicted Control  |      TN      |     FN
#___________________|______________|__________
#Predicted Case     |      FP      |     TP

	a <- 0
	b <- 0
	c <- 0
	d <- 0

	for ( i in 1:length(pred) ) {
	  
	  if ( all_test[i] == 0 & pred[i] == 0 ) {
	   a <- a + 1
	  }

	  if ( all_test[i] == 0 & pred[i] == 1 ) {
	   b <- b + 1
	  }
	  
	  if ( all_test[i] == 1 & pred[i] == 0 ) {
	   c <- c + 1
	  }
	  
	  if ( all_test[i] == 1 & pred[i] == 1 ) {
	   d <- d + 1
	  }
	  
	}	

	fm <- matrix(c(a,b,c,d), byrow=TRUE, ncol=2, nrow=2) #Create the matrix
	rownames(fm) <- c("0","1")
	colnames(fm) <- c("Control", "Case")
	return(fm)
	}

get_metrics <- function(fn, fp, tn, tp, nobs){
#Function that takes 5 real numbers an an input.
	#Arguments:
	#FN: Number of False Negatives of a prediction
	#FP: Number of False Positives of a prediction
	#TN: Number of True Negatives of a prediction
	#TP: Number of True Positives of a prediction
	#nobs: Sample-size of the test set
#The function returns a named vector with in depth
#metrics of the prediction (PPV, Sensitivity ...)
	
	sens <- tp/(tp+fn)
	spec <- tn/(tn+fp)
	
	if((tp+fp)==0){
		ppv=0
	} else {
		ppv <- tp/(tp+fp)
	}
	
	err <- (fp+fn)/(tn+tp+fn+fp)
	ci_low <- ifelse((err - (1.96*sqrt((err*(1-err))/nobs)))<0, 0, (err - (1.96*sqrt((err*(1-err))/nobs))))
	ci_high <- err + (1.96*sqrt((err*(1-err))/nobs))
	metrics <- c(tn, fp, fn, tp, sens, spec, ppv, err, ci_low, ci_high, nobs)
	names(metrics) <- c("TN", "FP", "FN", "TP", "Sensitivity", "Specificity", "PPV", "Error", "CI_Low", 
	                    "CI_High", "Test-set N")
	
	return(metrics)
}

get_metrics_conf <- function(confmx, nobs){
	#Function that takes a confusion matrix (output of the 
	#conf_matrix function) as an input and the sample-size of
	#the test set and returns a named vector with in depth
	#metrics of the prediction (PPV, Sensitivity ...)
		
	tp <- confmx[2,2]
	tn <- confmx[1,1]
	fp <- confmx[1,2]
	fn <- confmx[2,1]
	
	sens <- tp/(tp+fn)
	spec <- tn/(tn+fp)
	ppv <- tp/(tp+fp)
	err <- (fp+fn)/(tn+tp+fn+fp)
	ci_low <- ifelse((err - (1.96*sqrt((err*(1-err))/nobs)))<0, 0, (err - (1.96*sqrt((err*(1-err))/nobs))))
	ci_high <- err + (1.96*sqrt((err*(1-err))/nobs))
	metrics <- c(tn, fp, fn, tp, sens, spec, ppv, err, ci_low, ci_high, nobs)
	names(metrics) <- c("TN", "FP", "FN", "TP", "Sensitivity", "Specificity", "PPV", "Error", 
	                    "CI_Low", "CI_High", "Test-set N")
	
	return(metrics)
}

###Predict the result

lasso_prediction <- function(genes, split_factor, ci, ki, lambda_method, save=TRUE, prefix, id, path){
	#Function that splits the given gene expression data to a training
	#and a test set, it fits a lasso regression model to the training set
	#(binomial distribution) by using the 'glmnet' package, it cross-validates
    #it and uses this cross-validated model to predict the bionomial phenotype
	#of the test set. The function returns the detailed metrics of the prediction.
    #The Function can also upscale the training set by using the 'check_dist' and 
	#'simulate_cases' functions. 
	
	#Arguments: 
		#genes: 'important_genes' object (output of the decide_important function).
		#split_factor: Number from 0 to 1. Determines the splitting ratio of the data to
		#              training and test set. (e.g. split_factor=0.7: 70% of the data is
		#              being held for training, 30% for test).
		#ci: Number or the string "actual". Upscaling parameter: Determines how many simulated 
		#    Controls should be added in the existing cohort of original controls of the training set.
		#    If ci is number, it multiplies the number of the training set control samples (e.g. when 2 
		#    the number of Controls of the training set is duplicated). If "actual" the 
		#    controls are not being upscaled.
		#ki: Number or the string "actual". Upscaling parameter: Determines how many simulated 
		#    Non-Controls should be added in the existing cohort of original Non-Controls of the training set.
		#    If ci is number, it represents a fraction of the number of the control samples.
		#	 (e.g. when 2 the number of Non-Controls of the training set should be double of the corresponding
		#    number of Control samples in the training set. If "actual" the controls are not being upscaled.
		#lamda_method: Number or "min" or "1se". Lambda value that will be used for glmnet prediction. 
		#              See pred function in glmnet package documentation for more information.
		#save: Boolean. Save the rds object of the fitted model or not.
		#prefix: String. Prefix of the name of the saved rds object.
		#id: Number. Combined with prefix create the name of the saved rds object.
		#            name: prefix_id.rds
		#path: path to store the output rds file.
	
	xcontrols <- t(deseq_norm[genes, control_ids]) #Data of Control samples
	xcases <- t(deseq_norm[genes, case_ids]) #Data of Non-Control samples

	ncontrols <- nrow(xcontrols)
	ncases <- nrow(xcases)

	controls_sam_id <- sample(1:ncontrols, round(0.75*ncontrols)) #Ids held for training
	cases_sam_id <- sample(1:ncases, round(0.75*ncases))

	
	if(ci=="actual" | as.numeric(ci)==1){ #No upscaling of the Controls
		overall_controls_length <- length(controls_sam_id)
		train_controls <- xcontrols[controls_sam_id,]
		
	} else { #Upscaling of the Controls
		overall_controls_length <- ceiling(as.numeric(ci)*length(controls_sam_id))
		new_controls_length <- (overall_controls_length - length(controls_sam_id))
		controls_sam_id_for_sim <- rownames(xcontrols[controls_sam_id,])
		dists_control <- check_dist(genes, controls_sam_id_for_sim, prefix="control", 
		                 path="../results/24_05_2018/predict_cases", save=FALSE, plot=FALSE)
		simulated_controls <- simulate_cases(dists_control[[3]], new_controls_length, "control")
		train_controls <- rbind(xcontrols[controls_sam_id,], simulated_controls)
	}
	
	
	if(ki=="actual" | as.numeric(ki)==1){ #No upscaling of the Non-Controls
		overall_cases_length <- length(cases_sam_id)
		train_cases <- xcases[cases_sam_id,]
	} else { #Upscaling of the Non-Controls
		overall_cases_length <- ceiling(as.numeric(ki)*overall_controls_length)
		new_cases_length <- overall_cases_length - length(cases_sam_id)
		cases_sam_id_for_sim <- rownames(xcases[cases_sam_id,])
		dists_case <- check_dist(genes, cases_sam_id_for_sim, prefix="case", 
		              path="../results/24_05_2018/predict_cases", save=FALSE, plot=FALSE)
		simulated_cases <- simulate_cases(dists_case[[3]], new_cases_length, "case")
		train_cases <- rbind(xcases[cases_sam_id,], simulated_cases)
	}
		
		#Combine original and simulated samples of the training set
		x.train <- rbind(train_controls, train_cases) #Data
		y.train <- ifelse(finaly_filt$samples[rownames(x.train),]$group == "Control" | 
		           substr(rownames(x.train), 1, 7)=="control", 0, 1) #Phenotype
		y.train[is.na(y.train)] <- 1
    y.train <- factor(y.train, levels=c(0,1))
    print("-------------------------------------------")
    print("sf,ci,ki")
    print(c(split_factor, ci, ki))
    print("train_cases")		
    print(rownames(xcases[cases_sam_id,]))
    print("train_controls")
    print(rownames(xcontrols[controls_sam_id,]))
    #print(rownames(x.train))
    print("x.train")
    print(rownames(x.train))
    print("y.train")
    print(y.train)

    
		
		#Create the test set data
		x.test <- rbind(xcontrols[-controls_sam_id,], xcases[-cases_sam_id,]) 
    rownames(x.test) <- c(rownames(xcontrols)[-controls_sam_id], rownames(xcases)[-cases_sam_id])
    print(rownames(x.test))
    print("-------------------------------------------")
		#Fit the lasso regression
		fit.lasso <- glmnet(x.train, y.train, family="binomial", alpha=1)
		all_test <- ifelse(finaly_filt$samples[rownames(x.test),]$group=="Control", 0 ,1)
		cvfit = cv.glmnet(x.train, y.train, family = "binomial", 
		        type.measure = "class", alpha=1)
	
	print("OK Fitting")
	
	#Get the lambda value of the model
	if(lambda_method == "min"){
		lambda = cvfit$lambda.min
		lambda_m = "lambda.min"
		cofs = coef(cvfit, s="lambda.min")
	}
	if(lambda_method == "1se"){
		lambda = cvfit$lambda.1se
		lambda_m = "lambda.1se"
		cofs = coef(cvfit, s="lambda.1se")
	} 
	
	if(is.numeric(lambda_method)==TRUE){
		lambda=lambda_method
		lambda_m = lambda_method
		cofs = coef(cvfit, s=lambda_method)
		}
	
	#Predict the phenotype of the test-set
	pred <- predict(cvfit, newx = x.test, s = lambda_m, type = "class")
	pred1 <- predict(cvfit, newx = x.test, s = lambda_m, type = "link")
	pred2 <- predict(cvfit, newx = x.test, s = lambda_m, type = "coefficients")
	pred3 <- predict(cvfit, newx = x.test, s = lambda_m, type = "response")
	pred4 <- predict(cvfit, newx = x.test, s = lambda_m, type = "nonzero")
	row_sub <- apply(pred2, 1, function(row) all(row !=0 ))
	gett <- pred2[row_sub,]
	genesf <- names(gett)[2:length(gett)]

		
	print("OK Pred")
	ngenes <- length(cofs[,1][cofs[,1]!=0])-1
  
	
	fm <- conf_matrix(all_test, pred) #Get the confusion matrix of the prediction
	TN <- fm[1,1]
	FP <- fm[1,2]
	FN <- fm[2,1]
	TP <- fm[2,2]
	
	#Create the final matrix of the prediction metrics
	final_list <- c(id, get_metrics(tn=TN, fp=FP, fn=FN, tp=TP, nobs=nrow(x.test)), ngenes, paste(genesf, collapse="-"), lambda)
	final_matrix <- matrix(final_list, ncol=15, byrow=TRUE)
	colnames(final_matrix) <- c("model id", "TN", "FP", "FN", "TP", "Sensitivity", "Specificity", 
	                          "PPV", "Error", "CI_Low", "CI_High", "Test-set N", "no_of_genes","genes", "lambda")
	
	if(save==TRUE){
	prefix=prefix
		fname <- paste(path, prefix, "_", id, sep="")  
		saveRDS(x.test, paste(fname, "_xtest_", ".rds", sep=""))
		saveRDS(cvfit, paste(fname, "_cvfit_", ".rds", sep="")) 
	}
	return(final_matrix) #Return the final metrics matrix
	}

cross_val <- function(imp_genes, sf, cd, kk, lm, save=FALSE, prefix, id, path, times){
	#Function that runs the 'lasso_prediction' function for N=times times and averages the metrics results.
	#Arguments:
		#All the arguments of this function except of the 'times' argument correspond to the arguments of the
		#'lasso_prediction' function: imp_genes=genes, sf=split_factor, cd=ci, kk=ki, lm=lambda_method, save=save,
		#                             prefix=prefix, id=id, path=path
		#times: Real Number: Regulates the rounds of running the lasso_prediction function.
	#The function returns the averaged metrics table.
	tn <- c()
	fn <- c()
	tp <- c()
	fp <- c()
  whgenes <- c()
	sens <- c()
	spec <- c()
	ppv <- c()
	genes <- c()
	
	for(i in 1:times){
		las_pr <- lasso_prediction(genes=imp_genes, split_factor=sf, ci=cd, ki=kk, lambda_method=lm, 
		                           save=FALSE, prefix=id, id=id, path)
		
		tn[i] <- as.numeric(las_pr[,2])
		fn[i] <- as.numeric(las_pr[,4])
		tp[i] <- as.numeric(las_pr[,5])
		fp[i] <- as.numeric(las_pr[,3])
		sens[i] <- as.numeric(las_pr[,6])
		spec[i] <- as.numeric(las_pr[,7])
		ppv[i] <- as.numeric(las_pr[,8])
		genes[i] <- as.numeric(las_pr[,13])
    whgenes[i] <- list(strsplit(las_pr[,14], split="-")) 		
	}

whgenes <- names(sort(table(unlist(whgenes)),descending=TRUE))
fm2 <- c(sf, cd, kk, lm, mean(tn),mean(fp),mean(fn),mean(tp), mean(sens), mean(spec), mean(ppv), mean(genes), paste(whgenes, collapse=", "))
return(fm2)
}

combine_parameters <- function(imp_genez, prefixz, cscale, kscale, sf_range, timez){
	#Function that runs the 'cross_val' for every possible combination of splitting 
	#factors and training-set upscaling factors and returns a matrix of the prediction
	#metrics results of each round.
	#Arguments:
		#imp_genez: 'important_genes' object (output of the decide_important function).
		#prefixz: prefix to pass to the 'lasso_prediction' prefix argument.
		#sfrange: vector. list of numbers from 0 to 1. Each of them is being passed 
		#         to the cross_val functions's 'sf' argument.
		#cscale: vector. list of numbers. Each of them is being passed to 
		#        the 'cross_val' function's 'cd' argument.
		#kscale: vector. list of numbers. Each of them is being passed to 
		#        the 'cross_val' function's 'kk' argument.
		#lm: could be "min", "1se" or a number. Is being passed to the 
		#    'cross_val' function's 'lm' argument.
		#times: Number. Is being passed to the 'cross_val' function's 'times' argument.
	
	flist <- c()
	
	for(i in sf_range){
	
		for(j in cscale){
		  print(c(i,j,"actual","min"))
			cv_kac1 <- cross_val(imp_genes=imp_genez, sf=i, cd=j, kk="actual", lm="min", save=FALSE, prefix=prefixz, 
			                     id=1, times=timez)
      print(c(i,j,"actual","1se"))
			cv_kac2 <- cross_val(imp_genes=imp_genez, sf=i, cd=j, kk="actual", lm="1se", save=FALSE, prefix=prefixz, 
			                     id=1, times=timez)
			for(l in kscale){
	      print(c(i,"actual",l,"min"))
				cv_cac1 <- cross_val(imp_genes=imp_genez, sf=i, cd="actual", kk=l, lm="min", save=FALSE, prefix=prefixz, 
				                     id=1, times=timez)
        print(c(i,"actual",l,"1se"))
				cv_cac2 <- cross_val(imp_genes=imp_genez, sf=i, cd="actual", kk=l, lm="1se", save=FALSE, prefix=prefixz, 
				                     id=1, times=timez)
        print(c(i,j,l,"min"))
				cv1 <- cross_val(imp_genes=imp_genez, sf=i, cd=j, kk=l, lm="min", save=FALSE, prefix=prefixz, id=1, 
				                 times=timez)
        print(c(i,j,l,"1se"))
				cv2 <- cross_val(imp_genes=imp_genez, sf=i, cd=j, kk=l, lm="1se", save=FALSE, prefix=prefixz, id=1, 
				                 times=timez)
				flist <- c(flist, cv_cac1, cv_cac2, cv_kac1, cv_kac2, cv1, cv2)
			}
		}
	}
				
	fmat <- matrix(flist, ncol=13, byrow=TRUE)
	colnames(fmat) <- c("split_factor", "controls_scale_up", "cases_scale_up", "lambda_method", "TN", "FP", "FN", "TP", 
	                    "Sensitivity", "Specificity", "PPV", "no_of_genes", "genes")
	return(fmat)
}


run_the_cv <- function(par_comb, prefix, path){
	#Function that takes the output of the 'combine_parameters' function as its 
	#argument and it writes it into a CSV file. It also extracts the best model
    #of that object (maximum Sensitivity and PPV) and returns the prediction 
	#metrics of that model and its splitting and upscaling parameters.
	#Arguments:
		#par_comb: Output of the combine_parameters function.
		#prefix: String. It is being used to unify the name of the output CSV file.
		#path: path to write the ouotput CSV file.
	
	write.csv(par_comb, file=paste(path, "/combine_param_", prefix, ".csv", sep=""))
	
	#Extract the best model and its parameters
	bestppvs <- head(subset(par_comb[order(-as.numeric(par_comb[,11])),], par_comb[,9]>0.5)) 
	sp_factor <- as.numeric(bestppvs[which.max(bestppvs[,9]),][1]) 
	cscale <- bestppvs[which.max(bestppvs[,9]),][2]
	kscale <- bestppvs[which.max(bestppvs[,9]),][3]
	
	return(list(bestppvs, sp_factor, cscale, kscale)) #return the prediction metrics and the 
	                                                  #extracted parameters.
}

cross_val_lambdas <- function(genes, split_factor, ci, ki, lambda_method, save=FALSE, prefix="models_1806",
                              id=0, times, path, name){
	#Function that runs the 'lasso_prediction' function for N=times times and averages the metrics results.
	#Arguments:
		#All the arguments of this function except of the 'times' argument correspond to the arguments of the
		#'lasso_prediction' function: genes=genes, split_factor=split_factor, ci=ci, ki=ki, lambda_method=lambda_method, save=save,
		#                             prefix=prefix, id=id, path=path
		#times: Real Number: Regulates the rounds of running the lasso_prediction function.
	#The function returns a matrix of the prediction metrics of all the times 
	#that the 'lasso_prediction' function ran.
	
	flist <- c()
	for(i in 1:times){
	
		las_pr <- lasso_prediction(genes=genes, split_factor=split_factor, ci=ci, ki=ki, save=FALSE, 
		                           lambda_method=lambda_method, prefix=prefix, id=i, path=path)
		write.table(las_pr, file = paste(path, name,  sep = ""), sep=",", col.names = FALSE, append=TRUE)
		flist <- c(flist, las_pr)
		
		}
	
	fmat <- matrix(flist, ncol=15, byrow=TRUE)
	colnames(fmat) <- c("model id", "TN", "FP", "FN", "TP", "Sensitivity", "Specificity", "PPV", "Error", 
	                    "CI_Low", "CI_High", "Test-set N", "no_of_genes", "genes", "lambda")
	return(fmat)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

get_the_model_lambdas <- function(lamcv){
	#Function that takes the output of the 'cross_val_lambdas' function as its argument
	#and returns a named vector consisting of the lambda giving the best predictive model
	#(maximising PPV and sensitivity) and the average lambda of all the rounds of the cross
	#validation.
	#Arguments:
		#lamcv: Output of the 'cross_val_lambdas' function
		
	avglm <- mean(as.numeric(lamcv[,"lambda"])) #Average lambda
	
	ppv_ord <- lamcv[order(-as.numeric(lamcv[,"PPV"])),] #Sort the matrix in descending PPV
	                                                     #order
	best10 <- head(subset(ppv_ord, as.numeric(ppv_ord[,"Sensitivity"])>0.5), 10) #Maximise Sensitivity
	
	bestlm <- best10[which.max(as.numeric(best10[,"Sensitivity"])),]["lambda"] #Lambda giving the best
																		       #prediction accuracy
		
	lms <- c(avglm, bestlm)
	names(lms) <- c("Mean Lambda", "Best Lambda")
	
	return(lms) #Return the named vector
}
