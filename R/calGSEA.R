#------------------------------------
#' A Pathways Identification Function
#'
#' This function identifies enriched pathways for your dataset.
#' @param calExprOut The output of calExpr() function.
#' @param calDEGOut The output of calDEG() function.
#' @param DEGthr A vector contains cutoffs for DEGs including T scores and p value. Defaults 0 to t scores and 0.01 to p values.
#' @param Sample_1 A vector contains samples names for Group 1.
#' @param Sample_2 A vector contains samples names for Group 2.
#' @param OutFile The folder to save the output of GSEA analysis.
#' @param GeneSet The gene sets in gmt format. Defaults to MSigDB C2_all_v5_2_symbols.gmt
#' @keywords GSEA analysis.
#' @export 
#' @examples 
#' calGSEA(geneExp, DEG_df, DEGthr = c(0, 0.01), Sample_1 = c("s1", "s2", "s3"), Sample_2 = c("s4", "s5", "s6"), OutFile = "../GSEA/", GeneSet = "C2")

calGSEA <- function(calExprOut, calDEGOut, DEGthr = c(0, 0.01), Sample_1, Sample_2, OutFile, GeneSet = "C2"){
	expr <- calExprOut
	deg_df <- calDEGOut
	tscore_cut <- DEGthr[1]
	pval_cut <- DEGthr[2]
	sam1 <- as.vector(Sample_1)
	sam2 <- as.vector(Sample_2)
	my_out <- OutFile
	my_genSet <- GeneSet
	if(my_genSet == "C2"){
		my_genSet <- system.file("extdata", "c2.all.v5.2.symbols.gmt", package = "BioMethyl")
	}
	
	#----select genes for GSEA
	mygen <- rownames(subset(deg_df, Tscore > tscore_cut | Tscore < (-1) * tscore_cut | Pvalue < pval_cut))
	expr <- expr[mygen, ]

	dat1 <- expr[, sam1]
	dat2 <- expr[, sam2]

	mydat <- cbind(dat1, dat2)

	myoutf1 <- paste(my_out, "GSEA.gct", sep = "")
	myoutf2 <- paste(my_out, "GSEA.cls", sep = "")

	conOut <- file(myoutf1, "w")
	curLine <- "#1.2"
	writeLines(curLine, conOut)
	curLine <- paste(nrow(mydat), ncol(mydat), sep = "\t")
	writeLines(curLine, conOut)
	curLine <- paste(c("NAME", "DESCRIPTION", colnames(mydat)), collapse = "\t")
	writeLines(curLine, conOut)
	for(k in 1:nrow(mydat))
	{
		#cat("\r", k)
		curLine <- c(row.names(mydat)[k], "na", mydat[k, ])
		curLine <- paste(curLine, collapse = "\t")
		writeLines(curLine, conOut)	
	}
	close(conOut)

	conOut <- file(myoutf2, "w")
	curLine <- paste(c(ncol(mydat), 2, 1), collapse=" ")
	writeLines(curLine, conOut)
	curLine <- "# Group1 Group2"
	writeLines(curLine, conOut)
	curLine <- c(rep("Group1", ncol(dat1)), rep("Group2", ncol(dat2)))
	curLine <- paste(curLine, collapse = " ")
	writeLines(curLine, conOut)
	close(conOut)
	
	#system.file("extdata", "GSEA.1.0.R", package = "BioMethyl")
	#GSEA.program.location <- system.file("extdata", "GSEA.1.0.R", package = "BioMethyl")
	
	source <- as.function(source)
	source(system.file("extdata", "GSEA.1.0.R", package = "BioMethyl"), verbose=F, max.deparse.length=9999)
	
	if(length(sam1) < 10 & length(sam2) < 10){
		my_type <- "gene.labels"
	}else{
		my_type <- "sample.labels"
	}
	
	xx <- GSEA(input.ds = myoutf1, input.cls = myoutf2, gs.db = my_genSet, output.directory = my_out, doc.string = "GSEA.analysis",
	non.interactive.run = F, reshuffling.type = my_type, nperm = 1000, weighted.score.type = 0, nom.p.val.threshold = -1, fwer.p.val.threshold = -1,
	adjust.FDR.q.val = F, gs.size.threshold.min = 25, gs.size.threshold.max = 500, reverse.sign = T, preproc.type = 0, random.seed = 1111,
	perm.type = 1, fraction = 1.0, replace = F, save.intermediate.results = F, OLD.GSEA = F, use.fast.enrichment.routine = T)
	return (xx)
}
