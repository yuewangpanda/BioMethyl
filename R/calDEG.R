#-----------------------------------------------------------------
#' A Differentially Expressed Genes (DEGs) Identification Function
#'
#' This function identifies DEGs based on the inferred gene expression data.
#' @param calExprOut The output of calExpr() function.
#' @param Sample_1 A vector contains samples names for Group 1.
#' @param Sample_2 A vector contains samples names for Group 2.
#' @param SaveOut Do you need to save DEGs result in txt format? Defaults to FALSE.
#' @param OutFile The Path_of_output is needed when SaveOut is TRUE.
#' @keywords DEGs Identification.
#' @export
#' @examples
#' calDEG()

calDEG <- function(calExprOut, Sample_1, Sample_2, SaveOut = FALSE, OutFile){
	expr <- calExprOut
	sam1 <- as.vector(Sample_1)
	sam2 <- as.vector(Sample_2)
	
	expr_1 <- expr[, colnames(expr) %in% sam1]
	expr_2 <- expr[, colnames(expr) %in% sam2]

	genes <- rownames(expr)

	tscore <- pval <- mean_sam1 <- mean_sam2 <- NULL
	for(gg in genes){
		dat1 <- as.vector(expr_1[gg, ])
		dat2 <- as.vector(expr_2[gg, ])
		len1 <- length(unique(unlist(dat1)))
		len2 <- length(unique(unlist(dat2)))
		if(len1 == 1 & len2 == 1){
			tscore <- c(tscore, 0)
			pval <- c(pval, 1)
			mean_sam1 <- c(mean_sam1, 0)
			mean_sam2 <- c(mean_sam2, 0)
		}else{
			myt <- t.test(dat1, dat2)
			tscore <- c(tscore, myt[["statistic"]])
			pval <- c(pval, myt[["p.value"]])
			mean_sam1 <- c(mean_sam1, myt[["estimate"]][1])
			mean_sam2 <- c(mean_sam2, myt[["estimate"]][2])
		}
	}

	df <- data.frame(Gene = genes, Tscore = tscore, Pvalue = pval, Mean_Sam1 = mean_sam1, Mean_Sam2 = mean_sam2, stringsAsFactors = F)
	rownames(df) <- genes
	df <- df[order(df$Pvalue), ]
	df$FDR <- p.adjust(df$Pvalue, method = "BH")
	if(SaveOut == T){
		my_out <- OutFile
		write.table(myres, file = my_out, col.names = T, sep = "\t", quote = F)
	}
	return (df)
}

