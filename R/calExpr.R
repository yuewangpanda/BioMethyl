#----------------------------------------
#' A Gene expression Calculation Function
#'
#' This function infers gene expression based on DNA methylation data.
#' @param filterMethyData_out The output of filterMethyData() function.
#' @param CancerType Which model?.
#' @param Example Use 500 genes as an example. Defaults to FALSE.
#' @param SaveOut Do you need to save the expression profile in txt format? Defaults to FALSE.
#' @param OutFile The Path_of_output is needed when SaveOut is TRUE.
#' @keywords Gene expression calculation
#' @export 
#' @examples
#' calExpr()

calExpr <- function(filterMethyData_out, CancerType, Example = FALSE, SaveOut = FALSE, OutFile){
	my_meth <- as.matrix(filterMethyData_out)
	my_can <- CancerType
	
	if(my_can == "ACC"){
		myList <- acc_model
	}else if (my_can == "BLCA"){		
		myList <- blca_model
	}else if (my_can == "BRCA"){
		myList <- brca_model
	}else if (my_can == "CESC"){
		myList <- cesc_model
	}else if (my_can == "CHOL"){
		myList <- chol_model
	}else if (my_can == "COAD"){
		myList <- coad_model
	}else if (my_can == "COADREAD"){
		myList <- coadread_model
	}else if (my_can == "DLBC"){
		myList <- dlbc_model
	}else if (my_can == "ESCA"){
		myList <- esca_model
	}else if (my_can == "GBM"){
		myList <- gbm_model
	}else if (my_can == "GBMLGG"){
		myList <- gbmlgg_model
	}else if (my_can == "HNSC"){
		myList <- hnsc_model
	}else if (my_can == "KICH"){
		myList <- kich_model
	}else if (my_can == "KIPAN"){
		myList <- kipan_model
	}else if (my_can == "KIRC"){
		myList <- kirc_model
	}else if (my_can == "KIRP"){
		myList <- kirp_model
	}else if (my_can == "LAML"){
		myList <- laml_model
	}else if (my_can == "LGG"){
		myList <- lgg_model
	}else if (my_can == "LIHC"){
		myList <- lihc_model
	}else if (my_can == "LUAD"){
		myList <- luad_model
	}else if (my_can == "LUSC"){
		myList <- lusc_model
	}else if (my_can == "MESO"){
		myList <- meso_model
	}else if (my_can == "OV"){
		myList <- ov_model
	}else if (my_can == "PAAD"){
		myList <- paad_model
	}else if (my_can == "PCPG"){
		myList <- pcpg_model
	}else if (my_can == "PRAD"){
		myList <- prad_model
	}else if (my_can == "READ"){
		myList <- read_model
	}else if (my_can == "SARC"){
		myList <- sarc_model
	}else if (my_can == "SKCM"){
		myList <- skcm_model
	}else if (my_can == "STAD"){
		myList <- stad_model
	}else if (my_can == "STES"){
		myList <- stes_model
	}else if (my_can == "TGCT"){
		myList <- tgct_model
	}else if (my_can == "THCA"){
		myList <- thca_model
	}else if (my_can == "THYM"){
		myList <- thym_model
	}else if (my_can == "UCEC"){
		myList <- ucec_model
	}else if (my_can == "UCS"){
		myList <- ucs_model
	}else if (my_can == "UVM"){
		myList <- uvm_model
	}else if (my_can == "NonCancer"){
		myList <- normal_model
	}else{
		print ("Please input correct cancer type name!")
	}
	if(Example == T){
		myList <- head(myList, n = 500)
	}
	
	expr_list <- lapply(names(myList), function(gene) {
			tmp <- myList[[gene]]
			coef <- tmp[2:length(tmp)]
			comCpG <- intersect(names(coef), row.names(my_meth))
			meth_dat <- my_meth[comCpG, ]
			coef <- coef[comCpG]
			return (as.vector(coef %*% meth_dat) + tmp[1])
	})
	myres <- do.call(rbind.data.frame, expr_list)
	colnames(myres) <- colnames(my_meth)
	row.names(myres) <- sapply(strsplit(names(myList), "_"), "[", 1)

	if(SaveOut == T){
		my_out <- OutFile
		write.table(myres, file = my_out, col.names = T, row.names = T, sep = "\t", quote = F)
	}
	return(myres)
}
