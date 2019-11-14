#---------------------------------
#' A Model Recommendation Function
#'
#' This function gives you an appropriate model according to your DNA methylation data.
#' @param DNA methylation dataset after filtering.
#' @keywords Model selection
#' @export
#' @examples
#' referCancerType()

referCancerType <- function(filterMethyData_out){
	user_df <- filterMethyData_out
	can <- corr <- NULL
	for(i in 1:length(names(cancer_centroid))){
		cancer_cpg <- names(cancer_centroid[[i]])
		my_cpg <- rownames(user_df)
		com_cpg <- intersect(cancer_cpg, my_cpg)

		cancer_cen <- cancer_centroid[[i]][com_cpg]
		my_cpg <- user_df[names(cancer_cen), ]
		my_cpg <- apply(my_cpg, 1, mean, na.rm = T)
		
		mycor <- cor(cancer_cen, my_cpg, method = "s")
		if(mycor == "NA"){
			mycor <- 0
		}else{
			mycor <- round(mycor, 4)
		}
		mypr <- paste("The correlation between your methylation data and", names(cancer_centroid)[i], "is", mycor, sep = " ")
		print (mypr)
		can <- c(can, names(cancer_centroid)[i])
		corr <- c(corr, mycor)
	}
	res <- max(corr)
	if (res < 0.3){
		mypr <- paste("The", can[which.max(corr)], " model is the best choice. But you still need careful consideration!!!", sep = " ")
		print (mypr)
	}else{
		mypr <- paste("The", can[which.max(corr)], " model is the best choice!!!", sep = " ")
		print (mypr)
	}
}
