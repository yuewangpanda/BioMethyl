#-------------------------------------------
#' A DNA Methylation Data Filtering Function
#'
#' This function allows you to pre-process your DNA methylation data in two steps.
#' First, removing the cpg sites have na beta values in more than half samples.
#' Second, using k nearest neighbours to replace the rest missing values.
#' @param Your raw DNA methylation dataset.
#' @keywords DNA methylation data
#' @export
#' @examples
#' filterMethyData()

filterMethyData <- function(MethData){
	if(!("ENmix" %in% rownames(installed.packages()))){
		print ("Installation of ENmix R package...")
 		tryCatch(source("https://bioconductor.org/biocLite.R"), error=function(err) source("http://bioconductor.org/biocLite.R"))
 		biocLite("ENmix")
    }
    library("ENmix")
    
    dat <- as.data.frame(MethData)
    tmp <- sapply(1:nrow(dat), function(x) {sum(is.na(dat[x, ]))})	
	cutoff <- round(ncol(dat)/2, 0)
	mydat <- dat[tmp < cutoff, ]
	mydat <- as.matrix(mydat)
	mydat <- ENmix::rm.outlier(mydat, byrow = T, rmcr = T, rthre = 0.5, cthre = 0.5, impute = T, imputebyrow = TRUE)
	return (mydat)
}
