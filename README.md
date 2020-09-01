# BioMethyl R package

## Introduction
[BioMethyl](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz137/5364266) is an R package for Biological interpretation of DNA methylation data in the TCGA cancers context. 
For each cancer type, we trained linear regression models for each gene to calculate the associations between its gene expression and its covered CpG sites and recorded the associations (weights) in our BioMethyl.
When using new DNA methylation data as input, the BioMethyl infers samples' gene expression which can be applied to gene set enrichment analysis (GSEA) and gene set analysis (Fisher's exact test).
 

## Installation
1. Clone the package from GitHub using git clone https://github.com/yuewangpanda/BioMethyl.git package_path. An alternative way to download the package is from [here](https://morgan1.dartmouth.edu/~f002nfh/BioMethyl/). <br/>
If you are familoar with R, you can download the R scripts from the R/ folder above and run the example with your own data. <br/>
2. Open R and install the package using ```install.packages("package_path/BioMethyl_1.1.tar.gz", repos = NULL)```.<br/>
3. ```library("BioMethyl")```.<br/>

## Example data
DNA methylation data of 10 ER+ and 10 ER- TCGA breast cancer samples are provided as example data in the package. 
```{r}
data(MethData)
```
## BioMethyl Usage
(1) **DNA methylation data pre-process.**
BioMethyl removes CpG sites that have missing values in more than half samples and imputes the rest missing values by integrating the [ENmix](https://bioconductor.org/packages/release/bioc/html/ENmix.html) R package, with default parameters.
```{r}
dat <- filterMethyData(MethData)
```
(2) **Gene expression calculation.**
BioMethyl calculates the gene expression from the processed DNA methylation data with our pre-trained weights.
The first parameter is the processed DNA methylation data.
The second parameter is the [abbreviation](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) of TCGA cancer type. Due to the sample size limitations, BioMethy now works for ACC, BLCA, CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBM, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS and UVM. If you have the DNA methylation for a non-cancer disease, set this parameter as "NonCancer", the non-cancer model trained with TCGA normal samples will be called.
The third parameter is whether this is an example. Default setting is FALSE.
The fourth parameter is whether you need to save the expression profile in txt format. Defaults to FALSE. 
If you want to a text file as output, the last parameter Path_of_output is needed.
```{r}
myExpr <- calExpr(dat, "BRCA", Example = T, SaveOut = F, "")
```
(3) **Differentially expressed genes identification.** 
BioMethyl identifies DEGs from the inferred gene expression for downstream analysis. T-test was applied for identification.
The first parameter is the gene expression matirx.
The second and third parameters are vectors contains samples names of two groups for comparison.
The fourth parameter is whether this is an example. Default setting is FALSE.
The fourth parameter is whether you need to save the DEGs result as a text file. Defaults to FALSE. 
If you want to a text file as output, the last parameter Path_of_output is needed.
```{r}
samp_1 <- colnames(myExpr)[1:10]
samp_2 <- colnames(myExpr)[11:20]
mydf <- calDEG(myExpr, Sample_1 = samp_1, Sample_2 = samp_2, SaveOut = F, "")
```
(4) **GSEA analysis.** 
BioMethyl integrates [GSEA R package](http://software.broadinstitute.org/gsea/downloads.jsp) to conduct pathway enrichment analysis. 
GMT format gene set is needed for BioMethyl package.
The first parameter is the gene expression matirx.
The second parameter is DEG result.
The third parameter is a vector contains cutoffs for DEGs including T scores and p value. Defaults 0 to t scores and 0.01 to p values.
The fourth and fifth parameter samples names for two groups.
The sixth parameter is path to save the output of GSEA analysis.
The seventh parameter is the gene sets in gmt format. Defaults to MSigDB C2_all_v5_2_symbols.gmt, accessing by set to "C2".
```{r}
xx <- calGSEA(myExpr, mydf, DEGthr = c(0, 0.01), Sample_1 = samp_1, Sample_2 = samp_2, OutFile = "", GeneSet = "C2")
```
(5) **Cancer type recommendation.** 
By applying a centroid method, BioMethyl suggests a suitable cancer type model having the best similarity with TCGA cancers for a DNA methylation data when it is not clear. 
```{r}
data(cancer_centroid)
referCancerType(dat)
```
## Contact
Please contact Yue Wang at yue.wang@dartmouth.edu or Chao Cheng at Chao.Cheng@dartmouth.edu if you have any questions.
