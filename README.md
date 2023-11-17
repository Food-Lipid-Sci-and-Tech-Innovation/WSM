# WSM
# Parameter Description
DataPath: Set the directory containing metabolomics data and proteomics data
MetaboData: Metabolomics data
ProteinData：Proteomics data
Met_SgnValue：The threshold value for significant metabolite
Pro_SgnValue：The threshold value for significant protein
SgnType："pvalue" or "FDR"
CorValue：The threshold value for correlation coefficient
Cor_p：The threshold value for the significance of correlation analysis
Org：KEGG Organisms can be checked in "as.data.frame(keggList("organism"))"
pvalue_CutOff & qvalue_CutOff：The threshold value for KEGG enrichment analysis
# Data Format
Both metabolomics and proteins data were saved with ".xlsx” format where sheet1 contained biomolecular information and sheet2 contained group information.
# Version of Packages
R version 4.3.1
library(openxlsx) : 4.2.5.2
library(tidyverse): 2.0.0
library(psych): 2.3.6
library(reshape2): 1.4.4
library(MetaboAnalystR): 3.0.3
library(clusterProfiler): 4.8.1
library(org.Mm.eg.db): 3.17.0
library(httr): 1.4.6
library(KEGGREST): 1.36.0
