rm(list = ls())

## R packages loading
library(openxlsx)
library(tidyverse)
library(psych)
library(reshape2)
library(MetaboAnalystR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(httr)
library(KEGGREST)

## Functions
WSM <- function(DataPath, 
                MetaboData, ProteinData,
                Met_SgnValue = 0.05, Pro_SgnValue = 0.05, SgnType = "pvalue",
                CorValue = 0.6, Cor_p = 0.05,
                Org = "mmu",
                pvalue_CutOff = 1, qvalue_CutOff = 1)
{
# Parameter Settings  
  DataPath <- "C:/Users/sjc/Desktop/WSM"
  MetaboData <- "Simultaneous Metabolomics Data.xlsx"
  ProteinData <- "Simultaneous Proteomics Data.xlsx"
  Met_SgnValue <- 0.05
  Pro_SgnValue <- 0.05
  SgnType <- "pvalue"
  CorValue <- 0.6
  Cor_p <- 0.05
  Org <- "mmu"
  pvalue_CutOff <- 1
  qvalue_CutOff <- 1
  setwd(DataPath)
  
# Recognition of affected metabolites and proteins
# Affected metabolites
  Sim_met_df <- read.xlsx(MetaboData, sheet = 1, rowNames = T)
  group_met <- read.xlsx(MetaboData, sheet = 2, rowNames = T)
  
  pvalue_met <- c()
  log2FC_met <- c()
  
  for (i in 1:nrow(Sim_met_df))
  {
    pvalue_met <- c(pvalue_met,
                    t.test(as.numeric(Sim_met_df[i,which(group_met==unique(group_met$Group)[1])]),
                           as.numeric(Sim_met_df[i,which(group_met==unique(group_met$Group)[2])]),
                           var.equal = T)$p.value)
    log2FC_met <- c(log2FC_met,
                    log2((mean(as.numeric(Sim_met_df[i,which(group_met==unique(group_met$Group)[1])])))/
                         (mean(as.numeric(Sim_met_df[i,which(group_met==unique(group_met$Group)[2])]))))) 
  }
  FDR_met <- p.adjust(pvalue_met, "BH")
  Sim_significance_met <- cbind(Sim_met_df, log2FC_met, pvalue_met, FDR_met)
  if(SgnType == "pvalue")
  {
    Sim_affected_met <- Sim_significance_met[which(Sim_significance_met$pvalue_met<=Met_SgnValue),]
    Affected_metabolites <- t(Sim_met_df[which(pvalue_met<=Met_SgnValue),])
  }else if(SgnType == "FDR")
  {
    Sim_affected_met <- Sim_significance_met[which(Sim_significance_met$FDR_met<=Met_SgnValue),]
    Affected_metabolites <- t(Sim_met_df[which(FDR_met<=Met_SgnValue),])
  }
  write.xlsx(Sim_affected_met, file = "1 Affected Metabolites.xlsx", colNames = T, rowNames = T)
# Affected proteins
  Sim_pro_df <- read.xlsx(ProteinData, sheet = 1,rowNames = T)
  group_pro <- read.xlsx(ProteinData, sheet = 2,rowNames = T)
  
  pvalue_pro <- c()
  log2FC_pro <- c()
  
  for (i in 1:nrow(Sim_pro_df)) 
  {
    pvalue_pro <- c(pvalue_pro,
                    t.test(as.numeric(Sim_pro_df[i,which(group_pro==unique(group_pro$Group)[1])]),
                           as.numeric(Sim_pro_df[i,which(group_pro==unique(group_pro$Group)[2])]),
                           var.equal = T)$p.value)
    log2FC_pro <- c(log2FC_pro,
                    log2((mean(as.numeric(Sim_pro_df[i,which(group_pro==unique(group_pro$Group)[1])])))/
                         (mean(as.numeric(Sim_pro_df[i,which(group_pro==unique(group_pro$Group)[2])])))))
  }
  FDR_pro <- p.adjust(pvalue_pro, "BH")
  Sim_significance_pro <- cbind(Sim_pro_df, log2FC_pro, pvalue_pro, FDR_pro)
  if(SgnType == "pvalue")
  {
    Sim_affected_pro <- Sim_significance_pro[which(Sim_significance_pro$pvalue_pro<=Pro_SgnValue),]
    Affected_proteins <- t(Sim_pro_df[which(pvalue_pro<=Pro_SgnValue),])
  }else if(SgnType == "FDR")
  {
    Sim_affected_pro <- Sim_significance_pro[which(Sim_significance_pro$FDR_pro<=Pro_SgnValue),]
    Affected_proteins <- t(Sim_pro_df[which(FDR_pro<=Pro_SgnValue),])
  }
  write.xlsx(Sim_affected_pro, file = "2 Affected Proteins.xlsx", colNames = T, rowNames =T)
  
# Construction of mathematical correlation  
# Analysis method including: Pearson, Spearman and Kendall
  correlation <- corr.test(Affected_metabolites, Affected_proteins, method = "spearman")
  coefficient <- correlation$r
  pvalue_spearman <- correlation$p
  correlation_matrix <- data.frame(melt(coefficient, value.name = "correlation"),
                                   "pvalue_spearman" = as.vector(pvalue_spearman))
  Mathematical_correlation <- correlation_matrix[which(abs(correlation_matrix$correlation)>=CorValue &
                                                       correlation_matrix$pvalue_spearman<=Cor_p),]
  colnames(Mathematical_correlation) <- c("Metabolites","Proteins","Correlation Coefficient","pvalue_spearman")
  write.xlsx(Mathematical_correlation, file = "3 Mathematical correlation.xlsx")
  
# Construction of biological correlation
# Enrichment analysis of affected metabolites
# Convert metabolite ID to KEGG compound ID
  Affected_met_name <- rownames(Sim_affected_met)
  tosend <- list(queryList = Affected_met_name, inputType = "name")
  call <- "https://www.kegg.jp/kegg/rest/keggapi.html"
  query_metabolites <- httr::POST(call, body = tosend, encode = "json")
  #query_metabolites$status_code # 200 is OK, others mean some errors
  query_metabolites_text <- content(query_metabolites, "text", encoding = "UTF-8")
  query_metabolites_json <- rjson::fromJSON(query_metabolites_text)
  KEGGID <- as.data.frame(matrix(query_metabolites_json$KEGG))
  names(KEGGID)[1] <- "KEGG_ID"
  queryName <- as.vector(query_metabolites_json$Query)
  Affected_met_name <- data.frame(queryName,KEGGID)
  write.xlsx(Affected_met_name,"4 Affected metabolites mapping.xlsx") # Please remove "NA" character string
  kegg_data <- na.omit(Affected_met_name) # Filter metabolites without KEGG ID
# Construction of background data
  kegg_pathway <- as.data.frame(keggList("pathway", Org))
  kegg_pathway$"KEGG_pathway_ID" <- rownames(kegg_pathway)
  names(kegg_pathway)[1] <- "KEGG_pathway_name"
  kegg_pathway <- kegg_pathway[,c(2,1)]
  TERM2NAME_bg <- kegg_pathway
  TERM2NAME <- data.frame() # Filter pathway without metabolites
  for (i in 1:length(TERM2NAME_bg$KEGG_pathway_ID)) 
  {
    pathway_i <- TERM2NAME_bg$KEGG_pathway_ID[i]
    pathway <- keggGet(pathway_i)
    cpb <- as.data.frame(pathway[[1]]$COMPOUND)
    if(length(cpb) == 0)
    {
      TERM2NAME_bg <- TERM2NAME_bg[TERM2NAME_bg$KEGG_pathway_ID != pathway_i,]
    }else
    {
      TERM2NAME <- rbind(TERM2NAME, TERM2NAME_bg[TERM2NAME_bg$KEGG_pathway_ID == pathway_i,])
    }
  }
  saveRDS(TERM2NAME, file = "mmu_pathway.rds")
  
  TERM2NAME <- readRDS("mmu_pathway.rds")
  TERM2meta <- data.frame()
  for (i in 1:length(TERM2NAME$KEGG_pathway_ID)) 
  {
    pathway <- keggGet(TERM2NAME$KEGG_pathway_ID[i])
    cpb <- as.data.frame(pathway[[1]]$COMPOUND)
    cpb$"metabolite_ID" <- rownames(cpb)
    cpb$KEGG_pathway_ID <- pathway[[1]]$ENTRY[1]
    cpb <- cpb[,-1]
    TERM2meta <- rbind(TERM2meta,cpb)
  }
  TERM2meta <- TERM2meta[,c(2,1)]
  MetEnrichment <- enricher(gene = kegg_data$KEGG_ID,
                            TERM2GENE = TERM2meta,
                            TERM2NAME = TERM2NAME,
                            pvalueCutoff = pvalue_CutOff,
                            qvalueCutoff = qvalue_CutOff)
  result_data <- MetEnrichment@result
  names(result_data)[c(3,8)] <- c("MetRatio","MetID")
  write.xlsx(result_data, file = "5 Enrichment result of affected metabolites.xlsx")
# Enrichment analysis of affected proteins
# Convert protein ID to KEGG enzyme ID (convert ENTREZID to SYMBOL) 
  Affected_pro_name <- rownames(Sim_affected_pro)
  query_results_proteins <- bitr(Affected_pro_name, fromType = 'UNIPROT', toType = c('ENTREZID','SYMBOL'),
                                 OrgDb = 'org.Mm.eg.db')
  write.xlsx(query_results_proteins, file = "6 Affected proteins mapping.xlsx")
# Construction of background data  
  kegg_pathway <- enrichKEGG(gene = query_results_proteins$ENTREZID,
                             organism = Org,
                             keyType = "kegg",
                             pAdjustMethod = "fdr",
                             pvalueCutoff = pvalue_CutOff,
                             qvalueCutoff = qvalue_CutOff)
  protein_result <- setReadable(kegg_pathway, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID')
  ProEnrichment <- protein_result@result
  names(ProEnrichment)[c(3,8)] <- c("ProRatio","ProID")
  write.xlsx(ProEnrichment, file = "7 Enrichment result of affected proteins.xlsx")
# Organize data format
  MetEnrichment_result_df <- result_data
  rownames(MetEnrichment_result_df) <- MetEnrichment_result_df$ID
  num_Metab <- paste("Met", 1:max(MetEnrichment_result_df$Count), sep = "")
  Metab <- MetEnrichment_result_df[8]
  Metab <- separate(Metab, MetID, into = num_Metab, sep = "/")
  MetEnrichment_result <- cbind(MetEnrichment_result_df[c(1:7)], Metab)
  Metdf <- MetEnrichment_result[c(1,8:ncol(MetEnrichment_result))]
  long_Metdf <- melt(Metdf, id.vars = "ID")
  long_Metdf <- long_Metdf[-2]
  Met_pathway <- na.omit(long_Metdf)
  names(Met_pathway) <- c("Pathway ID","Metabolite ID")
  write.xlsx(Met_pathway, file = "8 Metabolite matching Pathway.xlsx")
  
  ProEnrichment_result_df <- ProEnrichment
  rownames(ProEnrichment_result_df) <- ProEnrichment_result_df$ID
  num_Pro <- paste("Pro", 1:max(ProEnrichment_result_df$Count), sep = "")
  Protein <- ProEnrichment_result_df[8]
  Protein <- separate(Protein, ProID, into = num_Pro, sep = "/")
  ProEnrichment_result <- cbind(ProEnrichment_result_df[c(1:7)], Protein)
  Proteindf <- ProEnrichment_result[c(1,8:ncol(ProEnrichment_result))]
  long_Prodf <- melt(Proteindf, id.vars = "ID")
  long_Prodf <- long_Prodf[-2]
  Pro_pathway <- na.omit(long_Prodf)
  names(Pro_pathway) <- c("Pathway ID","Protein ID")
  write.xlsx(Pro_pathway, file = "9 Protein matching Pathway.xlsx")
# Integrate the metabolomics KEGG analysis result and proteomics KEGG analysis result
  Biological_correlation <- merge(Met_pathway, Pro_pathway, by = 'Pathway.ID')
  write.xlsx(Biological_correlation, file = "10 Biological correlation.xlsx")
  
# Construction of integrated correlation network
  Metmapping <- Affected_met_name
  Promapping <- query_results_proteins
  kegg_data <- na.omit(Metmapping)
  math_df <- Mathematical_correlation[c(1,2)]
  
  mathcor1 <- data.frame()
  for (m in 1:nrow(math_df)) 
  {
    for (n in 1:nrow(kegg_data)) 
    {
      if(math_df[m,1]==kegg_data[n,1])
      {
        mathcor1_m <- data.frame(math_df[m,], kegg_data[n,])
        mathcor1 <- rbind(mathcor1, mathcor1_m)
      }
    }
  }
  mathcor_met <- mathcor1[,c(4,2)]
  
  mathcor2 <- data.frame()
  for (a in 1:nrow(mathcor_met)) 
  {
    for (b in 1:nrow(Promapping)) 
    {
      if(mathcor_met[a,2]==Promapping[b,1])
      {
        mathcor2_a <- data.frame(mathcor_met[a,], Promapping[b,])
        mathcor2 <- rbind(mathcor2, mathcor2_a)
      }
    }
  }
  mathcor_met_pro <- mathcor2[,c(1,5)]
  names(mathcor_met_pro) <- c("Metabolite.ID", "Protein.ID")
  trans_Mathematical_correlation <- mathcor_met_pro
  
  trans_Mathematical_correlation$combine <- paste(trans_Mathematical_correlation$Metabolite.ID,
                                                  trans_Mathematical_correlation$Protein.ID,
                                                  sep = "_")
  Biological_correlation$combine <- paste(Biological_correlation$Metabolite.ID,
                                          Biological_correlation$Protein.ID,
                                          sep = "_")
  Correlation_Network <- merge(trans_Mathematical_correlation, Biological_correlation, by = "combine")
  Correlation_Network <- Correlation_Network[,c(1:4)]
  write.xlsx(Correlation_Network, file = "11 Correlation Network.xlsx")
}

## Application of WSM function  
WSM(DataPath = "C:/Users/sjc/Desktop/WSM",
    MetaboData = "Simultaneous Metabolomics Data.xlsx",
    ProteinData = "Simultaneous Proteomics Data.xlsx",
    Met_SgnValue = 0.05, Pro_SgnValue = 0.05, SgnType = "pvalue",
    CorValue = 0.6, Cor_p = 0.05,
    Org="mmu",
    pvalue_CutOff = 1, qvalue_CutOff = 1)
