## Clear variables
rm(list = ls())

## Initialization

setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy

FileName <- c("Xena_TCGA_PanCancer_GE.xena")
FileName2 <- c("Xena_TCGA_PanCancer_Pheno")
FileName3 <- c("Xena_TCGA_PanCancer_Pheno2.tsv")

Target_gene_name_Multi <- c("RHOA_related_StageI&II(ajcc)")
Target_gene_name <- c("RHOA","ROCK1","ROCK2","FN1")
SetVersion <- c("_20210518V1")

ResultFolderName <- paste0("/",Target_gene_name_Multi,SetVersion) ## Generate output folder automatically
dir.create(paste0(PathName,ResultFolderName))

## Import genetic data file
GeneExp_Ori <- read.delim(paste0(PathName,"/",FileName),  # 資料檔名
                          header=T,          # 資料中的第一列，作為欄位名稱
                          sep="\t")
# GeneExp_Ori <- read.table(paste0(PathName,"/Xena_TCGA_LGG_GE"),  # 資料檔名
#                           header=F,          # 資料中的第一列，作為欄位名稱
#                           sep="")

## Import phenotype data file
Phenotype_Ori <- read.delim(paste0(PathName,"/",FileName2),  # 資料檔名
                            header=T,          # 資料中的第一列，作為欄位名稱
                            sep="\t")

Phenotype_Ori2 <- read.delim(paste0(PathName,"/",FileName3),  # 資料檔名
                            header=T,          # 資料中的第一列，作為欄位名稱
                            sep="\t")

###################################### Note ######################################
# paste0 ==> concatenate strings without any separation/delimiter
# paste("Hello", "World", sep = "-") ==> concatenate strings with seperator "-"
##################################################################################


GeneExp <- GeneExp_Ori 
GeneExp[,1] <- as.character(GeneExp[,1])

GeneExp <- GeneExp[!duplicated(GeneExp$sample),]
# GeneExp[!duplicated(GeneExp[,c('sample')]),] 
# GeneExp[!duplicated(GeneExp[,1]),] 
# library("data.table") 
# setDT(GeneExp)[, .SD[1], by = .(sample)] 

# unique(GeneExp, by=c("sample"))

row.names(GeneExp) <- GeneExp[,1]
#GeneExp <- GeneExp[1:length(GeneExp[,1]), 2:length(GeneExp[1,])]
GeneExp <- GeneExp[, -1]

# colnames(GeneExp) <- GeneExp[1,]
# GeneExp <- GeneExp[-1, ]

######## Filter Data by Phenotype ########
Phenotype_Stage <-  Phenotype_Ori[Phenotype_Ori$clinical_stage %in% c("Stage I","Stage II"),]
Phenotype_Stage <-  Phenotype_Ori[Phenotype_Ori$ajcc_pathologic_tumor_stage %in% c("Stage I","Stage II"),]
Phenotype_Stage2 <- Phenotype_Stage
Phenotype_Stage2$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage2$X_PATIENT)
Phenotype_Stage2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage2$sample)

Phenotype_Stage_All <- Phenotype_Ori
# Try stage
Phenotype_Stage_All <- Phenotype_Stage2 
Phenotype_Stage_All$X_PATIENT <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage_All$X_PATIENT)
Phenotype_Stage_All$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_Stage_All$sample)


PT = "Primary Tumor"
STN = "Solid Tissue Normal"
Phenotype_PT <-  Phenotype_Ori2[Phenotype_Ori2$sample_type %in% PT,]
Phenotype_PT2 <- Phenotype_PT
Phenotype_PT2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_PT2$sample)

Phenotype_STN <-  Phenotype_Ori2[Phenotype_Ori2$sample_type %in% STN,]
Phenotype_STN2 <- Phenotype_STN
Phenotype_STN2$sample <- gsub(pattern = '-',replacement = '.',x = Phenotype_STN2$sample)

#
GeneExp_Stage_PT <-  GeneExp[,colnames(GeneExp) %in% Phenotype_PT2$sample]
GeneExp_Stage_STN <-  GeneExp[,colnames(GeneExp) %in% Phenotype_STN2$sample]
# GeneExp_Stage_STN2 <-  GeneExp_Stage_STN[,colnames(GeneExp_Stage_STN) %in% Phenotype_Stage2$sample]

# PSList_Test <- intersect(Phenotype_PT2$sample_type_id,Phenotype_STN2$sample_type_id)
# PSList <- intersect(Phenotype_PT2$sample_type_id,Phenotype_STN2$sample_type_id)

## 
Phenotype_STN_fromGE <- Phenotype_Stage_All[Phenotype_Stage_All$sample %in% colnames(GeneExp_Stage_STN),]
Phenotype_PT_fromGE <- Phenotype_Stage_All[Phenotype_Stage_All$sample %in% colnames(GeneExp_Stage_PT),]

PNList <- intersect(Phenotype_PT_fromGE$X_PATIENT,Phenotype_STN_fromGE$X_PATIENT)

Phenotype_STN_fromGE_intersect <- Phenotype_STN_fromGE[Phenotype_STN_fromGE$X_PATIENT %in% PNList,]
Phenotype_PT_fromGE_intersect <- Phenotype_PT_fromGE[Phenotype_PT_fromGE$X_PATIENT %in% PNList,]
Phenotype_Stage_All_PNintersect <- Phenotype_Stage_All[Phenotype_Stage_All$X_PATIENT %in% PNList, ]
  
# Final Normal and PrimaryTumor
GeneExp_Stage_PT2 <-  GeneExp_Stage_PT[,colnames(GeneExp_Stage_PT) %in% Phenotype_PT_fromGE_intersect$sample]
GeneExp_Stage_STN2 <-  GeneExp_Stage_STN[,colnames(GeneExp_Stage_STN) %in% Phenotype_STN_fromGE_intersect$sample]

###### Try ######
# library("dplyr")
# Phenotype_PT3 <- left_join(Phenotype_PT2,Phenotype_Stage_All,by="sample")
# Phenotype_STN3 <- left_join(Phenotype_STN2,Phenotype_Stage_All,by="sample")
# 
# PSList <- as.data.frame(intersect(Phenotype_PT3$X_PATIENT,Phenotype_STN3$X_PATIENT))
# #PSList <- intersect(Phenotype_PT3$X_PATIENT,Phenotype_STN3$X_PATIENT)
# colnames(PSList) <-c("X_PATIENT")
# 
# PSList2 <- Phenotype_Stage_All[Phenotype_Stage_All$X_PATIENT %in% PSList$X_PATIENT,]
# 
# Phenotype_PT4 <- Phenotype_PT3[Phenotype_PT3$X_PATIENT %in% PSList2$X_PATIENT,]
# Phenotype_STN4 <- Phenotype_STN3[Phenotype_STN3$X_PATIENT %in% PSList2$X_PATIENT,]
# 
# 
# GeneExp_Stage_PT2 <-  GeneExp_Stage_PT[,colnames(GeneExp_Stage_PT) %in% Phenotype_PT4$sample]
# GeneExp_Stage_STN2 <- GeneExp_Stage_STN[,colnames(GeneExp_Stage_STN) %in% Phenotype_STN4$sample]

##########################################


## Find Gene

GeneExp_RhoA_PT <- GeneExp_Stage_PT2[rownames(GeneExp_Stage_PT2) == "RHOA",]
GeneExp_RhoA_PT[2,] <- colnames(GeneExp_RhoA_PT)
GeneExp_RhoA_PT <- as.data.frame(t(GeneExp_RhoA_PT))
colnames(GeneExp_RhoA_PT)[2] <- c("sample")

GeneExp_RhoA_STN <- GeneExp_Stage_STN2[rownames(GeneExp_Stage_STN2) == "RHOA",]
GeneExp_RhoA_STN[2,] <- colnames(GeneExp_RhoA_STN)
GeneExp_RhoA_STN <- as.data.frame(t(GeneExp_RhoA_STN))
colnames(GeneExp_RhoA_STN)[2] <- c("sample")

library("dplyr")
GeneExp_RhoA_PT <- left_join(GeneExp_RhoA_PT,Phenotype_Stage_All_PNintersect,by="sample")
GeneExp_RhoA_PT <- GeneExp_RhoA_PT[,1:3]
GeneExp_RhoA_STN <- left_join(GeneExp_RhoA_STN,Phenotype_Stage_All_PNintersect,by="sample")
GeneExp_RhoA_STN <- GeneExp_RhoA_STN[,1:3]
GeneExp_RhoA_All <- left_join(GeneExp_RhoA_PT,GeneExp_RhoA_STN,by="X_PATIENT")
GeneExp_RhoA_All$RHOA.x <- as.numeric(GeneExp_RhoA_All$RHOA.x)
GeneExp_RhoA_All$RHOA.y <- as.numeric(GeneExp_RhoA_All$RHOA.y)

GeneExp_RhoA_All_Candidates <- GeneExp_RhoA_All[GeneExp_RhoA_All$RHOA.x > GeneExp_RhoA_All$RHOA.y,]

### Test ###
# library("dplyr")
# GeneExp_RhoA_2Group <- full_join(GeneExp_RhoA_PT,GeneExp_RhoA_STN,by="sample")
# colnames(GeneExp_RhoA_2Group)[1] <- c("RhoA.PT")
# colnames(GeneExp_RhoA_2Group)[1] <- c("RhoA.STN")
# 
# GeneExp_RhoA_All <- left_join(GeneExp_RhoA_2Group,Phenotype_Stage_All_PNintersect,by="sample")
# # ttt <- intersect(colnames(GeneExp_Stage_PT2),colnames(GeneExp_Stage_STN2))

GeneExp_Stage_PT_Final <-  GeneExp_Stage_PT2[,colnames(GeneExp_Stage_PT2) %in% GeneExp_RhoA_All_Candidates$sample.x]
GeneExp_Stage_STN_Final <-  GeneExp_Stage_STN2[,colnames(GeneExp_Stage_STN2) %in% GeneExp_RhoA_All_Candidates$sample.y]

######## Old version #########
# 
# 
# # load package 'data.table' 
# library(data.table)
# 
# # Extract data with Target_gene_name
# Target_gene <- GeneExp[Target_gene_name,]
# 
# # load package 'dplyr'
# library(dplyr) # Basic data manupilation tools
# Target_gene_Mean <- rowMeans(data.matrix(Target_gene))
# # Target_gene_SD <- sd(data.matrix(Target_gene))
# library(matrixStats)
# Target_gene_SD <- rowSds(data.matrix(Target_gene))
# 
# # GeneExp_High <- GeneExp[,GeneExp[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD]
# # GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD]
# 
# 
# # Test OK
# # TTT <- as.data.frame(Target_gene_Mean+Target_gene_SD)
# # TTT2 <- as.data.frame(Target_gene_Mean-Target_gene_SD)
# # GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= TTT[Target_gene_name[1],]]
# # GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name[1],] <= TTT2[Target_gene_name[1],]]
# 
# MpSD <- as.data.frame(Target_gene_Mean+Target_gene_SD)
# MmSD <- as.data.frame(Target_gene_Mean-Target_gene_SD)
# Mean <- as.data.frame(Target_gene_Mean)
# 
# # GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= MpSD[Target_gene_name[1],] 
# #                         & GeneExp[Target_gene_name[2],] >= MpSD[Target_gene_name[2],] 
# #                         & GeneExp[Target_gene_name[3],] >= MpSD[Target_gene_name[3],]
# #                         & GeneExp[Target_gene_name[4],] >= MpSD[Target_gene_name[4],]]
# # GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name[1],] <= MmSD[Target_gene_name[1],]]
# 
# GeneExp_High <- GeneExp[,GeneExp[Target_gene_name[1],] >= Mean[Target_gene_name[1],] 
#                         & GeneExp[Target_gene_name[2],] >= Mean[Target_gene_name[2],] 
#                         & GeneExp[Target_gene_name[3],] >= Mean[Target_gene_name[3],]
#                         & GeneExp[Target_gene_name[4],] >= Mean[Target_gene_name[4],]]
# 
# GeneExp_Low<- GeneExp[,GeneExp[Target_gene_name[1],] <= Mean[Target_gene_name[1],] 
#                         & GeneExp[Target_gene_name[2],] <= Mean[Target_gene_name[2],] 
#                         & GeneExp[Target_gene_name[3],] <= Mean[Target_gene_name[3],]
#                         & GeneExp[Target_gene_name[4],] <= Mean[Target_gene_name[4],]]
######## Old version #########

GeneExp_High <- GeneExp_Stage_PT_Final 
GeneExp_Low <- GeneExp_Stage_STN_Final 


# Count the numbers
GeneExp_High_Num <- length(GeneExp_High[1,])
GeneExp_Low_Num <- length(GeneExp_Low[1,])
GeneExp_Gene_Num <- length(GeneExp[,1])
Sample_Num <- GeneExp_High_Num+GeneExp_Low_Num



GeneExp_High_Group <- matrix(c(0), nrow =1, ncol =GeneExp_High_Num)
GeneExp_Low_Group <- matrix(c(1), nrow =1, ncol =GeneExp_Low_Num)


GeneExp_Sum <- cbind(NAME=row.names(GeneExp),Description=matrix(c("na"), nrow =GeneExp_Gene_Num, ncol =1),GeneExp_High,GeneExp_Low)


GSEA_GeneExp <- data.frame(t(colnames(GeneExp_Sum)), stringsAsFactors=FALSE)
colnames(GSEA_GeneExp) <- GSEA_GeneExp

#https://blog.csdn.net/sinat_26917383/article/details/50676894
library("plyr")  
GSEA_GeneExp <- rbind.fill(GSEA_GeneExp,GeneExp_Sum)

#########################################


GSEASetting <- data.frame(NAME = c("#1.2",GeneExp_Gene_Num),Description = c('',Sample_Num))
GSEA_GeneExp<-rbind.fill(GSEASetting,GSEA_GeneExp)

write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.csv"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = ',')
write.table(GSEA_GeneExp,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,"_collapsed.gct"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = '\t')
#########################################


Pheno_Line1 <- c(Sample_Num,2,1)
Pheno_Line2 <- c(paste0("#",Target_gene_name_Multi,"_PT"),paste0(Target_gene_name_Multi,"_STN"))
Pheno_Line3 <- c(GeneExp_High_Group,GeneExp_Low_Group)
Pheno_sum <- rbind.fill(data.frame(t(Pheno_Line1)),data.frame(t(Pheno_Line2), stringsAsFactors=FALSE),data.frame(t(Pheno_Line3)))
write.table(Pheno_sum,file=paste0(PathName,ResultFolderName,"/",FileName,"_",Target_gene_name_Multi,".cls"),quote = FALSE,row.names = FALSE, na = "",col.names = FALSE)
#########################################