library(data.table)
library(dplyr)
library(limma)  # bioconductor 

## For volcano plot
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel) # CRAN. For "geom_text_repel" function

###########################################################################################
dir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
AMRFile<-"amredit.csv"
###########################################################################################
####################################################################################
## Step1. Read amr data and extract stool samples. Then, split it to Rifaximin and Placebo groups.
####################################################################################
amrData=fread(AMRFile, header=TRUE, stringsAsFactors=FALSE ); dim(amrData); amrData[1:3,1:3]  # 185  583
#             ID    AAC_2_Ib AAC_2_Ic
# 1: P_10_d01_SA        0        0
# 2: P_10_d01_ST        0        0
# 3: P_10_d30_SA        0        0

## Make a new column, "Subtype" to see Rifaximin/Placebo and SA/ST
#sapply(stringr::str_split(amrData$ID, pattern="_"), function(X) {X[1] })  # "P"  or "R"
#sapply(stringr::str_split(amrData$ID, pattern="_"), "]", 1)
RifaPlacebo <- sapply(stringr::str_split(amrData$ID, pattern="_"), "[", 1)  # "P"  or "R"
SAST <-  sapply(stringr::str_split(amrData$ID, pattern="_"), "[", 4)      # "SA" or "ST

amr_Subtype <- amrData %>% dplyr::mutate(Subtype=paste0(RifaPlacebo,"_",SAST ) )
table(amr_Subtype$Subtype)
# P_SA P_ST R_SA R_ST 
# 47   42   48   48
View(amr_Subtype)

############################################################################################
## Step4. Make Fraction barplot - the samples should be ordered by Day
## Plot R code: https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2
############################################################################################
amr_Placebo_ST_Day <- amr_Subtype %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_ST_Day)   # 42  583 by d01, d30, d90

## Make tidyr data format
amr_Placebo_ST_Day_tidyr <- amr_Placebo_ST_Day %>% tidyr::gather(AMRGene, Abundance, AAC_2_Ib:vph);
head(amr_Placebo_ST_Day_tidyr); dim(amr_Placebo_ST_Day_tidyr) # 24444 3
## Remove rows of Abudance==0
amr_Placebo_ST_Day_tidyr_NoZero <- amr_Placebo_ST_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Placebo_ST_Day_tidyr_NoZero) # 4480 3

## Count the frequency of the AMR genes. Some genes appear with low frequency. I should exclude the low-frequency AMRGenes. 
GeneFrequency <- table(amr_Placebo_ST_Day_tidyr_NoZero$AMRGene)
View(GeneFrequency)
length(GeneFrequency)  # there are all 412 AMRGenes.
max(GeneFrequency); min(GeneFrequency) #max 42  # min is 1 which means some AMRGenes appear just once 
GeneFrequency_High <- GeneFrequency[GeneFrequency>35]; length(GeneFrequency_High)  # 15 genes left. I should filter in these genes

## Filter in just 15 AMRGenes of high frequency.
amr_Placebo_ST_Day_tidyr_HigFreq <- amr_Placebo_ST_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
## Check how many kinds of genes are contained
length(table(amr_Placebo_ST_Day_tidyr_HigFreq$AMRGene))  # 15, Ok this small number is good to make fraction barplot.   

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
ID_D30 <- unique(grep("d30", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 13 samples
ID_D90 <- unique(grep("d90", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 11 samples

amr_Placebo_ST_Day_tidyr_HigFreq$ID <- factor(amr_Placebo_ST_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))



MyBarplot <- ggplot(amr_Placebo_ST_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot

MyBarplot_9CellType_Placebo <- MyBarplot + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Placebo group, Stool"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_Placebo

## Save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_Placebo, file="FractionBarplot_Placebo_ST_Day01Day30Day90_W14H6_20220508_v1.0.pdf", width=14, height=6)


####################################  Rifaximin STOOL ############################################################
##########################   Extract only "R_ST" group samples ###################################################
amr_Rifaximin_ST_Day <- amr_Subtype %>% dplyr::filter(grepl("R_ST", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Rifaximin_ST_Day) # 48  583  only R_ST by d01, d30, d90

## Make tidyr data format
amr_Rifaximin_ST_Day_tidyr <- amr_Rifaximin_ST_Day %>% tidyr::gather(AMRGene, Abundance, AAC_2_Ib:vph);
head(amr_Rifaximin_ST_Day_tidyr); dim(amr_Rifaximin_ST_Day_tidyr) # 27936  3
## Remove rows of Abudance==0
amr_Rifaximin_ST_Day_tidyr_NoZero <- amr_Rifaximin_ST_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Rifaximin_ST_Day_tidyr_NoZero) # 4722 3

## Count the frequency of the AMR genes. Some genes appear with low frequency. Exclude the low-frequency AMRGenes. 
GeneFrequency <- table(amr_Rifaximin_ST_Day_tidyr_NoZero$AMRGene)
View(GeneFrequency)
length(GeneFrequency)  # there are all 419 AMRGenes.
max(GeneFrequency); min(GeneFrequency) #max 48  # min is 1 which means some AMRGenes appear just once 
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. I should filter in these genes

## Filter in just 15 AMRGenes of high frequency.
amr_Rifaximin_ST_Day_tidyr_HigFreq <- amr_Rifaximin_ST_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
## Check how many kinds of genes are contained
length(table(amr_Rifaximin_ST_Day_tidyr_HigFreq$AMRGene))  # 14, Ok this small number is good to make fraction barplot.   
#unique(amr_Rifaximin_ST_Day_tidyr_HigFreq$AMRGene)

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
#grep("d01", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)
ID_D30 <- unique(grep("d30", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 16 samples
ID_D90 <- unique(grep("d90", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 14 samples

amr_Rifaximin_ST_Day_tidyr_HigFreq$ID <- factor(amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))

MyBarplot_R<- ggplot(amr_Rifaximin_ST_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot_R

## Change the bar colors manually and adjust the axis text.

MyBarplot_9CellType_Rifa <- MyBarplot_R + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  labs(title="Relative AMRgene abundance in Rifaximin group(sample n=48) over day01, day30, day90")+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Rifaximin group, Stool"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_Rifa

## Save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_Rifa, file="FractionBarplot_Rifa_ST_Day01Day30Day90.pdf", width=14, height=6)


##############################################  SA Placebo ####################################################
View(amr_Subtype)
amr_Placebo_SA_Day <- amr_Subtype %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_SA_Day)   # 47  583 by d01, d30, d90
View(amr_Placebo_SA_Day)
colnames(amr_Placebo_SA_Day)
## Make tidyr data format
amr_Placebo_SA_Day_tidyr <- amr_Placebo_SA_Day %>% tidyr::gather(AMRGene, Abundance, AAC_2_Ib:vph);
head(amr_Placebo_SA_Day_tidyr); dim(amr_Placebo_SA_Day_tidyr) # 27354     3
View(amr_Placebo_SA_Day_tidyr)
## Remove rows of Abudance==0
amr_Placebo_SA_Day_tidyr_NoZero <- amr_Placebo_SA_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Placebo_SA_Day_tidyr_NoZero) # 2979    3
View(amr_Placebo_SA_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency. I should exclude the low-frequency AMRGenes. 
GeneFrequency <- table(amr_Placebo_SA_Day_tidyr_NoZero$AMRGene)
length(GeneFrequency)  # there are all 319 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 47  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. I should filter in these genes
View(GeneFrequency_High)

## Filter in just 15 AMRGenes of high frequency.
amr_Placebo_SA_Day_tidyr_HigFreq <- amr_Placebo_SA_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
## Check how many kinds of genes are contained
length(table(amr_Placebo_SA_Day_tidyr_HigFreq$AMRGene))  # 14, this small number is good to make fraction barplot.   

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Placebo_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 19 samples
ID_D30 <- unique(grep("d30", amr_Placebo_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 16 samples
ID_D90 <- unique(grep("d90", amr_Placebo_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 12 samples

amr_Placebo_SA_Day_tidyr_HigFreq$ID <- factor(amr_Placebo_SA_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))
View(amr_Placebo_SA_Day_tidyr_HigFreq)

MyBarplot_P_SA <- ggplot(amr_Placebo_SA_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 

MyBarplot_P_SA

## Change the bar colors manually and adjust the axis text.

MyBarplot_9CellType_P_SA <- MyBarplot_P_SA + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Placebo group, Saliva"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_P_SA

## Save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_P_SA, file="FractionBarplot_DES_Placebo_SA.pdf", width=14, height=6)

############################################## SA Rifaximin    ##########################################################
View(amr_Subtype)
amr_Rifaximin_SA_Day <- amr_Subtype %>% dplyr::filter(grepl("R_SA", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Rifaximin_SA_Day)   # 48  583 by d01, d30, d90
View(amr_Rifaximin_SA_Day)
colnames(amr_Rifaximin_SA_Day)
## Make tidyr data format
amr_Rifaximin_SA_Day_tidyr <- amr_Rifaximin_SA_Day %>% tidyr::gather(AMRGene, Abundance, AAC_2_Ib:vph);
head(amr_Rifaximin_SA_Day_tidyr); dim(amr_Rifaximin_SA_Day_tidyr) # 27936  3
View(amr_Rifaximin_SA_Day_tidyr)

## Remove rows of Abudance==0
amr_Rifaximin_SA_Day_tidyr_NoZero <- amr_Rifaximin_SA_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Rifaximin_SA_Day_tidyr_NoZero) # 2926  3
View(amr_Rifaximin_SA_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency. Exclude the low-frequency AMRGenes. 
GeneFrequency <- table(amr_Rifaximin_SA_Day_tidyr_NoZero$AMRGene)
View(GeneFrequency)
length(GeneFrequency)  # there are all 302 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 47  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. I should filter in these genes
View(GeneFrequency_High)
## Filter in just 15 AMRGenes of high frequency.
amr_Rifaximin_SA_Day_tidyr_HigFreq <- amr_Rifaximin_SA_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
View(amr_Rifaximin_SA_Day_tidyr_HigFreq)
## Check how many kinds of genes are contained
length(table(amr_Rifaximin_SA_Day_tidyr_HigFreq$AMRGene))  # 14,this small number is good to make fraction barplot.   

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
ID_D30 <- unique(grep("d30", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 16 samples
ID_D90 <- unique(grep("d90", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 14 samples

amr_Rifaximin_SA_Day_tidyr_HigFreq$ID <- factor(amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))
View(amr_Rifaximin_SA_Day_tidyr_HigFreq)
MyBarplot_R_SA <- ggplot(amr_Rifaximin_SA_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot_R_SA

## Change the bar colors manually and adjust the axis text.

MyBarplot_9CellType_R_SA <- MyBarplot_R_SA + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Rifaximin group, Saliva"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_R_SA

## ave the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_R_ST, file="FractionBarplot_DES_Rifaximin_SA.pdf", width=14, height=6)


#############################################################################################################################################