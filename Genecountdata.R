library(data.table)
library(dplyr)
library(DESeq2)  # bioconductor 

## For volcano plot
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel) # CRAN. For "geom_text_repel" function

###########################################################################################
AMRFile<-"nonRelativeAmr.csv"
###########################################################################################
####################################################################################
## Step1. Read amr data and extract stool samples. Then, split it to Rifaximin and Placebo groups.
####################################################################################
amrData=fread(AMRFile, header=TRUE, stringsAsFactors=FALSE ); dim(amrData); amrData[1:3,1:3]  # 185  794
#             ID AAC(2')-IIb AAC(2')-Ib
# 1: P_10_d01_SA        0        0
# 2: P_10_d01_ST        0        0
# 3: P_10_d30_SA        0        0

## Check duplicated gene names
table(duplicated(colnames(amrData)))  # FALSE: 794  There is no duplicated gene symbols


## Remove special characters in column names
colnames(amrData) <- gsub("\\(|\\)","",colnames(amrData))  # 
colnames(amrData) <- gsub("-|\\/|\\.|\\'| ","_",colnames(amrData))

##  Check duplicated gene names again
table(duplicated(colnames(amrData)))  # FALSE: 794  There is no duplicated gene symbols


## Convert decimal numbers to integers
amrData_Integer <- ceiling(amrData[,2:ncol(amrData)])
rownames(amrData_Integer) <- amrData$ID 
amrData_Integer <- amrData_Integer %>% tibble::rownames_to_column("ID")

## Remove genes that have count zero for all samples
GeneZeroValue <- colnames(amrData_Integer)[2:ncol(amrData_Integer)][apply(amrData_Integer[,2:ncol(amrData_Integer)], 2, sum) == 0]; length(GeneZeroValue)  # 211
amrData_Integer_NoZeroGene <- amrData_Integer %>% dplyr::select_if(!colnames(amrData_Integer) %in% GeneZeroValue ); dim(amrData_Integer_NoZeroGene)  # 185  583



## Make a new column, "Subtype" to see Rifaximin/Placebo and SA/ST
#sapply(stringr::str_split(amrData_Integer_NoZeroGene$ID, pattern="_"), function(X) {X[1] })  # "P"  or "R"
RifaPlacebo <- sapply(stringr::str_split(amrData_Integer_NoZeroGene$ID, pattern="_"), "[", 1)  # "P"  or "R"
SAST <-  sapply(stringr::str_split(amrData_Integer_NoZeroGene$ID, pattern="_"), "[", 4)      # "SA" or "ST

amr_Subtype <- amrData_Integer_NoZeroGene %>% dplyr::mutate(Subtype=paste0(RifaPlacebo,"_",SAST ) )
table(amr_Subtype$Subtype)
# P_SA P_ST R_SA R_ST 
# 47   42   48   48



########################################################### Stool subsets
## Extract only "P_ST" or "R_ST" samples of DayO1 only 
amr_Placebo_ST_Day01 <- amr_Subtype %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_ST_Day01) # 18  583  only P_ST Day01 

amr_Rifaximin_ST_Day01 <- amr_Subtype %>% dplyr::filter(grepl("R_ST", Subtype) & grepl("d01", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_ST_Day01) # 18  583  only R_ST Day01   

## Extract only "P_ST" or "R_ST" samples of Day30 only 
amr_Placebo_ST_Day30 <- amr_Subtype %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d30", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_ST_Day30) # 13  583  only P_ST Day30 

amr_Rifaximin_ST_Day30 <- amr_Subtype %>% dplyr::filter(grepl("R_ST", Subtype) & grepl("d30", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_ST_Day30) # 16  583  only R_ST Day30   

## Extract only "P_ST" or "R_ST" samples of Day90 only 
amr_Placebo_ST_Day90 <- amr_Subtype %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_ST_Day90) # 11  583  only P_ST Day90 

amr_Rifaximin_ST_Day90 <- amr_Subtype %>% dplyr::filter(grepl("R_ST", Subtype) & grepl("d90", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_ST_Day90) # 14  583  only R_ST Day90 


########################################################## Salva subsets
## Extract only "P_SA" or "R_SA" samples of DayO1 only 
amr_Placebo_SA_Day01 <- amr_Subtype %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_SA_Day01) # 19  583  only P_SA Day01 

amr_Rifaximin_SA_Day01 <- amr_Subtype %>% dplyr::filter(grepl("R_SA", Subtype) & grepl("d01", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_SA_Day01) # 18  583  only R_SA Day01   

## Extract only "P_SA" or "R_SA" samples of Day30 only 
amr_Placebo_SA_Day30 <- amr_Subtype %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d30", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_SA_Day30) # 16  583  only P_SA Day30 

amr_Rifaximin_SA_Day30 <- amr_Subtype %>% dplyr::filter(grepl("R_SA", Subtype) & grepl("d30", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_SA_Day30) # 16  583  only R_SA Day30   

## Extract only "P_SA" or "R_SA" samples of Day90 only 
amr_Placebo_SA_Day90 <- amr_Subtype %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_SA_Day90) # 12  583  only P_SA Day90 

amr_Rifaximin_SA_Day90 <- amr_Subtype %>% dplyr::filter(grepl("R_SA", Subtype) & grepl("d90", ID)) %>% dplyr::select(-Subtype); 
dim(amr_Rifaximin_SA_Day90) # 14  583  only R_SA Day90 


####################################################################################
## Step2. Run DESeq2 pairwise
####################################################################################

func_DESeq2 <- function(PlaceboData, RifaximinData ) {
      PlaceboRifaximin_Day <- rbind(PlaceboData, RifaximinData); #dim(PlaceboRifaximin_Day)  # 36 794
      DataType <- paste0(sapply(stringr::str_split(PlaceboRifaximin_Day$ID, pattern="_"), "[", 1),"_", 
                         sapply(stringr::str_split(PlaceboRifaximin_Day$ID, pattern="_"), "[", 3),"_",
                         sapply(stringr::str_split(PlaceboRifaximin_Day$ID, pattern="_"), "[", 4))
      table(DataType); # P_d01_ST: 18,  R_d01_ST: 18 
      
      ## DESeq2 annotation data
      DataTypeAnnotation <- PlaceboRifaximin_Day %>% dplyr::select(ID) %>% dplyr::mutate(DataType=DataType)
      
      ## DESeq2 input data
      DESeq2_Input <- PlaceboRifaximin_Day %>% tibble::column_to_rownames("ID") %>% t %>% data.frame %>% tibble::rownames_to_column("GeneSymbol")
      ## Add 1 as a pseudo data
      DESeq2_Input_Pseudo <-cbind( (DESeq2_Input[,"GeneSymbol"]),  (DESeq2_Input[,2:ncol(DESeq2_Input)] + 1) )
      colnames(DESeq2_Input_Pseudo)[1]<-"GeneSymbol"
      
      
      ## Run DESeq2
      DESeqData<-DESeqDataSetFromMatrix(countData=DESeq2_Input_Pseudo, colData=DataTypeAnnotation,design=~(DataType), tidy=TRUE)
      DESeqModel<-DESeq(DESeqData)  # DESeq model
      
      Result_DESeqModel<-results(DESeqModel, tidy=TRUE)
      colnames(Result_DESeqModel)[1]<-"GeneSymb"
      #change PeakName to PeakGeneStartEnd
      head(Result_DESeqModel)
      #print(Result_DESeqModel)
      
      table(is.na(Result_DESeqModel$padj))  #  TRUE 382
      Result_DESeqModel_NoNA <- Result_DESeqModel[!is.na(Result_DESeqModel$padj),]
      print(dim(Result_DESeqModel_NoNA)) # 411  7
      print(Result_DESeqModel_NoNA)
      return(Result_DESeqModel_NoNA)
}    
      
DESeq2Out_Stool_Day01 <- func_DESeq2(amr_Placebo_ST_Day01, amr_Rifaximin_ST_Day01)
dim(DESeq2Out_Stool_Day01); View(DESeq2Out_Stool_Day01)
DESeq2Out_Stool_Day30 <- func_DESeq2(amr_Placebo_ST_Day30, amr_Rifaximin_ST_Day30)
DESeq2Out_Stool_Day90 <- func_DESeq2(amr_Placebo_ST_Day90, amr_Rifaximin_ST_Day90)

DESeq2Out_Saliva_Day01 <- func_DESeq2(amr_Placebo_SA_Day01, amr_Rifaximin_SA_Day01)
dim(DESeq2Out_Saliva_Day01)
DESeq2Out_Saliva_Day30 <- func_DESeq2(amr_Placebo_SA_Day30, amr_Rifaximin_SA_Day30)
DESeq2Out_Saliva_Day90 <- func_DESeq2(amr_Placebo_SA_Day90, amr_Rifaximin_SA_Day90)

############################################################################################
## Step3A-1. Process the DESeq2 out datasets and rbind
############################################################################################

## Stool data: rbind for day01, day30, day90 WRST summary data. 
DESeq2Out_Stool_Day01_Proc <- DESeq2Out_Stool_Day01 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day01")
DESeq2Out_Stool_Day30_Proc <- DESeq2Out_Stool_Day30 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day30")
DESeq2Out_Stool_Day90_Proc <- DESeq2Out_Stool_Day90 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day90")

StoolData_DESeq2Summary_DayAll <- rbind(DESeq2Out_Stool_Day01_Proc, DESeq2Out_Stool_Day30_Proc, DESeq2Out_Stool_Day90_Proc)
dim(StoolData_DESeq2Summary_DayAll);table(StoolData_DESeq2Summary_DayAll$DataSource_Day)   # 663    4
##print(StoolData_DESeq2Summary_DayAll)

## Saliva data: rbind for day01, day30, day90 WRST summary data. 
DESeq2Out_Saliva_Day01_Proc <- DESeq2Out_Saliva_Day01 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day01")
DESeq2Out_Saliva_Day30_Proc <- DESeq2Out_Saliva_Day30 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day30")
DESeq2Out_Saliva_Day90_Proc <- DESeq2Out_Saliva_Day90 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day90")

SalivaData_DESeq2Summary_DayAll <- rbind(DESeq2Out_Saliva_Day01_Proc, DESeq2Out_Saliva_Day30_Proc, DESeq2Out_Saliva_Day90_Proc)
dim(SalivaData_DESeq2Summary_DayAll);table(SalivaData_DESeq2Summary_DayAll$DataSource_Day)   # 1746  4
View(SalivaData_DESeq2Summary_DayAll)

############################################################################################
## Step3A-2. Filter in genes of significant p-val < 0.001. And make dotplot
############################################################################################
StoolData_DESeq2Summary_Sort <- StoolData_DESeq2Summary_DayAll[with(StoolData_DESeq2Summary_DayAll, order(pvalue ,decreasing=FALSE)),] %>%
                dplyr::filter(pvalue <= 1e-03); length(StoolData_DESeq2Summary_Sort$GeneSymb)  # 63

StoolData_DESeq2Summary_DayAll_Sig <- StoolData_DESeq2Summary_DayAll %>% dplyr::filter(StoolData_DESeq2Summary_DayAll$GeneSymb %in% StoolData_DESeq2Summary_Sort$GeneSymb) 
length(StoolData_DESeq2Summary_DayAll_Sig$GeneSymb); dim(StoolData_DESeq2Summary_DayAll_Sig) # 142 4     
unique(StoolData_DESeq2Summary_DayAll_Sig$GeneSymb)

SalivaData_DESeq2Summary_Sort <- SalivaData_DESeq2Summary_DayAll[with(SalivaData_DESeq2Summary_DayAll, order(pvalue ,decreasing=FALSE)),] %>%
                dplyr::filter(pvalue <= 1e-02); length(SalivaData_DESeq2Summary_Sort$GeneSymb) #13

SalivaData_DESeq2Summary_DayAll_Sig <- SalivaData_DESeq2Summary_DayAll %>% dplyr::filter(SalivaData_DESeq2Summary_DayAll$GeneSymb %in% SalivaData_DESeq2Summary_Sort$GeneSymb)
length(SalivaData_DESeq2Summary_DayAll_Sig$GeneSymb); dim(SalivaData_DESeq2Summary_DayAll_Sig) # 33 4     
unique(SalivaData_DESeq2Summary_DayAll_Sig$GeneSymb)

MyDotplotInput_Comb_ST <- StoolData_DESeq2Summary_DayAll_Sig ;dim(MyDotplotInput_Comb_ST)  # 142  4
MyDotplotInput_Comb_SA <- SalivaData_DESeq2Summary_DayAll_Sig ;dim(MyDotplotInput_Comb_SA)  # 33  4
#View(MyDotplotInput_Comb)  
#sort(unique(MyDotplotInput_Comb$DataSource_Day), decreasing=FALSE)


SampleLevel_SA <- c("Saliva_Day01", "Saliva_Day30", "Saliva_Day90") 
SampleLevel_ST <- c("Stool_Day01", "Stool_Day30","Stool_Day90")
MyDotplotInput_Comb_ST$DataSource_Day <- factor(MyDotplotInput_Comb_ST$DataSource_Day, levels=SampleLevel_ST)
MyDotplotInput_Comb_SA$DataSource_Day <- factor(MyDotplotInput_Comb_SA$DataSource_Day, levels=SampleLevel_SA)
View(MyDotplotInput_Comb_ST)
View(MyDotplotInput_Comb_SA)

## Make dotplot

#filename = 'plot.pdf';width = 8;height = 10;means_separator = '\t';pvalues_separator = '\t';output_extension = '.pdf'
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

MyDotlot_ST <- ggplot(MyDotplotInput_Comb_ST, aes(x=GeneSymb, y=DataSource_Day)) +
  geom_point(aes(size=-log10(pvalue),color=log2FoldChange)) +
  scale_color_gradientn('Log2 Fold Change', colors=my_palette) +
  labs(title=paste0("Log2 Fold Change between Placebo and Rifaximin groups with DESeq2, Stool"))+
  theme_bw() +   # coord_flip() +
  theme(plot.title=element_text(size=17),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
MyDotlot_ST

ggsave(MyDotlot, filename=paste0("OutDotplot_DESeq2FCpval_AMR_Placebo_vs_Rifaximin_W15H10.pdf"), width=15, height=10, limitsize=F)


###########   Saliva 

my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

MyDotlot_SA <- ggplot(MyDotplotInput_Comb_SA, aes(x=GeneSymb, y=DataSource_Day)) +
  geom_point(aes(size=-log10(pvalue),color=log2FoldChange)) +
  scale_color_gradientn('Log2 Fold Change', colors=my_palette) +
  labs(title=paste0("Log2 Fold Change between Placebo and Rifaximin groups with DESeq2, Saliva"))+
  theme_bw() +   # coord_flip() +
  theme(plot.title=element_text(size=17),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

MyDotlot_SA


############################################################################################
## Step3B. Make volcano plot and annotate significant genes. 
############################################################################################
#####################################  ST  01  #######################################################

DESeq2Out_Stool_Day01_Proc <- DESeq2Out_Stool_Day01 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day01")
StoolData01_DESeq2Summary_Sort <- DESeq2Out_Stool_Day01_Proc[with(DESeq2Out_Stool_Day01_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(StoolData01_DESeq2Summary_Sort$GeneSymb)  # 65

StoolData01_DESeq2Summary_Sig <- DESeq2Out_Stool_Day01_Proc %>% dplyr::filter(DESeq2Out_Stool_Day01_Proc$GeneSymb %in% StoolData01_DESeq2Summary_Sort$GeneSymb) 
length(StoolData01_DESeq2Summary_Sig$GeneSymb); dim(StoolData01_DESeq2Summary_Sig) # 65  4
View(StoolData01_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
StoolData01_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(StoolData01_DESeq2Summary_Sig$log2FoldChange) > 2.0 & StoolData01_DESeq2Summary_Sig$pvalue < 5e-02, StoolData01_DESeq2Summary_Sig$GeneSymb, "")
StoolData01_DESeq2Summary_Sig$Activity=ifelse(StoolData01_DESeq2Summary_Sig$pvalue < 5e-02 & abs(StoolData01_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                              ifelse(StoolData01_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

StoolData01_DESeq2Summary_Sig$YAxis <- -log(StoolData01_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(StoolData01_DESeq2Summary_Sig$log2FoldChange)  # -2.35404
max(StoolData01_DESeq2Summary_Sig$log2FoldChange)  # 3.712837

## To see the y-axis range
min(StoolData01_DESeq2Summary_Sig$YAxis)  # 1.301275
max(StoolData01_DESeq2Summary_Sig$YAxis)  # 5.703792

View(StoolData01_DESeq2Summary_Sig)

p_ST_01 <- ggplot(data=StoolData01_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=3.5, show.legend=F) +
  scale_color_manual(values=c("dark blue", "black","red"))+
  xlim(c(-5, 5)) + ylim(c(0, 8.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=18) vs Rifaximin (n=18) on Day01, Stool")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_ST_01 <-p_ST_01 + geom_text_repel(box.padding=0.5,  max.overlaps = Inf)
MyPlot_ST_01
unique(StoolData01_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaST_vs_PalaceboST_Day01.pdf", width=8, height=8)
#################################################################################################################################
##############################################  ST   30  ##################################################
DESeq2Out_Stool_Day30_Proc <- DESeq2Out_Stool_Day30 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day30")
dim(DESeq2Out_Stool_Day30_Proc)  # 221   4
StoolData30_DESeq2Summary_Sort <- DESeq2Out_Stool_Day30_Proc[with(DESeq2Out_Stool_Day30_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(StoolData30_DESeq2Summary_Sort$GeneSymb)  # 87

StoolData30_DESeq2Summary_Sig <- DESeq2Out_Stool_Day30_Proc %>% dplyr::filter(DESeq2Out_Stool_Day30_Proc$GeneSymb %in% StoolData30_DESeq2Summary_Sort$GeneSymb) 
length(StoolData30_DESeq2Summary_Sig$GeneSymb); dim(StoolData30_DESeq2Summary_Sig) # 87  4
View(StoolData30_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
StoolData30_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(StoolData30_DESeq2Summary_Sig$log2FoldChange) > 2.0 & StoolData30_DESeq2Summary_Sig$pvalue < 5e-02, StoolData30_DESeq2Summary_Sig$GeneSymb, "")
StoolData30_DESeq2Summary_Sig$Activity=ifelse(StoolData30_DESeq2Summary_Sig$pvalue < 5e-02 & abs(StoolData30_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                              ifelse(StoolData30_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

StoolData30_DESeq2Summary_Sig$YAxis <- -log(StoolData30_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(StoolData30_DESeq2Summary_Sig$log2FoldChange)  # -2.163346
max(StoolData30_DESeq2Summary_Sig$log2FoldChange)  # 5.14372

## To see the y-axis range
min(StoolData30_DESeq2Summary_Sig$YAxis)  # 1.303896
max(StoolData30_DESeq2Summary_Sig$YAxis)  # 6.803912

View(StoolData30_DESeq2Summary_Sig)

p_ST_30 <- ggplot(data=StoolData30_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=2.2, show.legend=F) +
  scale_color_manual(values=c("dark blue", "black","red"))+
  xlim(c(-9, 9)) + ylim(c(0, 9.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=13) vs Rifaximin (n=16) on Day30, Stool")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_ST_30 <-p_ST_30 + geom_text_repel(box.padding=0.5, max.overlaps = Inf, size=2.6, force=4.8, min.segment.length = 0.6, force_pull=0.1, max.time=2.5)
MyPlot_ST_30
unique(StoolData_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaST_vs_PalaceboST_Day30.pdf", width=8, height=8)

############################################  ST  90  ###########################################################

DESeq2Out_Stool_Day90_Proc <- DESeq2Out_Stool_Day90 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Stool_Day90")
dim(DESeq2Out_Stool_Day90_Proc)  #  176   4
StoolData90_DESeq2Summary_Sort <- DESeq2Out_Stool_Day90_Proc[with(DESeq2Out_Stool_Day90_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(StoolData90_DESeq2Summary_Sort$GeneSymb)  # 76

StoolData90_DESeq2Summary_Sig <- DESeq2Out_Stool_Day90_Proc %>% dplyr::filter(DESeq2Out_Stool_Day90_Proc$GeneSymb %in% StoolData90_DESeq2Summary_Sort$GeneSymb) 
length(StoolData90_DESeq2Summary_Sig$GeneSymb); dim(StoolData90_DESeq2Summary_Sig) # 76  4
View(StoolData90_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
StoolData90_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(StoolData90_DESeq2Summary_Sig$log2FoldChange) > 2.0 & StoolData90_DESeq2Summary_Sig$pvalue < 5e-02, StoolData90_DESeq2Summary_Sig$GeneSymb, "")
StoolData90_DESeq2Summary_Sig$Activity=ifelse(StoolData90_DESeq2Summary_Sig$pvalue < 5e-02 & abs(StoolData90_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                              ifelse(StoolData90_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

StoolData90_DESeq2Summary_Sig$YAxis <- -log(StoolData90_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(StoolData90_DESeq2Summary_Sig$log2FoldChange)  # -3.117044
max(StoolData90_DESeq2Summary_Sig$log2FoldChange)  # 4.379913

## To see the y-axis range
min(StoolData90_DESeq2Summary_Sig$YAxis)  # 1.363148
max(StoolData90_DESeq2Summary_Sig$YAxis)  # 5.223238

View(StoolData90_DESeq2Summary_Sig)
## I tired 'xlim(c(-0.05, 0.05))' but the plot doesn't look pretty. 
p_ST_90 <- ggplot(data=StoolData90_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=2.5, show.legend=F) +
  scale_color_manual(values=c("dark blue","black","red"))+
  xlim(c(-8, 8)) + ylim(c(0, 7.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=11) vs Rifaximin (n=14) on Day90, Stool")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_ST_90 <-p_ST_90 + geom_text_repel(box.padding=0.5, max.overlaps = Inf, size=2.6, force=4.0, min.segment.length = 0.7, force_pull=0.3, max.time=5.5)
MyPlot_ST_90
unique(StoolData90_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaST_vs_PalaceboST_Day90.pdf", width=8, height=8)

##############################################  SA  01 ###################################################################

DESeq2Out_Saliva_Day01_Proc <- DESeq2Out_Saliva_Day01 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day01")
SalivaData01_DESeq2Summary_Sort <- DESeq2Out_Saliva_Day01_Proc[with(DESeq2Out_Saliva_Day01_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(SalivaData01_DESeq2Summary_Sort$GeneSymb)  # 22

SalivaData01_DESeq2Summary_Sig <- DESeq2Out_Saliva_Day01_Proc %>% dplyr::filter(DESeq2Out_Saliva_Day01_Proc$GeneSymb %in% SalivaData01_DESeq2Summary_Sort$GeneSymb) 
length(SalivaData01_DESeq2Summary_Sig$GeneSymb); dim(SalivaData01_DESeq2Summary_Sig) # 22  4
View(SalivaData01_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
SalivaData01_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(SalivaData01_DESeq2Summary_Sig$log2FoldChange) > 2.0 & SalivaData01_DESeq2Summary_Sig$pvalue < 5e-02, SalivaData01_DESeq2Summary_Sig$GeneSymb, "")
SalivaData01_DESeq2Summary_Sig$Activity=ifelse(SalivaData01_DESeq2Summary_Sig$pvalue < 5e-02 & abs(SalivaData01_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                               ifelse(SalivaData01_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

SalivaData01_DESeq2Summary_Sig$YAxis <- -log(SalivaData01_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(SalivaData01_DESeq2Summary_Sig$log2FoldChange)  # -2.266025
max(SalivaData01_DESeq2Summary_Sig$log2FoldChange)  # 1.810847

## To see the y-axis range
min(SalivaData01_DESeq2Summary_Sig$YAxis)  # 1.317024
max(SalivaData01_DESeq2Summary_Sig$YAxis)  # 3.805208

View(SalivaData01_DESeq2Summary_Sig)
## I tired 'xlim(c(-0.05, 0.05))' but the plot doesn't look pretty. 
p_SA_01 <- ggplot(data=SalivaData01_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=2.5, show.legend=F) +
  scale_color_manual(values=c("dark blue", "black","red"))+
  xlim(c(-5, 5)) + ylim(c(0, 5.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=19) vs Rifaximin (n=18) on Day01, Saliva")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_SA_01 <-p_SA_01 + geom_text_repel(box.padding=0.5, max.overlaps = Inf)
MyPlot_SA_01
unique(StoolData01_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaSA_vs_PalaceboSA_Day01.pdf", width=8, height=8)
#################################################################################################################################
##############################################   SA  30   ########################################################################

DESeq2Out_Saliva_Day30_Proc <- DESeq2Out_Saliva_Day30 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day30")
SalivaData30_DESeq2Summary_Sort <- DESeq2Out_Saliva_Day30_Proc[with(DESeq2Out_Saliva_Day30_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(SalivaData30_DESeq2Summary_Sort$GeneSymb)  # 23

SalivaData30_DESeq2Summary_Sig <- DESeq2Out_Saliva_Day30_Proc %>% dplyr::filter(DESeq2Out_Saliva_Day30_Proc$GeneSymb %in% SalivaData30_DESeq2Summary_Sort$GeneSymb) 
length(SalivaData30_DESeq2Summary_Sig$GeneSymb); dim(SalivaData30_DESeq2Summary_Sig) # 23  4
View(SalivaData30_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
SalivaData30_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(SalivaData30_DESeq2Summary_Sig$log2FoldChange) > 2.0 & SalivaData30_DESeq2Summary_Sig$pvalue < 5e-02, SalivaData30_DESeq2Summary_Sig$GeneSymb, "")
SalivaData30_DESeq2Summary_Sig$Activity=ifelse(SalivaData30_DESeq2Summary_Sig$pvalue < 5e-02 & abs(SalivaData30_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                               ifelse(SalivaData30_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

SalivaData30_DESeq2Summary_Sig$YAxis <- -log(SalivaData30_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(SalivaData30_DESeq2Summary_Sig$log2FoldChange)  # -2.046292
max(SalivaData30_DESeq2Summary_Sig$log2FoldChange)  # 2.152002

## To see the y-axis range
min(SalivaData30_DESeq2Summary_Sig$YAxis)  # 1.30166
max(SalivaData30_DESeq2Summary_Sig$YAxis)  # 4.168613

View(SalivaData30_DESeq2Summary_Sig)
## I tired 'xlim(c(-0.05, 0.05))' but the plot doesn't look pretty. 
p_SA_30 <- ggplot(data=SalivaData30_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=3.5, show.legend=F) +
  scale_color_manual(values=c("dark blue", "black","red"))+
  xlim(c(-5, 5)) + ylim(c(0, 5.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=16) vs Rifaximin (n=16) on Day30, Saliva")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_SA_30 <-p_SA_30 + geom_text_repel(box.padding=0.5,  max.overlaps = Inf)
MyPlot_SA_30
unique(StoolData01_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaSA_vs_PalaceboSA_Day30.pdf", width=8, height=8)

##############################################  SA  90  #########################################################################

DESeq2Out_Saliva_Day90_Proc <- DESeq2Out_Saliva_Day90 %>% dplyr::select(-c(baseMean,lfcSE,stat,padj)) %>% dplyr::mutate(DataSource_Day="Saliva_Day90")
SalivaData90_DESeq2Summary_Sort <- DESeq2Out_Saliva_Day90_Proc[with(DESeq2Out_Saliva_Day90_Proc, order(pvalue ,decreasing=FALSE)),] %>%
  dplyr::filter(pvalue <= 0.9); length(SalivaData90_DESeq2Summary_Sort$GeneSymb)  # 10

SalivaData90_DESeq2Summary_Sig <- DESeq2Out_Saliva_Day90_Proc %>% dplyr::filter(DESeq2Out_Saliva_Day90_Proc$GeneSymb %in% SalivaData90_DESeq2Summary_Sort$GeneSymb) 
length(SalivaData90_DESeq2Summary_Sig$GeneSymb); dim(SalivaData90_DESeq2Summary_Sig) # 10  4
View(SalivaData90_DESeq2Summary_Sig)

## ggrepel annotation by p-value 0.01
SalivaData90_DESeq2Summary_Sig$GeneAnnotation <- ifelse(abs(SalivaData90_DESeq2Summary_Sig$log2FoldChange) > 2.0 & SalivaData90_DESeq2Summary_Sig$pvalue < 5e-02, SalivaData90_DESeq2Summary_Sig$GeneSymb, "")
SalivaData90_DESeq2Summary_Sig$Activity=ifelse(SalivaData90_DESeq2Summary_Sig$pvalue < 5e-02 & abs(SalivaData90_DESeq2Summary_Sig$log2FoldChange) >= 2.0, 
                                               ifelse(SalivaData90_DESeq2Summary_Sig$log2FoldChange > 2.0 ,'Up','Down'), 'Stable')

SalivaData90_DESeq2Summary_Sig$YAxis <- -log(SalivaData90_DESeq2Summary_Sig$pvalue,10)

## To see the x-axis range
min(SalivaData90_DESeq2Summary_Sig$log2FoldChange)  # -1.452689
max(SalivaData90_DESeq2Summary_Sig$log2FoldChange)  # 1.582645

## To see the y-axis range
min(SalivaData90_DESeq2Summary_Sig$YAxis)  # 1.334206
max(SalivaData90_DESeq2Summary_Sig$YAxis)  # 1.782731

View(SalivaData90_DESeq2Summary_Sig)

p_SA_90 <- ggplot(data=SalivaData90_DESeq2Summary_Sig,aes(x=log2FoldChange, y=YAxis, colour=Activity,label=GeneAnnotation)) +
  geom_point(alpha=0.4, size=3.5, show.legend=F) +
  scale_color_manual(values=c("dark blue", "black","red"))+
  xlim(c(-4, 4)) + ylim(c(0, 4.0)) +
  geom_vline(xintercept=c(-2.0, 2.0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept=-log(5e-02, 10), lty=4, col="black", lwd=0.8) +
  labs(x="log2(fold change)", y="-log10(p-value)",title="Placebo (n=12) vs Rifaximin (n=14) on Day90, Saliva")  + theme_bw() +
  theme(legend.position="none", 
        axis.text.x=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=20,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.x=element_text(colour="black",size=22,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=22,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.2),
        axis.line.x=element_line(color="black", size =1.2),
        axis.ticks=element_line(colour="black",size=1.7),
        axis.ticks.length=unit(.22, "cm"), text=element_text(size=20))

MyPlot_SA_90 <-p_SA_90 + geom_text_repel(box.padding=0.5,  max.overlaps = Inf)
MyPlot_SA_90
unique(StoolData01_DESeq2Summary_Sig$GeneAnnotation)
ggsave(MyPlot, file="VolcanoPlot_DESeq2_RifaSA_vs_PalaceboSA_Day90.pdf", width=8, height=8)

#################################################################################################################################
#################################################################################################################################
## Step4. Make Fraction barplot - the samples should be ordered by Day
## Plot R code: https://urldefense.com/v3/__https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2__;!!NHLzug!Oej5-ovJpXHABep3nXT1D7tNoh5smwRQ12wY10BHcuO85H8699IOiIeDqx2mSKi9r3xbh1UUZuk04v4pEIoDalqa$ 
#################################################################################################################################
######################################################  ST Placebo  #############################################################
View(amr_Subtype)
amr_Placebo_ST_Day <- amr_Subtype %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_ST_Day)   # 42  583 by d01, d30, d90
View(amr_Placebo_ST_Day)
colnames(amr_Placebo_ST_Day)
## Make tidyr data format
amr_Placebo_ST_Day_tidyr <- amr_Placebo_ST_Day %>% tidyr::gather(AMRGene, Abundance, AAC2__Ib:vph);
head(amr_Placebo_ST_Day_tidyr); dim(amr_Placebo_ST_Day_tidyr) # 24444 3
View(amr_Placebo_ST_Day_tidyr)
## Remove rows of Abudance==0
amr_Placebo_ST_Day_tidyr_NoZero <- amr_Placebo_ST_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Placebo_ST_Day_tidyr_NoZero) # 4480 3
View(amr_Placebo_ST_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency. 
GeneFrequency <- table(amr_Placebo_ST_Day_tidyr_NoZero$AMRGene)
length(GeneFrequency)  # there are all 412 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 42  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>35]; length(GeneFrequency_High)  # 15 genes left.
View(GeneFrequency_High)
## Filter in just 15 AMRGenes of high frequency.
amr_Placebo_ST_Day_tidyr_HigFreq <- amr_Placebo_ST_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
## Check how many kinds of genes are contained
length(table(amr_Placebo_ST_Day_tidyr_HigFreq$AMRGene))  # 15, Ok this small number is good to make fraction barplot.   
View(amr_Placebo_ST_Day_tidyr_HigFreq)

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
ID_D30 <- unique(grep("d30", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 13 samples
ID_D90 <- unique(grep("d90", amr_Placebo_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 11 samples

amr_Placebo_ST_Day_tidyr_HigFreq$ID <- factor(amr_Placebo_ST_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))

MyBarplot_P_ST <- ggplot(amr_Placebo_ST_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot_P_ST

## I want to change the bar colors manually and adjust the axis text.

MyBarplot_9CellType_P_ST <- MyBarplot_P_ST + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Placebo group, Stool"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_P_ST

## Let's save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_P_ST, file="FractionBarplot_DES_Placebo_ST.pdf", width=14, height=6)


############################################## ST Rifaximin ##########################################################
View(amr_Subtype)
amr_Rifaximin_ST_Day <- amr_Subtype %>% dplyr::filter(grepl("R_ST", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Rifaximin_ST_Day)   # 48  583 by d01, d30, d90
View(amr_Rifaximin_ST_Day)
colnames(amr_Rifaximin_ST_Day)
## Make tidyr data format
amr_Rifaximin_ST_Day_tidyr <- amr_Rifaximin_ST_Day %>% tidyr::gather(AMRGene, Abundance, AAC2__Ib:vph);
head(amr_Rifaximin_ST_Day_tidyr); dim(amr_Rifaximin_ST_Day_tidyr) # 27936  3
View(amr_Rifaximin_ST_Day_tidyr)

## Remove rows of Abudance==0
amr_Rifaximin_ST_Day_tidyr_NoZero <- amr_Rifaximin_ST_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Rifaximin_ST_Day_tidyr_NoZero) # 4722    3
View(amr_Rifaximin_ST_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency.
GeneFrequency <- table(amr_Rifaximin_ST_Day_tidyr_NoZero$AMRGene)
View(GeneFrequency)
length(GeneFrequency)  # there are all 419 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 48  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. 
View(GeneFrequency_High)
## Filter in just 15 AMRGenes of high frequency.
amr_Rifaximin_ST_Day_tidyr_HigFreq <- amr_Rifaximin_ST_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
View(amr_Rifaximin_ST_Day_tidyr_HigFreq)
## Check how many kinds of genes are contained
length(table(amr_Rifaximin_ST_Day_tidyr_HigFreq$AMRGene))  # 14,  this small number is good to make fraction barplot.   

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
ID_D30 <- unique(grep("d30", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 16 samples
ID_D90 <- unique(grep("d90", amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 14 samples

amr_Rifaximin_ST_Day_tidyr_HigFreq$ID <- factor(amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))


MyBarplot_R_ST <- ggplot(amr_Rifaximin_ST_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot_R_ST


MyBarplot_9CellType_R_ST <- MyBarplot_R_ST + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Rifaximin group, Stool"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_R_ST

## Let's save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_R_ST, file="FractionBarplot_DES_Rifaximin_ST.pdf", width=14, height=6)

##############################################  SA Placebo ####################################################
View(amr_Subtype)
amr_Placebo_SA_Day <- amr_Subtype %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Placebo_SA_Day)   # 47  583 by d01, d30, d90
View(amr_Placebo_SA_Day)
colnames(amr_Placebo_SA_Day)
## Make tidyr data format
amr_Placebo_SA_Day_tidyr <- amr_Placebo_SA_Day %>% tidyr::gather(AMRGene, Abundance, AAC2__Ib:vph);
head(amr_Placebo_SA_Day_tidyr); dim(amr_Placebo_SA_Day_tidyr) # 27354     3
View(amr_Placebo_SA_Day_tidyr)
## Remove rows of Abudance==0
amr_Placebo_SA_Day_tidyr_NoZero <- amr_Placebo_SA_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Placebo_SA_Day_tidyr_NoZero) # 2979    3
View(amr_Placebo_SA_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency. 
GeneFrequency <- table(amr_Placebo_SA_Day_tidyr_NoZero$AMRGene)
length(GeneFrequency)  # there are all 319 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 47  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. 
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
## Let's see MyBarplot
MyBarplot_P_SA



MyBarplot_9CellType_P_SA <- MyBarplot_P_SA + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Placebo group, Saliva"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_P_SA

## Let's save the fraction barplot to PDF file. 
ggsave(MyBarplot_9CellType_P_SA, file="FractionBarplot_DES_Placebo_SA.pdf", width=14, height=6)

############################################## SA Rifaximin    ##########################################################
View(amr_Subtype)
amr_Rifaximin_SA_Day <- amr_Subtype %>% dplyr::filter(grepl("R_SA", Subtype)  & grepl("d01|d30|d90", ID)) %>% dplyr::select(-Subtype)
dim(amr_Rifaximin_SA_Day)   # 48  583 by d01, d30, d90
View(amr_Rifaximin_SA_Day)
colnames(amr_Rifaximin_SA_Day)
## Make tidyr data format
amr_Rifaximin_SA_Day_tidyr <- amr_Rifaximin_SA_Day %>% tidyr::gather(AMRGene, Abundance, AAC2__Ib:vph);
head(amr_Rifaximin_SA_Day_tidyr); dim(amr_Rifaximin_ST_Day_tidyr) # 27936  3
View(amr_Rifaximin_SA_Day_tidyr)

## Remove rows of Abudance==0
amr_Rifaximin_SA_Day_tidyr_NoZero <- amr_Rifaximin_SA_Day_tidyr %>% dplyr::filter(Abundance!=0);
dim(amr_Rifaximin_SA_Day_tidyr_NoZero) # 2926  3
View(amr_Rifaximin_SA_Day_tidyr_NoZero)

## Count the frequency of the AMR genes. Some genes appear with low frequency. 
GeneFrequency <- table(amr_Rifaximin_SA_Day_tidyr_NoZero$AMRGene)
View(GeneFrequency)
length(GeneFrequency)  # there are all 302 AMRGenes.
max(GeneFrequency); min(GeneFrequency)   # 47  1  
GeneFrequency_High <- GeneFrequency[GeneFrequency>37]; length(GeneFrequency_High)  # 14 genes left. 
View(GeneFrequency_High)
## Filter in just 15 AMRGenes of high frequency.
amr_Rifaximin_SA_Day_tidyr_HigFreq <- amr_Rifaximin_SA_Day_tidyr_NoZero %>% dplyr::filter(AMRGene %in% names(GeneFrequency_High))
View(amr_Rifaximin_SA_Day_tidyr_HigFreq)
## Check how many kinds of genes are contained
length(table(amr_Rifaximin_SA_Day_tidyr_HigFreq$AMRGene))  # 14, this small number is good to make fraction barplot.   

## Sort the data by "ID" column
ID_D01 <- unique(grep("d01", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D01)  ## 18 samples
ID_D30 <- unique(grep("d30", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D30)  ## 16 samples
ID_D90 <- unique(grep("d90", amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, value=TRUE)); length(ID_D90)  ## 14 samples

amr_Rifaximin_SA_Day_tidyr_HigFreq$ID <- factor(amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, levels=c(ID_D01,ID_D30,ID_D90))
View(amr_Rifaximin_SA_Day_tidyr_HigFreq)
MyBarplot_R_SA <- ggplot(amr_Rifaximin_SA_Day_tidyr_HigFreq, aes(fill=AMRGene, y=Abundance, x=ID))+ geom_bar(position="fill", stat="identity") 
## Let's see MyBarplot
MyBarplot_R_SA


MyBarplot_9CellType_R_SA <- MyBarplot_R_SA + 
  scale_fill_manual(values=c("#5BCDE4","yellow","grey","#08954F", "#C78F0B","blue","brown", "#CC99FF","#CCCC00","#336600",
                             "orange","#FF3399", "#FFCCCC","#E0E0E0","#FFE5CC")) + theme_bw() + 
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) + 
  labs(title=paste0("Relative abundance of high-frequency AMR genes in Rifaximin group, Saliva"), x="Sample ID",y="AMR gene abaudance",face="bold")+
  theme(axis.text.x=element_text(size=10,colour="black",angle=90,hjust=0.95,vjust=0.2), axis.text.y=element_text(size=10,colour="black", angle=00,hjust=0.95,vjust=0.2), 
        axis.title.y=element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold")) + scale_y_continuous(expand = c(0,0))

## See the barplot of different bar colors
MyBarplot_9CellType_R_SA

ggsave(MyBarplot_9CellType_R_ST, file="FractionBarplot_DES_Rifaximin_SA.pdf", width=14, height=6)

############################################################################################################################
## Step5. Prepare the boxplot input data 
##############################################  ST   Placebo   #############################################################
## Print "amr_Placebo_ST_Day_tidyr_HigFreq" again and check the data
head(amr_Placebo_ST_Day_tidyr_HigFreq)    # This data has "ID", "AMRGene", "Abundance" columns. 

## From the ID column,extract onl y "d01", "d30", "d90", and make "Day013090" column.  
## And then,  make "Day013090_AMRGene" column that has "dOO_BahA"   Please see the line 196 below
Day013090_Data_P <- sapply(stringr::str_split(amr_Placebo_ST_Day_tidyr_HigFreq$ID, pattern="_"), "[", 3); 
head(Day013090_Data_P)  # "d01" "d30" "d01" "d30" "d90" "d01"

amr_Placebo_ST_Day_Processed <- amr_Placebo_ST_Day_tidyr_HigFreq %>% dplyr::mutate(Day013090=Day013090_Data_P,
                                                                                   Day013090_AMRGene=paste0(Day013090_Data_P,"_", amr_Placebo_ST_Day_tidyr_HigFreq$AMRGene))
head(amr_Placebo_ST_Day_Processed)

###########################################  ST  Rifaximin  ###############################################################
## Print "amr_Placebo_ST_Day_tidyr_HigFreq" again and check the data
head(amr_Rifaximin_ST_Day_tidyr_HigFreq)    # This data has "ID", "AMRGene", "Abundance" columns. 

## From the ID column, extract only "d01", "d30", "d90", and make "Day013090" column. 
## And then,  make "Day013090_AMRGene" column that has "dOO_BahA"   
Day013090_Data_R <- sapply(stringr::str_split(amr_Rifaximin_ST_Day_tidyr_HigFreq$ID, pattern="_"), "[", 3); 
head(Day013090_Data_R)  # "d01" "d30" "d01" "d30" "d90" "d01"

amr_Rifaximin_ST_Day_Processed <- amr_Rifaximin_ST_Day_tidyr_HigFreq %>% dplyr::mutate(Day013090=Day013090_Data_R,
                                                                                       Day013090_AMRGene=paste0(Day013090_Data_R,"_", amr_Rifaximin_ST_Day_tidyr_HigFreq$AMRGene))
head(amr_Rifaximin_ST_Day_Processed)
View(amr_Rifaximin_ST_Day_Processed)

##############################################  SA   Placebo   #############################################################
## Print "amr_Placebo_ST_Day_tidyr_HigFreq" again and check the data
head(amr_Placebo_SA_Day_tidyr_HigFreq)    # This data has "ID", "AMRGene", "Abundance" columns. 

## From the ID column, extract only "d01", "d30", "d90", and make "Day013090" column.  
## And then, I should make "Day013090_AMRGene" column that has "dOO_BahA"   
Day013090_Data_P_SA <- sapply(stringr::str_split(amr_Placebo_SA_Day_tidyr_HigFreq$ID, pattern="_"), "[", 3); 
head(Day013090_Data_P_SA)  # "d01" "d30" "d01" "d30" "d90" "d01"

amr_Placebo_SA_Day_Processed <- amr_Placebo_SA_Day_tidyr_HigFreq %>% dplyr::mutate(Day013090=Day013090_Data_P_SA,
                                                                                   Day013090_AMRGene=paste0(Day013090_Data_P_SA,"_", amr_Placebo_SA_Day_tidyr_HigFreq$AMRGene))
head(amr_Placebo_SA_Day_Processed)

########################################### SA  Rifaximin  ###############################################################
## Print "amr_Placebo_ST_Day_tidyr_HigFreq" again and check the data
head(amr_Rifaximin_SA_Day_tidyr_HigFreq)    # This data has "ID", "AMRGene", "Abundance" columns. 

## From the ID column, only "d01", "d30", "d90", and make "Day013090" column.  
## And then, make "Day013090_AMRGene" column that has "dOO_BahA"   Please see the line 196 below
Day013090_Data_R_SA <- sapply(stringr::str_split(amr_Rifaximin_SA_Day_tidyr_HigFreq$ID, pattern="_"), "[", 3); 
head(Day013090_Data_R_SA)  # "d01" "d30" "d01" "d30" "d90" "d01"

amr_Rifaximin_SA_Day_Processed <- amr_Rifaximin_SA_Day_tidyr_HigFreq %>% dplyr::mutate(Day013090=Day013090_Data_R_SA,
                                                                                       Day013090_AMRGene=paste0(Day013090_Data_R_SA,"_", amr_Rifaximin_SA_Day_tidyr_HigFreq$AMRGene))
head(amr_Rifaximin_SA_Day_Processed)
View(amr_Rifaximin_SA_Day_Processed)


##############################################################################################################################
################ Step6. Make Boxplot with the gene abundance fraction  - Day01 vs Day30 vs Day90  ############################
## Lecture for boxplot:   https://urldefense.com/v3/__https://stackoverflow.com/questions/22808219/ggplot2-how-do-you-control-point-color-with-geom-jitter-and-geom-point__;!!NHLzug!Oej5-ovJpXHABep3nXT1D7tNoh5smwRQ12wY10BHcuO85H8699IOiIeDqx2mSKi9r3xbh1UUZuk04v4pENuXt5GD$ 
##############################################################################################################################
################################################   Placebo ST   ##############################################################
## Check the min and max of 'Abudnance' data because need to set the min and max in the y-axis of the boxplot. 
MinYaxis<-floor(min(as.numeric(amr_Placebo_ST_Day_Processed[,"Abundance"])));  print(MinYaxis) # 1
MaxYaxis<-round(max(as.numeric(amr_Placebo_ST_Day_Processed[,"Abundance"])), digits=1); print(MaxYaxis)  #17748


## Let's change the X-axis (Day013090_AMRGene) order
AMRGene_Unique <- unique(amr_Placebo_ST_Day_Processed$AMRGene) %>% sort #alphabetical order
print(AMRGene_Unique) #15

Day_AMRGene_Comb <- c()
for(EachAMRGene in AMRGene_Unique) { 
  # EachAMRGene <- AMRGene_Unique[1]
  amr_Placebo_ST_Day_Subset <- amr_Placebo_ST_Day_Processed %>% dplyr::filter(AMRGene == EachAMRGene); print(table(amr_Placebo_ST_Day_Subset$AMRGene))
  Day_AMRGene_Sort <- sort(amr_Placebo_ST_Day_Subset$Day013090_AMRGene) #day 01->30->90 within a specific gene
  print(Day_AMRGene_Sort)
  Day_AMRGene_Comb <- c(Day_AMRGene_Comb,Day_AMRGene_Sort )
}     
View(Day_AMRGene_Comb )

amr_Placebo_ST_Day_Processed$Day013090_AMRGene <- factor(amr_Placebo_ST_Day_Processed$Day013090_AMRGene, levels=unique(Day_AMRGene_Comb))  #unique(Day_AMRGene_Comb)
unique(sort(amr_Placebo_ST_Day_Processed$Day013090_AMRGene)) #45
unique(amr_Placebo_ST_Day_Processed$Day013090_AMRGene)
View(amr_Placebo_ST_Day_Processed)

BoxPlot_P_ST <-ggplot(amr_Placebo_ST_Day_Processed,aes(x=Day013090_AMRGene, y=Abundance , fill=Day013090 ))+ 
  stat_boxplot(geom='errorbar', width=0.5, lwd=0.5) + geom_boxplot(color="black", varwidth=FALSE, lwd=0.5)+  # geom_boxplot(fill="white")  #lwd=1.2  To make box lines thicker
  scale_y_continuous(expand=c(0,0), limits=c(MinYaxis, MaxYaxis))+    
  scale_fill_manual(values=rep(c("#009900", "#E69F00", "#FF3399"),15), name="Days after treatment")+  # This needs "aes(fill=Day013090 )"
  labs(title=paste0("High-frequency AMR gene abundance in Placebo group between Day01 vs Day03 vs Day90, Stool"), x="",y="AMR Gene Abundance",face="bold")+
  theme_classic()
BoxPlot_P_ST

MyPlot_P_ST<-BoxPlot_P_ST + # geom_point(aes(fill=MacrophageRatio), pch=21,size=3, position=position_jitterdodge()) + 
  theme(axis.text.x=element_text(colour="Black",size=9,angle=90,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.0),
        axis.line.x=element_line(color="black", size =1.0),
        axis.ticks=element_line(colour="black",size=1.0),
        axis.ticks.length=unit(.22, "cm")) 

# theme(legend.position="none") + 
MyPlot_P_ST  
ggsave(MyPlot,width=8, height=8, dpi=300, filename=paste0("OutFractionBoxplot_Placebo_ST_Day01Day30Day90_W8H8.pdf"), useDingbats=FALSE)

################################### ST  Rifaximin  #########################################################
## Check the min and max of 'Abudnance' data because need to set the min and max in the y-axis of the boxplot. 
MinYaxis_R<-floor(min(as.numeric(amr_Rifaximin_ST_Day_Processed[,"Abundance"])));  print(MinYaxis_R) # 1
MaxYaxis_R<-round(max(as.numeric(amr_Rifaximin_ST_Day_Processed[,"Abundance"])), digits=1); print(MaxYaxis_R) #20778


## Let's change the X-axis (Day013090_AMRGene) order
AMRGene_Unique_R <- unique(amr_Rifaximin_ST_Day_Processed$AMRGene) %>% sort #(alphabetical order)
AMRGene_Unique_R
unique(amr_Rifaximin_ST_Day_Processed$AMRGene) #14

Day_AMRGene_Comb_R <- c()
for(EachAMRGene in AMRGene_Unique_R) { 
  # EachAMRGene <- AMRGene_Unique[1]
  amr_Rifaximin_ST_Day_Subset <- amr_Rifaximin_ST_Day_Processed %>% dplyr::filter(AMRGene == EachAMRGene); a<-table(amr_Rifaximin_ST_Day_Subset$AMRGene)
  #print(a)
  Day_AMRGene_Sort_R <- sort(amr_Rifaximin_ST_Day_Subset$Day013090_AMRGene) #sort:d01>d30>d90
  print(Day_AMRGene_Sort_R)
  Day_AMRGene_Comb_R <- c(Day_AMRGene_Comb_R,Day_AMRGene_Sort_R )
  #print(Day_AMRGene_Comb)
}     
print(Day_AMRGene_Comb_R)    
#length(Day_AMRGene_Comb_R) #590
amr_Rifaximin_ST_Day_Processed$Day013090_AMRGene <- factor(amr_Rifaximin_ST_Day_Processed$Day013090_AMRGene, levels=unique(Day_AMRGene_Comb_R)) 
print(amr_Rifaximin_ST_Day_Processed)
View(amr_Rifaximin_ST_Day_Processed)
#u<-unique(Day_AMRGene_Comb_R)                             
#u #42

BoxPlot_R_ST <-ggplot(amr_Rifaximin_ST_Day_Processed,aes(x=Day013090_AMRGene, y=Abundance , fill=Day013090 ))+ 
  stat_boxplot(geom='errorbar', width=0.5, lwd=0.5) + geom_boxplot(color="black", varwidth=FALSE, lwd=0.5)+  # geom_boxplot(fill="white")  #lwd=1.2  To make box lines thicker
  scale_y_continuous(expand=c(0,0), limits=c(MinYaxis_R, MaxYaxis_R))+    
  scale_fill_manual(values=rep(c("#009900", "#E69F00", "#FF3399"),15), name="Days after treatment")+  # This needs "aes(fill=Day013090 )"
  labs(title=paste0("High-frequency AMR gene abundance in Rifaximin group between Day01 vs Day03 vs Day90, Stool"), x="",y="AMR Gene Abundance",face="bold")+
  theme_classic()
BoxPlot_R_ST

MyPlot_R_ST<-BoxPlot_R_ST + # geom_point(aes(fill=MacrophageRatio), pch=21,size=3, position=position_jitterdodge()) + 
  theme(axis.text.x=element_text(colour="Black",size=9,angle=90,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.0),
        axis.line.x=element_line(color="black", size =1.0),
        axis.ticks=element_line(colour="black",size=1.0),
        axis.ticks.length=unit(.22, "cm")) 
MyPlot_R_ST    
ggsave(MyPlot, height=5,width=8, dpi=300, filename=paste0("OutFractionBoxplot_R_ST_Day01Day30Day90.pdf"), useDingbats=FALSE)
#  theme(legend.position="none") +   # or pch=21 
######################################   SA Placebo   ####################################################################
MinYaxis<-floor(min(as.numeric(amr_Placebo_SA_Day_Processed[,"Abundance"])));  print(MinYaxis) # 1
MaxYaxis<-round(max(as.numeric(amr_Placebo_SA_Day_Processed[,"Abundance"])), digits=1); print(MaxYaxis)  #2960

## Let's change the X-axis (Day013090_AMRGene) order
AMRGene_Unique <- unique(amr_Placebo_SA_Day_Processed$AMRGene) %>% sort #alphabetical order
print(AMRGene_Unique) #14

Day_AMRGene_Comb <- c()
for(EachAMRGene in AMRGene_Unique) { 
  # EachAMRGene <- AMRGene_Unique[1]
  amr_Placebo_SA_Day_Subset <- amr_Placebo_SA_Day_Processed %>% dplyr::filter(AMRGene == EachAMRGene); table(amr_Placebo_SA_Day_Subset$AMRGene)
  Day_AMRGene_Sort <- sort(amr_Placebo_SA_Day_Subset$Day013090_AMRGene) #day 01->30->90 within a specific gene
  print(Day_AMRGene_Sort)
  Day_AMRGene_Comb <- c(Day_AMRGene_Comb,Day_AMRGene_Sort )
}     
#View(Day_AMRGene_Comb )

amr_Placebo_SA_Day_Processed$Day013090_AMRGene <- factor(amr_Placebo_SA_Day_Processed$Day013090_AMRGene, levels=unique(Day_AMRGene_Comb)) 
unique(sort(amr_Placebo_SA_Day_Processed$Day013090_AMRGene)) #42
View(amr_Placebo_SA_Day_Processed) 

BoxPlot_P_SA <-ggplot(amr_Placebo_SA_Day_Processed,aes(x=Day013090_AMRGene, y=Abundance , fill=Day013090 ))+ 
  stat_boxplot(geom='errorbar', width=0.5, lwd=0.5) + geom_boxplot(color="black", varwidth=FALSE, lwd=0.5)+  # geom_boxplot(fill="white")  #lwd=1.2  To make box lines thicker
  scale_y_continuous(expand=c(0,0), limits=c(MinYaxis, MaxYaxis))+    
  scale_fill_manual(values=rep(c("#009900", "#E69F00", "#FF3399"),15), name="Days after treatment")+  # This needs "aes(fill=Day013090 )"
  labs(title=paste0("High-frequency AMR gene abundance in Placebo group between Day01 vs Day03 vs Day90, Saliva"), x="",y="AMR Gene Abundance",face="bold")+
  theme_classic()
BoxPlot_P_SA

MyPlot_P_SA<-BoxPlot_P_SA + # geom_point(aes(fill=MacrophageRatio), pch=21,size=3, position=position_jitterdodge()) + 
  theme(axis.text.x=element_text(colour="Black",size=9,angle=90,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.0),
        axis.line.x=element_line(color="black", size =1.0),
        axis.ticks=element_line(colour="black",size=1.0),
        axis.ticks.length=unit(.22, "cm")) 
MyPlot_P_SA
ggsave(MyPlot,width=8, height=8, dpi=300, filename=paste0("OutFractionBoxplot_Placebo_SA_Day01Day30Day90_W8H8.pdf"), useDingbats=FALSE)
#  theme(legend.position="none") +   # or pch=21  
###################################    SA Rifaximin  #########################################################
## Check the min and max of 'Abudnance' data because need to set the min and max in the y-axis of the boxplot. 
MinYaxis_R<-floor(min(as.numeric(amr_Rifaximin_SA_Day_Processed[,"Abundance"])));  print(MinYaxis_R) # 1
MaxYaxis_R<-round(max(as.numeric(amr_Rifaximin_SA_Day_Processed[,"Abundance"])), digits=1); print(MaxYaxis_R) #4127

## Let's change the X-axis (Day013090_AMRGene) order
AMRGene_Unique_R <- unique(amr_Rifaximin_SA_Day_Processed$AMRGene) %>% sort #(alphabetical order)
AMRGene_Unique_R
unique(amr_Rifaximin_ST_Day_Processed$AMRGene) #14

Day_AMRGene_Comb_R <- c()
for(EachAMRGene in AMRGene_Unique_R) { 
  # EachAMRGene <- AMRGene_Unique[1]
  amr_Rifaximin_SA_Day_Subset <- amr_Rifaximin_SA_Day_Processed %>% dplyr::filter(AMRGene == EachAMRGene); a<-table(amr_Rifaximin_SA_Day_Subset$AMRGene)
  print(a)
  Day_AMRGene_Sort_R <- sort(amr_Rifaximin_SA_Day_Subset$Day013090_AMRGene) #sort:d01>d30>d90
  print(Day_AMRGene_Sort_R)
  Day_AMRGene_Comb_R <- c(Day_AMRGene_Comb_R,Day_AMRGene_Sort_R )
  #print(Day_AMRGene_Comb)
}     
print(Day_AMRGene_Comb_R)    
#length(Day_AMRGene_Comb_R) #601
amr_Rifaximin_SA_Day_Processed$Day013090_AMRGene <- factor(amr_Rifaximin_SA_Day_Processed$Day013090_AMRGene, levels=unique(Day_AMRGene_Comb_R)) 
print(amr_Rifaximin_SA_Day_Processed)
View(amr_Rifaximin_SA_Day_Processed)
u<-unique(Day_AMRGene_Comb_R)
u #42

BoxPlot_R_SA <-ggplot(amr_Rifaximin_SA_Day_Processed,aes(x=Day013090_AMRGene, y=Abundance , fill=Day013090 ))+ 
  stat_boxplot(geom='errorbar', width=0.5, lwd=0.5) + geom_boxplot(color="black", varwidth=FALSE, lwd=0.5)+  # geom_boxplot(fill="white")  #lwd=1.2  To make box lines thicker
  scale_y_continuous(expand=c(0,0), limits=c(MinYaxis_R, MaxYaxis_R))+    
  scale_fill_manual(values=rep(c("#009900", "#E69F00", "#FF3399"),15), name="Days after treatment")+  # This needs "aes(fill=Day013090 )"
  labs(title=paste0("High-frequency AMR gene abundance in Rifaximin group between Day01 vs Day03 vs Day90, Saliva"), x="",y="AMR Gene Abundance",face="bold")+
  theme_classic()
BoxPlot_R_SA

MyPlot_R_SA<-BoxPlot_R_SA + # geom_point(aes(fill=MacrophageRatio), pch=21,size=3, position=position_jitterdodge()) + 
  theme(axis.text.x=element_text(colour="Black",size=9,angle=90,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", size=1.0),
        axis.line.x=element_line(color="black", size =1.0),
        axis.ticks=element_line(colour="black",size=1.0),
        axis.ticks.length=unit(.22, "cm")) 
MyPlot_R_SA    
ggsave(MyPlot, height=5,width=8, dpi=300, filename=paste0("OutFractionBoxplot_R_SA_Day01Day30Day90.pdf"), useDingbats=FALSE)
#  theme(legend.position="none") +   # or pch=21
############################################################################################
## Step7. t-test between d01 vs d30 or d01 vs d90.
############################################################################################
### WRST function
func_WRST<-function(listName, WRST_pair_list, TidyrData) {
  print(paste0("listName: ", listName))
  print(paste0("pairList1: ",WRST_pair_list[1]))
  print(paste0("pairList2: ",WRST_pair_list[2]))
  WRST_x<-TidyrData$Abundance[TidyrData$Day013090_AMRGene==WRST_pair_list[1]]
  WRST_y<-TidyrData$Abundance[TidyrData$Day013090_AMRGene==WRST_pair_list[2]]
  WRST_summary<-wilcox.test(WRST_x, WRST_y, paired=FALSE, alternative="two.sided")
  WRST_pval<-WRST_summary$p.value
  print(paste0("WRST pval:", WRST_pval))
  return(WRST_pval)
}

###########################################  Placebo ST   ######################################################
View(amr_Placebo_ST_Day_Processed)
unique(sort(amr_Placebo_ST_Day_Processed$AMRGene))

WRST_pair_list<-list("d01_acrB_vs_d30_acrB"=c("d01_acrB","d30_acrB"),
                      "d01_acrB_vs_d90_acrB"=c("d01_acrB","d90_acrB"),
                      
                      "d01_adeF_vs_d30_adeF"=c("d01_adeF","d30_adeF"),
                      "d01_adeF_vs_d90_adeF"=c("d01_adeF","d90_adeF"),
                      
                      "d01_BahA_vs_d30_BahA"=c("d01_BahA","d30_BahA"),
                      "d01_BahA_vs_d90_BahA"=c("d01_BahA","d90_BahA"),
                      
                      "d01_CblA_1_vs_d30_CblA_1"=c("d01_CblA_1","d30_CblA_1"),
                      "d01_CblA_1_vs_d90_CblA_1"=c("d01_CblA_1","d90_CblA_1"),
                      
                      "d01_CfxA3_vs_d30_CfxA3"=c("d01_CfxA3","d30_CfxA3"),
                      "d01_CfxA3_vs_d90_CfxA3"=c("d01_CfxA3","d90_CfxA3"),
                      
                      "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                      =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                      "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                      =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                      
                      "d01_mdtB_vs_d30_mdtB"=c("d01_mdtB","d30_mdtB"),
                      "d01_mdtB_vs_d90_mdtB"=c("d01_mdtB","d90_mdtB"),
                      
                      "d01_mdtC_vs_d30_mdtC"=c("d01_mdtC","d30_mdtC"),
                      "d01_mdtC_vs_d90_mdtC"=c("d01_mdtC","d90_mdtC"),
                      
                      "d01_mefA_vs_d30_mefA"=c("d01_mefA","d30_mefA"),
                      "d01_mefA_vs_d90_mefA"=c("d01_mefA","d90_mefA"),
                      
                      "d01_patA_vs_d30_patA"=c("d01_patA","d30_patA"),
                      "d01_patA_vs_d90_patA"=c("d01_patA","d90_patA"),
                      
                      "d01_poxtA_vs_d30_poxtA"=c("d01_poxtA","d30_poxtA"),
                      "d01_poxtA_vs_d90_poxtA"=c("d01_poxtA","d90_poxtA"),
                      
                      "d01_tetA46_vs_d30_tetA46"=c("d01_tetA46","d30_tetA46"),
                      "d01_tetA46_vs_d90_tetA46"=c("d01_tetA46","d90_tetA46"),
                     
                      "d01_tetQ_vs_d30_tetQ"=c("d01_tetQ","d30_tetQ"),
                      "d01_tetQ_vs_d90_tetQ"=c("d01_tetQ","d90_tetQ"),
                      
                      "d01_tetW_vs_d30_tetW"=c("d01_tetW","d30_tetW"),
                      "d01_tetW_vs_d90_tetW"=c("d01_tetW","d90_tetW"),
                      
                      "d01_vanI_vs_d30_vanI"=c("d01_vanI","d30_vanI"),
                      "d01_vanI_vs_d90_vanI"=c("d01_vanI","d90_vanI"))


#amr_Placebo_ST_Day_Processed$Abundance[amr_Placebo_ST_Day_Processed$Day013090_AMRGene=="d01_adeF"]
#names(ttest_pair_list)
#ttest_pair_list[1]

#TidyrData<-list(amr_Placebo_ST_Day_Processed)
#View(TidyrData)
WRST_summary<-mapply(func_WRST,listName=names(WRST_pair_list), WRST_pair_list=WRST_pair_list, 
                      TidyrData=list(amr_Placebo_ST_Day_Processed))
WRST_summary
OutFileName_WRST_GenSig_ST_P <- "WRSTSummary_Placebo_ST_Day01Day30Day90_15ARMGene.txt"
write.table(WRST_summary, OutFileName_WRST_GenSig_ST_P,col.names=F, sep="\t", quote=F)


###################################################  Rifaximin  ST  #####################################################
### t-test function                                                        
#View(amr_Rifaximin_ST_Day_Processed)
unique(sort(amr_Rifaximin_ST_Day_Processed$AMRGene))


WRST_pair_list_R_ST <-list( "d01_acrD_vs_d30_acrD"=c("d01_acrD","d30_acrD"),
                            "d01_acrD_vs_d90_acrD"=c("d01_acrD","d90_acrD"),
  
                            "d01_adeF_vs_d30_adeF"=c("d01_adeF","d30_adeF"),
                            "d01_adeF_vs_d90_adeF"=c("d01_adeF","d90_adeF"),
  
                            "d01_emrY_vs_d30_emrY"=c("d01_emrY","d30_emrY"),
                            "d01_emrY_vs_d90_emrY"=c("d01_emrY","d90_emrY"),
  
                            "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                            =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
  
                            "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                            =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
  
                            "d01_macA_vs_d30_macA"=c("d01_macA","d30_macA"),
                            "d01_macA_vs_d90_macA"=c("d01_macA","d90_macA"),
                         
                         "d01_mdtC_vs_d30_mdtC"=c("d01_mdtC","d30_mdtC"),
                         "d01_mdtC_vs_d90_mdtC"=c("d01_mdtC","d90_mdtC"),
             
                         "d01_mdtO_vs_d30_mdtO"=c("d01_mdtO","d30_mdtO"),
                         "d01_mdtO_vs_d90_mdtO"=c("d01_mdtO","d90_mdtO"),            
                         
                         "d01_mefA_vs_d30_mefA"=c("d01_mefA","d30_mefA"),
                         "d01_mefA_vs_d90_mefA"=c("d01_mefA","d90_mefA"),

                         "d01_patA_vs_d30_patA"=c("d01_patA","d30_patA"),
                         "d01_patA_vs_d90_patA"=c("d01_patA","d90_patA"),
                         
                         "d01_poxtA_vs_d30_poxtA"=c("d01_poxtA","d30_poxtA"),
                         "d01_poxtA_vs_d90_poxtA"=c("d01_poxtA","d90_poxtA"),
                         
                         "d01_tetA46_vs_d30_tetA46"=c("d01_tetA46","d30_tetA46"),
                         "d01_tetA46_vs_d90_tetA46"=c("d01_tetA46","d90_tetA46"),
                         
                         "d01_tetQ_vs_d30_tetQ"=c("d01_tetQ","d30_tetQ"),
                         "d01_tetQ_vs_d90_tetQ"=c("d01_tetQ","d90_tetQ"),
                         
                         "d01_tetW_vs_d30_tetW"=c("d01_tetW","d30_tetW"),
                         "d01_tetW_vs_d90_tetW"=c("d01_tetW","d90_tetW"),
                         
                         "d01_vanI_vs_d30_vanI"=c("d01_vanI","d30_vanI"),
                         "d01_vanI_vs_d90_vanI"=c("d01_vanI","d90_vanI"))


WRST_summary_R_ST<-mapply(func_WRST,listName=names(WRST_pair_list_R_ST), WRST_pair_list=WRST_pair_list_R_ST, 
                        TidyrData=list(amr_Rifaximin_ST_Day_Processed))
WRST_summary_R_ST  

OutFileName_WRST_GenSig_R_ST <- "WRSTSummary_Rifaximin_ST_Day01Day30Day90_14ARMGene.txt"
write.table(WRST_summary_R_ST, OutFileName_WRST_GenSig_R_ST,col.names=F, sep="\t", quote=F)

##################################   Placebo SA  ####################################################
#View(amr_Placebo_SA_Day_Processed)
unique(sort(amr_Placebo_SA_Day_Processed$AMRGene))

WRST_pair_list_P_SA <-list("d01_adeF_vs_d30_adeF"=c("d01_adeF","d30_adeF"),
                      "d01_adeF_vs_d90_adeF"=c("d01_adeF","d90_adeF"),
                      
                      "d01_AxyY_vs_d30_AxyY"=c("d01_AxyY","d30_AxyY"),
                      "d01_AxyY_vs_d90_AxyY"=c("d01_AxyY","d90_AxyY"),
                      
                      "d01_CfxA3_vs_d30_CfxA3"=c("d01_CfxA3","d30_CfxA3"),
                      "d01_CfxA3_vs_d90_CfxA3"=c("d01_CfxA3","d90_CfxA3"),
                      
                      "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                      =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                      "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                      =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                      
                      "d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics_vs_d30_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"
                      =c("d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics","d30_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"),
                      "d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics_vs_d90_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"
                      =c("d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics","d90_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"),
                      
                      "d01_mefA_vs_d30_mefA"=c("d01_mefA","d30_mefA"),
                      "d01_mefA_vs_d90_mefA"=c("d01_mefA","d90_mefA"),
                      
                      "d01_mphB_vs_d30_mphB"=c("d01_mphB","d30_mphB"),
                      "d01_mphB_vs_d90_mphB"=c("d01_mphB","d90_mphB"),
                      
                      "d01_OXA_390_vs_d30_OXA_390"=c("d01_OXA_390","d30_OXA_390"),
                      "d01_OXA_390_vs_d90_OXA_390"=c("d01_OXA_390","d90_OXA_390"),
                      
                      "d01_pmrA_vs_d30_pmrA"=c("d01_pmrA","d30_pmrA"),
                      "d01_pmrA_vs_d90_pmrA"=c("d01_pmrA","d90_pmrA"),
  
                      "d01_tetA46_vs_d30_tetA46"=c("d01_tetA46","d30_tetA46"),
                      "d01_tetA46_vs_d90_tetA46"=c("d01_tetA46","d90_tetA46"),
                      
                      "d01_tetB46_vs_d30_tetB46"=c("d01_tetB46","d30_tetB46"),
                      "d01_tetB46_vs_d90_tetB46"=c("d01_tetB46","d90_tetB46"),
                      
                      "d01_tetB60_vs_d30_tetB60"=c("d01_tetB60","d30_tetB60"),
                      "d01_tetB60_vs_d90_tetB60"=c("d01_tetB60","d90_tetB60"),
                      
                      "d01_tetQ_vs_d30_tetQ"=c("d01_tetQ","d30_tetQ"),
                      "d01_tetQ_vs_d90_tetQ"=c("d01_tetQ","d90_tetQ"),
                      
                      "d01_vanTG_vs_d30_vanTG"=c("d01_vanTG","d30_vanTG"),
                      "d01_vanTG_vs_d90_vanTG"=c("d01_vanTG","d90_vanTG"))


#amr_Placebo_ST_Day_Processed$Abundance[amr_Placebo_ST_Day_Processed$Day013090_AMRGene=="d01_adeF"]
#names(ttest_pair_list)
#ttest_pair_list[1]
#unique(amr_Placebo_ST_Day_Processed$AMRGene)
View(TidyrData)
WRST_summary_P_SA<-mapply(func_WRST,listName=names(WRST_pair_list_P_SA), WRST_pair_list=WRST_pair_list_P_SA, 
                      TidyrData=list(amr_Placebo_SA_Day_Processed))
WRST_summary_P_SA
OutFileName_WRST_GenSig_SA_P <- "WRSTSummary_Placebo_SA_Day01Day30Day90_14ARMGene.txt"
write.table(WRST_summary_P_SA, OutFileName_WRST_GenSig_SA_P,col.names=F, sep="\t", quote=F)


##########################################  Rifaximin SA   #############################################################
#View(amr_Rifaximin_SA_Day_Processed)
unique(amr_Rifaximin_SA_Day_Processed$AMRGene)

WRST_pair_list_R_SA <-list("d01_adeF_vs_d30_adeF"=c("d01_adeF","d30_adeF"),
                           "d01_adeF_vs_d90_adeF"=c("d01_adeF","d90_adeF"),
                           
                           "d01_AxyY_vs_d30_AxyY"=c("d01_AxyY","d30_AxyY"),
                           "d01_AxyY_vs_d90_AxyY"=c("d01_AxyY","d90_AxyY"),
                           
                           "d01_CfxA3_vs_d30_CfxA3"=c("d01_CfxA3","d30_CfxA3"),
                           "d01_CfxA3_vs_d90_CfxA3"=c("d01_CfxA3","d90_CfxA3"),
                           
                           "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                           =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d30_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                           "d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin_vs_d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"
                           =c("d01_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin","d90_Escherichia_coli_EF_Tu_mutants_conferring_resistance_to_Pulvomycin"),
                           
                           "d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics_vs_d30_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"
                           =c("d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics","d30_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"),
                           "d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics_vs_d90_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"
                           =c("d01_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics","d90_Haemophilus_influenzae_PBP3_conferring_resistance_to_beta_lactam_antibiotics"),
                           
                           "d01_mefA_vs_d30_mefA"=c("d01_mefA","d30_mefA"),
                           "d01_mefA_vs_d90_mefA"=c("d01_mefA","d90_mefA"),
                           
                           "d01_msbA_vs_d30_msbA"=c("d01_msbA","d30_msbA"),
                           "d01_msbA_vs_d90_msbA"=c("d01_msbA","d90_msbA"),
                           
                           "d01_OXA_390_vs_d30_OXA_390"=c("d01_OXA_390","d30_OXA_390"),
                           "d01_OXA_390_vs_d90_OXA_390"=c("d01_OXA_390","d90_OXA_390"),
                           
                           "d01_patB_vs_d30_patB"=c("d01_patB","d30_patB"),
                           "d01_patB_vs_d90_patB"=c("d01_patB","d90_patB"),
                           
                           "d01_pmrA_vs_d30_pmrA"=c("d01_pmrA","d30_pmrA"),
                           "d01_pmrA_vs_d90_pmrA"=c("d01_pmrA","d90_pmrA"),
                           
                           "d01_tetA46_vs_d30_tetA46"=c("d01_tetA46","d30_tetA46"),
                           "d01_tetA46_vs_d90_tetA46"=c("d01_tetA46","d90_tetA46"),
                           
                           "d01_tetB46_vs_d30_tetB46"=c("d01_tetB46","d30_tetB46"),
                           "d01_tetB46_vs_d90_tetB46"=c("d01_tetB46","d90_tetB46"),
                           
                           "d01_tetG_vs_d30_tetG"=c("d01_tetG","d30_tetG"),
                           "d01_tetG_vs_d90_tetG"=c("d01_tetG","d90_tetG"),
                           
                           "d01_tetQ_vs_d30_tetQ"=c("d01_tetQ","d30_tetQ"),
                           "d01_tetQ_vs_d90_tetQ"=c("d01_tetQ","d90_tetQ"))


#amr_Placebo_ST_Day_Processed$Abundance[amr_Placebo_ST_Day_Processed$Day013090_AMRGene=="d01_adeF"]
#names(ttest_pair_list)
#ttest_pair_list[1]
#unique(amr_Placebo_ST_Day_Processed$AMRGene)
View(TidyrData)
WRST_summary_R_SA<-mapply(func_WRST,listName=names(WRST_pair_list_R_SA), WRST_pair_list=WRST_pair_list_R_SA, 
                          TidyrData=list(amr_Rifaximin_SA_Day_Processed))
WRST_summary_R_SA
OutFileName_WRST_GenSig_SA_R <- "WRSTSummary_Rifaximin_SA_Day01Day30Day90.txt"
write.table(WRST_summary_R_SA, OutFileName_WRST_GenSig_SA_R,col.names=F, sep="\t", quote=F)

################   ########################     ##################     ##################################    #######################
############################################################################################
## Step8. Prepare lineplot input.    
## Lecture: https://urldefense.com/v3/__http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization__;!!NHLzug!Oej5-ovJpXHABep3nXT1D7tNoh5smwRQ12wY10BHcuO85H8699IOiIeDqx2mSKi9r3xbh1UUZuk04v4pEJHmDSE7$ 
############################################################################################

## Function to make lineplot input
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
??ddply

## Start from the line 31
amr_Subtype <- amrData %>% dplyr::mutate(Subtype=paste0(RifaPlacebo,"_",SAST ) )
table(amr_Subtype$Subtype)
View(amr_Subtype)
# P_SA P_ST R_SA R_ST 
# 47   42   48   48

## Lineplot will display a specific gene's abudance difference betwee Placebo and Rifaximin samples across day01, day30, day90
## need to extract P_ST and R_ST.  - I will make a lineplot on Stool data first. You can make it on Saliva data later. 
amr_Subtype_ST <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01|d30|d90", ID));
table(amr_Subtype_ST$Subtype) # P_ST:42, R_ST:48
amr_Subtype_ST$Subtype

## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_ST$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
amr_Subtype_ST_Day <- amr_Subtype_ST %>% dplyr::mutate(Day=Day)
table(amr_Subtype_ST_Day$Day)  # d01:36,  d30:29, d90:25 

## Select genes of interest to make lineplots. 
MyInterestGene <- c("tetQ")    #   see the volcanoplot again to select intersting genes

for(MyGene in MyInterestGene) {
  #MyGene <- MyInterestGene[1]
  amr_Subtype_ST_Day_Gene <- amr_Subtype_ST_Day %>% dplyr::select(print(MyGene), Subtype, Day)
  LineplotInput <- data_summary(amr_Subtype_ST_Day_Gene, varname=MyGene, 
                                groupnames=c("Subtype", "Day"))
  colnames(LineplotInput)[3]<-"Abundance"
  
  ## Find Y-axis max value
  Yaxis_Max <- max(LineplotInput$Abundance)
  
  # Use position_dodge to move overlapped errorbars horizontally
  MyLineplot <- ggplot(LineplotInput, aes(x=Day, y=Abundance, group=Subtype, shape=Subtype, linetype=Subtype, color=Subtype)) + 
    geom_errorbar(aes(ymin=Abundance, ymax=Abundance+sd), width=.3,position=position_dodge(0.05)) +
    geom_line(size=2) + geom_point(color="black", size=3)+  ylim(c(0, Yaxis_Max+0.2)) +
    #scale_color_brewer(palette="Paired")+
    labs(title=paste0("Plot of ",MyGene," abundance between Placebo and Rifaximin"),x="Day0, Day30, Day90", y="AMR gene abaudance")+
    theme_classic()+scale_color_manual(values=c('green','red'))
  
  ## Print MyLineplot in the Rstudio
  print(MyLineplot)  
  
  ##  save the fraction barplot to PDF file. 
  ggsave(MyLineplot, file=paste0("Lineplot_PlaceboRifaximin_ST_",MyGene,"_W5H5_20220517_v1.0.pdf"), width=5, height=5)
  
}

print(LineplotInput)
print(amr_Subtype_ST_Day)


