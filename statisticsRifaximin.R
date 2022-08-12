dim(amr_Subtype)  # 185  588
tail(colnames(amr_Subtype))
View(amr_Subtype)

#########  normal distribution test:Shapiro-Wilk test   ST, Shannon

amr_Subtype_ST<- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_ST$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_ST_Day <- amr_Subtype_ST %>% dplyr::mutate(Day=Day)
table(amr_Subtype_ST_Day$Day)  # d01:36,  d30:29, d90:25 
View(amr_Subtype_ST_Day)

amr_Subtype_ST_Day <- amr_Subtype_ST_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_ST_Day))                      
View(amr_Subtype_ST_Day)

############################    by subtypeday    ##########################################

ST_shapiro_shannon <-shapiro.test(amr_Subtype_ST_Day$Shannon)
ST_shapiro_shannon   # p-value = 2.275e-07 ---> non - normal
kruskal.test(Shannon~ Subtypeday, data=amr_Subtype_ST_Day) #p-value = 0.2696  --> non significant
#shannonaovST<-aov(Shannon~ Subtype +Day, data=amr_Subtype_ST_Day)
#summary(shannonaovST)   # Subtype  0.0558 / Day  0.6022  ---> non significant

##########################   by day01 vs day30 vs day90  ###############################

ST_shapiro_shannon <-shapiro.test(amr_Subtype_ST_Day$Shannon)
ST_shapiro_shannon  # p-value = 2.275e-07   ---> non- normal
kruskal.test(Shannon~ Day, data=amr_Subtype_ST_Day)  # p-value = 0.6478  --> non significant
#shannonaovSTday<-aov(Shannon~ Day, data=amr_Subtype_ST_Day)
#summary(shannonaovSTday)  #0.675----->non significant

####################   day01, day30, day90 in placebo  ##################################

View(amr_Subtype_ST_Day)
tail(colnames(amr_Subtype_ST_Day))

amr_Subtype_ST_Day_P<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_P)

ST_shapiro_shannon_P <-shapiro.test(amr_Subtype_ST_Day_P$Shannon)
ST_shapiro_shannon_P   # p-value = 3.647e-06 ---> non-normal

kruskal.test(Shannon~ Day, data=amr_Subtype_ST_Day_P) # p-value = 0.535  --> non significant

#shannonaovST_Pday<-aov(Shannon~ Day, data=amr_Subtype_ST_Day_P)
#summary(shannonaovST_Pday)   #0.48 ---> non significant

####################   day01, day30, day90 in rifaximin  ################################

amr_Subtype_ST_Day_R<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("R_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_R)

ST_shapiro_shannon_R <-shapiro.test(amr_Subtype_ST_Day_R$Shannon)
ST_shapiro_shannon_R   # p-value = 0.001654 ---> non-normal

kruskal.test(Shannon~ Day, data=amr_Subtype_ST_Day_R) #p-value = 0.9725 ---> non significant

#shannonaovST_Rday<-aov(Shannon~ Day, data=amr_Subtype_ST_Day_R)
#summary(shannonaovST_Rday)  #  0.949--> non significant

########################  placebo vs Rifaximin - day01, day30, day90   ##################

amr_Subtype_ST_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01", ID)); #36
View(amr_Subtype_ST_d01)
amr_Subtype_ST_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d30", ID)); #29
amr_Subtype_ST_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d90", ID)); #25

ST_shapiro_shannon_d01 <-shapiro.test(amr_Subtype_ST_d01$Shannon)
ST_shapiro_shannon_d01   #p-value = 0.0009864 --> non normal
ST_shapiro_shannon_d30 <-shapiro.test(amr_Subtype_ST_d30$Shannon)
ST_shapiro_shannon_d30   # p-value = 0.0006557  --> non- normal
ST_shapiro_shannon_d90 <-shapiro.test(amr_Subtype_ST_d90$Shannon)
ST_shapiro_shannon_d90   # p-value = 0.006636  --> non- normal

tail(colnames(amr_Subtype_ST_Day))

index_1_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_P
index_2_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_P
index_3_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_P  

index_1_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_R
index_2_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_R
index_3_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_R  

#######  d01  ######
shannon_1_ST_P<-amr_Subtype_ST_Day$Shannon[index_1_ST_P]
shannon_1_ST_R<-amr_Subtype_ST_Day$Shannon[index_1_ST_R]
wilcox_result_1<- wilcox.test(shannon_1_ST_P,shannon_1_ST_R)
wilcox_result_1   # p-value = 0.4429 --> not significant

#######  d30  #######
shannon_2_ST_P<-amr_Subtype_ST_Day$Shannon[index_2_ST_P]
shannon_2_ST_R<-amr_Subtype_ST_Day$Shannon[index_2_ST_R]
wilcox_result_2<-wilcox.test(shannon_2_ST_P,shannon_2_ST_R)
wilcox_result_2   #  p-value = 0.1318 --> not significant

#######  d90  ########
shannon_3_ST_P<-amr_Subtype_ST_Day$Shannon[index_3_ST_P]
shannon_3_ST_R<-amr_Subtype_ST_Day$Shannon[index_3_ST_R]
wilcox_result_3<-wilcox.test(shannon_3_ST_P,shannon_3_ST_R)
wilcox_result_3   #  p-value = 0.1656 --> not significant

#############################   SA shannon   ##########################################
######### normal distribution test:Shapiro-Wilk test   SA, Shannon
amr_Subtype_SA<- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_SA$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_SA_Day <- amr_Subtype_SA %>% dplyr::mutate(Day=Day)
table(amr_Subtype_SA_Day$Day)  # d01:37,  d30:32, d90:26 
View(amr_Subtype_SA_Day)

amr_Subtype_SA_Day <- amr_Subtype_SA_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_SA_Day))                      
View(amr_Subtype_SA_Day)

############################  SA Shannon by subtypeday   ##################################

SA_shapiro_shannon <-shapiro.test(amr_Subtype_SA_Day$Shannon)
SA_shapiro_shannon   # p-value = 0.07246 --->  normal

subtypeaov <- aov(Shannon ~ Subtypeday, data = amr_Subtype_SA_Day)
subtypeaov2 <- aov(Shannon ~ Subtype+Day, data = amr_Subtype_SA_Day)
summary(subtypeaov)   # p-value 0.3  ---> non significant
summary(subtypeaov2)  # subtype 0.3536/ Day 0.0924
##########################   SA Shannon by day01 vs day30 vs day90    ###################

ST_shapiro_shannon <-shapiro.test(amr_Subtype_SA_Day$Shannon)
ST_shapiro_shannon  # p-value = 0.07246 --->  normal

dayaovSA <- aov(Shannon ~ Day, data = amr_Subtype_SA_Day)
summary(dayaovSA)    # p-value=  0.1 ---> non significant

####################  SA Shannon day01, day30, day90 in placebo  #########################

View(amr_Subtype_SA_Day)
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_Day_P<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_P)

SA_shapiro_shannon_P <-shapiro.test(amr_Subtype_SA_Day_P$Shannon)
SA_shapiro_shannon_P   # p-value = 0.1749 ---> normal

dayaovSA_P <- aov(Shannon ~ Day, data = amr_Subtype_SA_Day_P)
summary(dayaovSA_P)     # p-value 0.352  ---> non-significant

####################  SA Shannon day01, day30, day90 in rifaximin  ########################

amr_Subtype_SA_Day_R<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("R_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_R)

SA_shapiro_shannon_R <-shapiro.test(amr_Subtype_SA_Day_R$Shannon)
SA_shapiro_shannon_R   # p-value = 0.2937 ---> normal

dayaovSA_R <- aov(Shannon ~ Day, data = amr_Subtype_SA_Day_R)
summary(dayaovSA_R) # 0.239  ---> non significant

#################  SA Shannon  placebo vs Rifaximin - day01, day30, day90   ##############
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01", ID)); #37
amr_Subtype_SA_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d30", ID)); #32
amr_Subtype_SA_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d90", ID)); #26

SA_shapiro_shannon_d01 <-shapiro.test(amr_Subtype_SA_d01$Shannon)
SA_shapiro_shannon_d01   #p-value = 0.3967 --> normal
SA_shapiro_shannon_d30 <-shapiro.test(amr_Subtype_SA_d30$Shannon)
SA_shapiro_shannon_d30   # p-value = 0.03546  --> non- normal
SA_shapiro_shannon_d90 <-shapiro.test(amr_Subtype_SA_d90$Shannon)
SA_shapiro_shannon_d90   # p-value = 0.4596  --> normal

index_1_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_P
index_2_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_P
index_3_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_P  

index_1_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_R
index_2_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_R
index_3_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_R  

#######  d01  ######
shannon_1_SA_P<-amr_Subtype_SA_Day$Shannon[index_1_SA_P]
shannon_1_SA_R<-amr_Subtype_SA_Day$Shannon[index_1_SA_R]
t_result_1_SA<-t.test(shannon_1_SA_P,shannon_1_SA_R)
t_result_1_SA   #  p-value = 0.3788 --> not significant

#######  d30  ######
shannon_2_SA_P<-amr_Subtype_SA_Day$Shannon[index_2_SA_P]
shannon_2_SA_R<-amr_Subtype_SA_Day$Shannon[index_2_SA_R]
wilcox_result_2_SA<-wilcox.test(shannon_2_SA_P,shannon_2_SA_R)
wilcox_result_2_SA   # p-value = 0.2099 --> not significant

#######  d90  ######
shannon_3_SA_P<-amr_Subtype_SA_Day$Shannon[index_3_SA_P]
shannon_3_SA_R<-amr_Subtype_SA_Day$Shannon[index_3_SA_R]
t_result_3_SA<-t.test(shannon_3_SA_P,shannon_3_SA_R)
t_result_3_SA   #  p-value = 0.9397 --> not significant

######### normal distribution test:Shapiro-Wilk test   ST, richness

amr_Subtype_ST<- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_ST$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_ST_Day <- amr_Subtype_ST %>% dplyr::mutate(Day=Day)
table(amr_Subtype_ST_Day$Day)  # d01:36,  d30:29, d90:25 
View(amr_Subtype_ST_Day)

amr_Subtype_ST_Day <- amr_Subtype_ST_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_ST_Day))                      
View(amr_Subtype_ST_Day)

############################  ST, richness  by subtypeday   ###############################

ST_shapiro_richness <-shapiro.test(amr_Subtype_ST_Day$richness)
ST_shapiro_richness   # p-value = 0.8434 ---> normal

subtypeaovrich <- aov(richness ~ Subtypeday, data = amr_Subtype_ST_Day)
summary(subtypeaovrich)  # p-value = 0.589--->non significant
subtypeaovrich2 <- aov(richness ~ Subtype + Day, data = amr_Subtype_ST_Day)
summary(subtypeaovrich2)  #Subtype 0.171/ Day 0.774

###################  ST, richness by day01 vs day30 vs day90  #############################

ST_shapiro_richness <-shapiro.test(amr_Subtype_ST_Day$richness)
ST_shapiro_richness  # p-value = 0.8434  --->  normal

dayaovrich <- aov(richness ~ Day, data = amr_Subtype_ST_Day)
summary(dayaovrich) # p-value =  0.742 --->non significant

####################   ST, richness   day01, day30, day90 in placebo  ####################

View(amr_Subtype_ST_Day)
tail(colnames(amr_Subtype_ST_Day))

amr_Subtype_ST_Day_P<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_P)

ST_shapiro_richness_P <-shapiro.test(amr_Subtype_ST_Day_P$richness)
ST_shapiro_richness_P   # p-value = 0.4231 --->  normal

dayaovrich_ST_P <- aov(richness ~ Day, data = amr_Subtype_ST_Day_P)
summary(dayaovrich_ST_P)  # p-value  0.819---> non significant

####################  ST, richness  day01, day30, day90 in rifaximin  #####################

amr_Subtype_ST_Day_R<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("R_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_R)

ST_shapiro_richness_R <-shapiro.test(amr_Subtype_ST_Day_R$richness)
ST_shapiro_richness_R   # p-value = 0.9729 --->  normal

dayaovrich_ST_R <- aov(richness ~ Day, data = amr_Subtype_ST_Day_R)
summary(dayaovrich_ST_R)   # p-value =  0.539  ---> non significant

########################  placebo vs Rifaximin - day01, day30, day90   ##################

amr_Subtype_ST_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01", ID)); #36
amr_Subtype_ST_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d30", ID)); #29
amr_Subtype_ST_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d90", ID)); #25

ST_shapiro_richness_d01 <-shapiro.test(amr_Subtype_ST_d01$richness)
ST_shapiro_richness_d01   # p-value = 0.5262 --> normal
ST_shapiro_richness_d30 <-shapiro.test(amr_Subtype_ST_d30$richness)
ST_shapiro_richness_d30   # p-value = 0.8545 --> normal
ST_shapiro_richness_d90 <-shapiro.test(amr_Subtype_ST_d90$richness)
ST_shapiro_richness_d90   # p-value = 0.3685 --> normal

tail(colnames(amr_Subtype_ST_Day))

index_1_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_P
index_2_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_P
index_3_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_P  

index_1_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_R
index_2_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_R
index_3_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_R  

#######  d01  ######
richness_1_ST_P<-amr_Subtype_ST_Day$richness[index_1_ST_P]
richness_1_ST_R<-amr_Subtype_ST_Day$richness[index_1_ST_R]
t_result_1_ST<- t.test(richness_1_ST_P, richness_1_ST_R)
t_result_1_ST   # p-value = 0.9717 --> non significant

#######  d30  #######
richness_2_ST_P<-amr_Subtype_ST_Day$richness[index_2_ST_P]
richness_2_ST_R<-amr_Subtype_ST_Day$richness[index_2_ST_R]
t_result_2_ST<-t.test(richness_2_ST_P, richness_2_ST_R)
t_result_2_ST   #  p-value = 0.1992 --> not significant

#######  d90  ########
richness_3_ST_P<-amr_Subtype_ST_Day$richness[index_3_ST_P]
richness_3_ST_R<-amr_Subtype_ST_Day$richness[index_3_ST_R]
t_result_3_ST<-t.test(richness_3_ST_P,richness_3_ST_R)
t_result_3_ST   #  p-value = 0.1243 --> not significant

################################  SA richness #################################
######### normal distribution test:Shapiro-Wilk test   SA, richness
amr_Subtype_SA<- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_SA$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_SA_Day <- amr_Subtype_SA %>% dplyr::mutate(Day=Day)
table(amr_Subtype_SA_Day$Day)  # d01:36,  d30:29, d90:25 
View(amr_Subtype_SA_Day)

amr_Subtype_SA_Day <- amr_Subtype_SA_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_SA_Day))                      
View(amr_Subtype_SA_Day)

############################  SA richness  by subtypeday  ###################################

SA_shapiro_richness <-shapiro.test(amr_Subtype_SA_Day$richness)
SA_shapiro_richness   # p-value = 0.1194 --->  normal

subtypeaovrich <- aov(richness ~ Subtypeday, data = amr_Subtype_SA_Day)
summary(subtypeaovrich)   # p-value  0.105  ---> non significant

subtypeaovrich2 <- aov(richness ~ Subtype+Day, data = amr_Subtype_SA_Day)
summary(subtypeaovrich2)   # p-value  Subtype  0.4976/ Day 0.0129 *

######################  SA richness by day01 vs day30 vs day90  ##########################

ST_shapiro_richness <-shapiro.test(amr_Subtype_SA_Day$richness)
ST_shapiro_richness  #p-value = 0.1194 --->  normal

dayaovSArich <- aov(richness~ Day, data = amr_Subtype_SA_Day)
summary(dayaovSArich)    # p-value=  0.0133 *--->  significant
TukeyHSD(dayaovSArich, which = "Day") #d30-d01 0.0161905*** / d90-d01 0.0786603 / d90-d30 0.8920816

###############    SA richness  day01, day30, day90 in placebo   #########################

View(amr_Subtype_SA_Day)
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_Day_P<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_P)

SA_shapiro_richness_P <-shapiro.test(amr_Subtype_SA_Day_P$richness)
SA_shapiro_richness_P   #  p-value = 0.2734 ---> normal

dayaovSA_P_rich <- aov(richness ~ Day, data = amr_Subtype_SA_Day_P)
summary(dayaovSA_P_rich)     # p-value 0.0695   ---> non-significant

####################  SA richness day01, day30, day90 in rifaximin  ########################

amr_Subtype_SA_Day_R<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("R_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_R)

SA_shapiro_richness_R <-shapiro.test(amr_Subtype_SA_Day_R$richness)
SA_shapiro_richness_R  # p-value = 0.5874 ---> normal

dayaovSA_R_rich <- aov(richness ~ Day, data = amr_Subtype_SA_Day_R)
summary(dayaovSA_R_rich) # 0.171  ---> non significant

#################  SA richness  placebo vs Rifaximin - day01, day30, day90   ##################
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01", ID)); #37
amr_Subtype_SA_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d30", ID)); #32
amr_Subtype_SA_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d90", ID)); #26

SA_shapiro_richness_d01 <-shapiro.test(amr_Subtype_SA_d01$richness)
SA_shapiro_richness_d01   #p-value = 0.02693 --> non- normal
SA_shapiro_richness_d30 <-shapiro.test(amr_Subtype_SA_d30$richness)
SA_shapiro_richness_d30   # p-value = 0.6726 -->  normal
SA_shapiro_richness_d90 <-shapiro.test(amr_Subtype_SA_d90$richness)
SA_shapiro_richness_d90   # p-value = 0.5585  --> normal

index_1_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_P
index_2_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_P
index_3_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_P  

index_1_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_R
index_2_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_R
index_3_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_R  

#######  d01  ######
richness_1_SA_P<-amr_Subtype_SA_Day$richness[index_1_SA_P]
richness_1_SA_R<-amr_Subtype_SA_Day$richness[index_1_SA_R]
wilcox_1_SA_rich<-wilcox.test(richness_1_SA_P,richness_1_SA_R)
wilcox_1_SA_rich   # p-value = 0.7727 --> not significant
#######  d30  ######
richness_2_SA_P<-amr_Subtype_SA_Day$richness[index_2_SA_P]
richness_2_SA_R<-amr_Subtype_SA_Day$richness[index_2_SA_R]
t_2_SA_rich<- t.test(richness_2_SA_P,richness_2_SA_R)
t_2_SA_rich   # p-value = 0.6687 --> not significant
#######  d90  ######
richness_3_SA_P<-amr_Subtype_SA_Day$richness[index_3_SA_P]
richness_3_SA_R<-amr_Subtype_SA_Day$richness[index_3_SA_R]
t_3_SA_rich<-t.test(richness_3_SA_P, richness_3_SA_R)
t_3_SA_rich   # p-value = 0.5958 --> not significant

######### normal distribution test:Shapiro-Wilk test ST, abundance 
amr_Subtype_ST<- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_ST$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_ST_Day <- amr_Subtype_ST %>% dplyr::mutate(Day=Day)
table(amr_Subtype_ST_Day$Day)  # d01:36,  d30:29, d90:25 
View(amr_Subtype_ST_Day)

amr_Subtype_ST_Day <- amr_Subtype_ST_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_ST_Day))                      
View(amr_Subtype_ST_Day)

############################ ST abundance by subtypeday  ######################################

ST_shapiro_rowsum <-shapiro.test(amr_Subtype_ST_Day$rowsum)
ST_shapiro_rowsum   # p-value = 1.949e-15 ---> non - normal

kruskal.test(rowsum~ Subtypeday, data=amr_Subtype_ST_Day) #p-value = 0.9835  --> non significant

##########################  ST abundance by day01 vs day30 vs day90  ###########################

ST_shapiro_rowsum <-shapiro.test(amr_Subtype_ST_Day$rowsum)
ST_shapiro_rowsum  # p-value = 1.949e-15   ---> non- normal

kruskal.test(rowsum~ Day, data=amr_Subtype_ST_Day)  # p-value = 0.8582  --> non significant

########################## ST abundance day01, day30, day90 in placebo   #########################
View(amr_Subtype_ST_Day)
tail(colnames(amr_Subtype_ST_Day))

amr_Subtype_ST_Day_P<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("P_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_P)

ST_shapiro_rowsum_P <-shapiro.test(amr_Subtype_ST_Day_P$rowsum)
ST_shapiro_rowsum_P   # p-value = 2.053e-09 ---> non-normal

kruskal.test(rowsum~ Day, data=amr_Subtype_ST_Day_P) # p-value = 0.8557  --> non significant

######################   ST abundance day01, day30, day90 in rifaximin  ###########################

amr_Subtype_ST_Day_R<- amr_Subtype_ST_Day %>% dplyr::filter(grepl("R_ST", Subtype)  & grepl("d01|d30|d90", ID)); #42
View(amr_Subtype_ST_Day_R)

ST_shapiro_rowsum_R <-shapiro.test(amr_Subtype_ST_Day_R$rowsum)
ST_shapiro_rowsum_R   # p-value = 3.904e-11 ---> non-normal

kruskal.test(rowsum~ Day, data=amr_Subtype_ST_Day_R) #p-value = 0.9517 ---> non significant

########################  placebo vs Rifaximin - day01, day30, day90   ##################

amr_Subtype_ST_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d01", ID)); #36
amr_Subtype_ST_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d30", ID)); #29
amr_Subtype_ST_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_ST|R_ST", Subtype)  & grepl("d90", ID)); #25

ST_shapiro_rowsum_d01 <-shapiro.test(amr_Subtype_ST_d01$rowsum)
ST_shapiro_rowsum_d01   #p-value = 3.108e-09 --> non normal
ST_shapiro_rowsum_d30 <-shapiro.test(amr_Subtype_ST_d30$rowsum)
ST_shapiro_rowsum_d30   # p-value = 4.431e-08  --> non- normal
ST_shapiro_rowsum_d90 <-shapiro.test(amr_Subtype_ST_d90$rowsum)
ST_shapiro_rowsum_d90   # p-value = 1.547e-08  --> non- normal

tail(colnames(amr_Subtype_ST_Day))

index_1_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_P
index_2_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_P
index_3_ST_P<-which(amr_Subtype_ST_Day$Subtype=="P_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_P  

index_1_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d01")
index_1_ST_R
index_2_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d30") 
index_2_ST_R
index_3_ST_R<-which(amr_Subtype_ST_Day$Subtype=="R_ST" & amr_Subtype_ST_Day$Day=="d90") 
index_3_ST_R  

#######  d01  ######
rowsum_1_ST_P<-amr_Subtype_ST_Day$rowsum[index_1_ST_P]
rowsum_1_ST_R<-amr_Subtype_ST_Day$rowsum[index_1_ST_R]
wilcox_result_1<- wilcox.test(rowsum_1_ST_P,rowsum_1_ST_R)
wilcox_result_1   # p-value = 0.5628 --> not significant

#######  d30  #######
rowsum_2_ST_P<-amr_Subtype_ST_Day$rowsum[index_2_ST_P]
rowsum_2_ST_R<-amr_Subtype_ST_Day$rowsum[index_2_ST_R]
wilcox_result_2<-wilcox.test(rowsum_2_ST_P, rowsum_2_ST_R)
wilcox_result_2   #  p-value = 0.9484 --> not significant

#######  d90  ########
rowsum_3_ST_P<-amr_Subtype_ST_Day$rowsum[index_3_ST_P]
rowsum_3_ST_R<-amr_Subtype_ST_Day$rowsum[index_3_ST_R]
wilcox_result_3<-wilcox.test(rowsum_3_ST_P, rowsum_3_ST_R)
wilcox_result_3   #  p-value = 0.9786 --> not significant

################################   SA abundance   #################################
######### normal distribution test:Shapiro-Wilk test   SA, abundance
amr_Subtype_SA<- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01|d30|d90", ID))
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype_SA$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_SA_Day <- amr_Subtype_SA %>% dplyr::mutate(Day=Day)
table(amr_Subtype_SA_Day$Day)  # d01:36,  d30:29, d90:25 
View(amr_Subtype_SA_Day)

amr_Subtype_SA_Day <- amr_Subtype_SA_Day %>% dplyr::mutate(Subtypeday=paste0(Subtype,"_",Day ) )
tail(colnames(amr_Subtype_SA_Day))                      
View(amr_Subtype_SA_Day)

############################  SA  abundance by subtypeday ###############################

SA_shapiro_rowsum <-shapiro.test(amr_Subtype_SA_Day$rowsum)
SA_shapiro_rowsum   # p-value = p-value < 2.2e-16 --->  non - normal

kruskal.test(rowsum~ Subtypeday , data=amr_Subtype_SA_Day)  #p-value = 0.7725 ---> non- significant

#####################   SA  abundance   by day01 vs day30 vs day90  #####################

ST_shapiro_rowsum <-shapiro.test(amr_Subtype_SA_Day$rowsum)
ST_shapiro_rowsum  #  p-value < 2.2e-16 ---> non normal

kruskal.test(rowsum~ Day, data=amr_Subtype_SA_Day)  # p-value = 0.4813  ---> non- significant

####################   day01, day30, day90 in placebo  SA  ################################

View(amr_Subtype_SA_Day)
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_Day_P<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("P_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_P)

SA_shapiro_rowsum_P <- shapiro.test(amr_Subtype_SA_Day_P$rowsum)
SA_shapiro_rowsum_P   # p-value = 1.442e-14 ---> non normal

kruskal.test(rowsum~ Day, data=amr_Subtype_SA_Day_P)  # p-value = 0.915 ---> non- significant

####################   day01, day30, day90 in rifaximin SA   ##############################

amr_Subtype_SA_Day_R<- amr_Subtype_SA_Day %>% dplyr::filter(grepl("R_SA", Subtype)  & grepl("d01|d30|d90", ID)); 
View(amr_Subtype_SA_Day_R)

SA_shapiro_rowsum_R <-shapiro.test(amr_Subtype_SA_Day_R$rowsum)
SA_shapiro_rowsum_R   # p-value = 4.969e-07 ---> non normal

kruskal.test(rowsum~ Day, data=amr_Subtype_SA_Day_R)   # p-value = 0.3166 ---> non- significant

########################  placebo vs Rifaximin - day01, day30, day90   ####################
tail(colnames(amr_Subtype_SA_Day))

amr_Subtype_SA_d01 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d01", ID)); #37
amr_Subtype_SA_d30 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d30", ID)); #32
amr_Subtype_SA_d90 <- amr_Subtype %>% dplyr::filter(grepl("P_SA|R_SA", Subtype)  & grepl("d90", ID)); #26

SA_shapiro_rowsum_d01 <-shapiro.test(amr_Subtype_SA_d01$rowsum)
SA_shapiro_rowsum_d01   #p-value = 2.607e-05 --> non normal
SA_shapiro_rowsum_d30 <-shapiro.test(amr_Subtype_SA_d30$rowsum)
SA_shapiro_rowsum_d30   # p-value = 5.054e-06 --> non- normal
SA_shapiro_rowsum_d90 <-shapiro.test(amr_Subtype_SA_d90$rowsum)
SA_shapiro_rowsum_d90   # p-value = 1.301e-10  --> non- normal

index_1_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_P
index_2_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_P
index_3_SA_P<-which(amr_Subtype_SA_Day$Subtype=="P_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_P  

index_1_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d01")
index_1_SA_R
index_2_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d30") 
index_2_SA_R
index_3_SA_R<-which(amr_Subtype_SA_Day$Subtype=="R_SA" & amr_Subtype_SA_Day$Day=="d90") 
index_3_SA_R  

#######  d01  ######
rowsum_1_SA_P<-amr_Subtype_SA_Day$rowsum[index_1_SA_P]
rowsum_1_SA_R<-amr_Subtype_SA_Day$rowsum[index_1_SA_R]
w_result_1_SA<-wilcox.test(rowsum_1_SA_P,rowsum_1_SA_R)
w_result_1_SA   #  p-value = 0.5578 --> not significant

#######  d30  ######
rowsum_2_SA_P<-amr_Subtype_SA_Day$rowsum[index_2_SA_P]
rowsum_2_SA_R<-amr_Subtype_SA_Day$rowsum[index_2_SA_R]
w_result_2_SA<-wilcox.test(rowsum_2_SA_P,rowsum_2_SA_R)
w_result_2_SA   # p-value = 0.5147  --> not significant

#######  d90  ######
rowsum_3_SA_P<-amr_Subtype_SA_Day$rowsum[index_3_SA_P]
rowsum_3_SA_R<-amr_Subtype_SA_Day$rowsum[index_3_SA_R]
w_result_3_SA<-wilcox.test(rowsum_3_SA_P,rowsum_3_SA_R)
w_result_3_SA   # p-value = 0.7424 --> non signficant

################## richness compare P_ST vs R_ST vs P_SA vs R_SA ###################
## Make "Day" column 
Day <- sapply(stringr::str_split(amr_Subtype$ID, pattern="_"), "[", 3)  # "d01", "d30" or "d90"
Day
amr_Subtype_Day <- amr_Subtype %>% dplyr::mutate(Day=Day)
table(amr_Subtype_Day$Day)  # d01: 73,  d30: 61, d90: 51 
tail(colnames(amr_Subtype_Day))
View(amr_Subtype_Day)
unique(amr_Subtype_Day$Subtype)

shapiro_richness_P <-shapiro.test(amr_Subtype_Day$richness)
shapiro_richness_P  # p-value = 0.007731 ---> non normal

kruskal.test(richness~ Subtype , data=amr_Subtype_Day) #p-value < 2.2e-16 ---> significant
pairwise.wilcox.test(amr_Subtype_Day$richness, amr_Subtype_Day$Subtype)
?pairwise.wilcox.test

index_7_ST_P<-which(amr_Subtype_Day$Subtype=="P_ST")
index_7_ST_P
index_8_ST_R<-which(amr_Subtype_Day$Subtype=="R_ST") 
index_8_ST_R
index_9_SA_P<-which(amr_Subtype_Day$Subtype=="P_SA") 
index_9_SA_P  
index_10_SA_R<-which(amr_Subtype_Day$Subtype=="R_SA") 
index_10_SA_R 

richness_ST_P<-amr_Subtype_Day$richness[index_7_ST_P]
richness_ST_R<-amr_Subtype_Day$richness[index_8_ST_R]
richness_SA_P<-amr_Subtype_Day$richness[index_9_SA_P]
richness_SA_R<-amr_Subtype_Day$richness[index_10_SA_R]

wilcox_STP_STR_richness<-wilcox.test(richness_ST_P,richness_ST_R)
wilcox_STP_STR_richness   # p-value = 0.186 --> not significant
wilcox_STP_SAP_richness<-wilcox.test(richness_ST_P,richness_SA_P)
wilcox_STP_SAP_richness  #p-value = 5.794e-12 --> significant *******
wilcox_STP_SAR_richness<-wilcox.test(richness_ST_P,richness_SA_R)
wilcox_STP_SAR_richness  # p-value = 4.388e-12 --> significant *********
wilcox_STR_SAP_richness<-wilcox.test(richness_ST_R,richness_SA_P)
wilcox_STR_SAP_richness   # p-value = 1.741e-08 ---> significant    *********
wilcox_STR_SAR_richness<-wilcox.test(richness_ST_R,richness_SA_R)
wilcox_STR_SAR_richness    # p-value = 6.056e-09---> significant   ***********
wilcox_SAP_SAR_richness<-wilcox.test(richness_SA_P,richness_SA_R)
wilcox_SAP_SAR_richness    # p-value = 0.5489---> non significant

################## Shannon compare P_ST vs R_ST vs P_SA vs R_SA ###################
tail(colnames(amr_Subtype_Day))
View(amr_Subtype_Day)

shapiro_Shannon <-shapiro.test(amr_Subtype_Day$Shannon)
shapiro_Shannon # p-value = 0.0006582 ---> non - normal

kruskal.test(Shannon ~ Subtype , data=amr_Subtype_Day) #p-value = 7.61e-06 ---> significant
pairwise.wilcox.test(amr_Subtype_Day$Shannon, amr_Subtype_Day$Subtype)

#       P_SA    P_ST    R_SA   
#  P_ST 2.7e-06 -       -      
#  R_SA 0.31470 0.00012 -      
#  R_ST 0.10530 0.10530 0.19720``

################## Abundance(mean) compare P_ST vs R_ST vs P_SA vs R_SA ###################
tail(colnames(amr_Subtype_Day))
View(amr_Subtype_Day)

shapiro_Shannon_sum <-shapiro.test(amr_Subtype_Day$rowsum)
shapiro_Shannon_sum # p-value < 2.2e-16 ---> non - normal

kruskal.test(rowsum ~ Subtype , data=amr_Subtype_Day) # p-value < 2.2e-16 ---> significant
pairwise.wilcox.test(amr_Subtype_Day$rowsum, amr_Subtype_Day$Subtype)

#     P_SA    P_ST    R_SA   
#P_ST 3.1e-14 -       -      
#R_SA 1       < 2e-16 -      
#R_ST < 2e-16 1       < 2e-16
