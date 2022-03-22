library(rio)
library(ggplot2)
CORraw<- import("./data/correlation.csv", format="csv")
cor.test(CORraw$X, CORraw$Y, method = "spearman")
############Fuopyram####################
FOD<-subset(CORraw,fun=="F"&trt=="OD")[,2]
FAM<-subset(CORraw,fun=="F"&trt=="AM")[,2]
FSP<-subset(CORraw,fun=="F"&trt=="SP")[,2]
#HTS vs amended
cor.test(FOD, FAM, method = "spearman")
#HTS vs spore
cor.test(FOD, FSP, method = "spearman")
#amended vs spore
cor.test(FAM, FSP, method = "spearman")


################Switch####################
SOD<-subset(CORraw,fun=="S"&trt=="OD")[,2]
SAM<-subset(CORraw,fun=="S"&trt=="AM")[,2]
SSP<-subset(CORraw,fun=="S"&trt=="SP")[,2]
#HTS vs amended
cor.test(SOD, SAM, method = "spearman")
#HTS vs spore
cor.test(SOD, SSP, method = "spearman")
#amended vs spore
cor.test(SAM, SSP, method = "spearman")

###############Tebuconazole#############
TOD<-subset(CORraw,fun=="T"&trt=="OD")[,2]
TAM<-subset(CORraw,fun=="T"&trt=="AM")[,2]
TSP<-subset(CORraw,fun=="T"&trt=="SP")[,2]
#HTS vs amended
cor.test(TOD, TAM, method = "spearman")
#HTS vs spore
cor.test(TOD, TSP, method = "spearman")
#amended vs spore
cor.test(TAM, TSP, method = "spearman")



