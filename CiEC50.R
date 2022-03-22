# Copyright (c) 2021 PinghuWu
# calculation of the EC50 for mircroplate test
# install.packages("drc")
library(drc)
library(rio)

########### microplate ############

##==Fluopyram==##
isolate_micro <- c("19F021","19F026","19F028","19F029","19F033","19F034","19F035","19F037","19F042","20F013","20F014","20F017","CS11_1a","CS11_10","CS11_11","CS11_12","CS11_13","CS11_14","CS11_15","CS11_16","CS11_8","CS12_1","CS13_1a","CS13_1_1","CS13_10","CS13_11","CS13_3","CS13_4","CS13_5","CS13_6","CS13_7","CS13_8","CS13_9","CS21_1a","CS21_10","CS21_11","CS21_2","CS21_3","CS21_4","CS21_5","CS21_6","CS21_7","CS21_8","CS21_9","F022","KSBJ10_1","KSBJ11_1","KSBJ13_1a","KSBJ13_1_1","KSBJ13_2","KSBJ13_3","KSBJ14_1","KSBJ14_2","KSBJ16_1","KSBJ17_1","KSBJ17_2","KSBJ22_1","KSBJ24_1","KSBJ25","KSBJ28","KSBJ29","KSBJ31a","KSBJ31_1","KSBJ33_1","KSBJ33_2","PTN3","PTN5","RT1","RT2","RT3","RT4","RT5","CS11_3","CS11_4","CS11_5","CS11_6","CS11_7","CS11_9","PTN1","F018")
#importfile
microAllraw<- import("./data/fun_F_all.csv", format="csv")

#Calculate relative growth rate by for-loop
micro_Tinhrall<- NULL
for (i in isolate_micro){
  
  for(j in 1:3){
    micro_T0raw <- subset(microAllraw,con==0&rep==j&iso==i)
    micro_Tavr <- sapply(micro_T0raw[6], mean)
    micro_Tavr<-c(1,1,micro_Tavr)
    micro_Tinh <- subset(microAllraw,rep==j&iso==i)
    micro_Tinh <- micro_Tinh[-c(2,4,5,7,8,9,10)]
    micro_Tinh <- data.matrix(micro_Tinh)
    micro_Tinh <- t(micro_Tinh)
    micro_Tinhr <- (micro_Tinh/micro_Tavr)
    B400 <- c(1,1,100)
    micro_Tinhr <- (micro_Tinhr*B400)
    micro_Tinhr <- t(micro_Tinhr)
    repiso<-rep(i,16)
    repiso<-data.matrix(repiso)
    micro_Tinhr<-cbind(repiso,micro_Tinhr)
    micro_Tinhrall <-rbind(micro_Tinhrall,micro_Tinhr)
  }
  
}
colnames(micro_Tinhrall)<-c("iso","rep","con","rgr")
write.table(micro_Tinhrall,file = "./data/fun_F_inh.csv")

#Calculate EC50
microTinhraw<- import("./data/fun_F_inh.csv", format="csv")
EC50.Tmic<-NULL
for(i in isolate_micro){
  micro_TEC50 <- subset(microTinhraw,iso==i)
  drm.temp <- drm(rgr ~ con, iso, data = micro_TEC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Tmic <- rbind(EC50.Tmic,EC50.temp)}

write.table(EC50.Tmic, file = "./result/EC50_F_LL4.csv")


##==Switch==##
isolate_micro <- c("19F021","19F026","19F028","19F029","19F033","19F035","19F034","19F037","19F042","20F013","20F014","20F017","CS11_1a","CS11_10","CS11_11","CS11_12","CS11_13","CS11_14","CS11_15","CS11_16","CS11_8","CS12_1","CS13_1a","CS13_1_1","CS13_10","CS13_11","CS13_3","CS13_4","CS13_5","CS13_6","CS13_7","CS13_8","CS13_9","CS21_1a","CS21_10","CS21_11","CS21_2","CS21_3","CS21_4","CS21_5","CS21_6","CS21_7","CS21_8","CS21_9","F022","KSBJ10_1","KSBJ11_1","KSBJ13_1a","KSBJ13_1_1","KSBJ13_2","KSBJ13_3","KSBJ14_1","KSBJ14_2","KSBJ16_1","KSBJ17_1","KSBJ17_2","KSBJ22_1","KSBJ24_1","KSBJ25","KSBJ28","KSBJ29","KSBJ31a","KSBJ31_1","KSBJ33_1","KSBJ33_2","PTN3","PTN5","RT1","RT2","RT3","RT4","RT5","CS11_3","CS11_4","CS11_5","CS11_6","CS11_7","CS11_9","PTN1","F018")

#importfile
microAllraw<- import("./data/fun_S_all.csv", format="csv")

#Calculate relative growth rate by for-loop
micro_Tinhrall<- NULL
for (i in isolate_micro){
  
  for(j in 1:3){
    micro_T0raw <- subset(microAllraw,con==0&rep==j&iso==i)
    micro_Tavr <- sapply(micro_T0raw[6], mean)
    micro_Tavr<-c(1,1,micro_Tavr)
    micro_Tinh <- subset(microAllraw,rep==j&iso==i)
    micro_Tinh <- micro_Tinh[-c(2,4,5,7,8,9,10)]
    micro_Tinh <- data.matrix(micro_Tinh)
    micro_Tinh <- t(micro_Tinh)
    micro_Tinhr <- (micro_Tinh/micro_Tavr)
    B400 <- c(1,1,100)
    micro_Tinhr <- (micro_Tinhr*B400)
    micro_Tinhr <- t(micro_Tinhr)
    repiso<-rep(i,16)
    repiso<-data.matrix(repiso)
    micro_Tinhr<-cbind(repiso,micro_Tinhr)
    micro_Tinhrall <-rbind(micro_Tinhrall,micro_Tinhr)
  }
  
}
colnames(micro_Tinhrall)<-c("iso","rep","con","rgr")
write.table(micro_Tinhrall,file = "./data/fun_S_inh.csv")

#Calculate EC50
microTinhraw<- import("./data/fun_S_inh.csv", format="csv")
EC50.Tmic<-NULL
for(i in isolate_micro){
  EC50.temp<-NULL
  micro_TEC50 <- subset(microTinhraw,iso==i)
  drm.temp <- try(drm(rgr ~ con, iso, data = micro_TEC50, fct = LL.4()))
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Tmic <- rbind(EC50.Tmic,EC50.temp)}

write.table(EC50.Tmic, file = "./result/EC50_S_LL.4.csv")

##==Tebuconazole==##
#define isolates
isolate_micro <- c("19F021","19F026","19F028","19F029","19F033","19F034","19F035","19F037","19F042","20F013a","20F014","20F017","CS11_1a","CS11_10","CS11_11","CS11_12","CS11_13","CS11_14","CS11_15","CS11_16","CS11_3","CS11_4","CS11_5","CS11_6","CS11_7","CS11_8","CS11_9","CS12_1","CS13_1a","CS13_1_1","CS13_10","CS13_11","CS13_3","CS13_4","CS13_5","CS13_6","CS13_7","CS13_8","CS13_9","CS21_1a","CS21_10","CS21_11","CS21_2","CS21_3","CS21_4","CS21_5","CS21_6","CS21_7","CS21_8","CS21_9","F022","KSBJ10_1","KSBJ11_1","KSBJ13_1a","KSBJ13_1_1","KSBJ13_2","KSBJ13_3","KSBJ14_1","KSBJ14_2","KSBJ16_1","KSBJ17_1","KSBJ17_2","KSBJ22_1","KSBJ24_1","KSBJ25","KSBJ28","KSBJ29","KSBJ31","KSBJ31_1","KSBJ33_1","KSBJ33_2","PTN1","PTN3","PTN4","PTN5","RT1","RT2","RT3","RT4","RT5","F018")
#importfile
microAllraw<- import("./data/fun_T_all.csv", format="csv")

#Calculate relative growth rate by for-loop
micro_Tinhrall<- NULL
for (i in isolate_micro){
  
  for(j in 1:3){
    micro_T0raw <- subset(microAllraw,con==0&rep==j&iso==i)
    micro_Tavr <- sapply(micro_T0raw[6], mean)
    micro_Tavr<-c(1,1,micro_Tavr)
    micro_Tinh <- subset(microAllraw,rep==j&iso==i)
    micro_Tinh <- micro_Tinh[-c(2,4,5,7,8,9,10)]
    micro_Tinh <- data.matrix(micro_Tinh)
    micro_Tinh <- t(micro_Tinh)
    micro_Tinhr <- (micro_Tinh/micro_Tavr)
    B400 <- c(1,1,100)
    micro_Tinhr <- (micro_Tinhr*B400)
    micro_Tinhr <- t(micro_Tinhr)
    repiso<-rep(i,16)
    repiso<-data.matrix(repiso)
    micro_Tinhr<-cbind(repiso,micro_Tinhr)
    micro_Tinhrall <-rbind(micro_Tinhrall,micro_Tinhr)
  }
  
}
colnames(micro_Tinhrall)<-c("iso","rep","con","rgr")
write.table(micro_Tinhrall,file = "./data/fun_T_inh.csv")

#Calculate EC50
microTinhraw<- import("./data/fun_T_inh.csv", format="csv")
EC50.Tmic<-NULL
for(i in isolate_micro){
  micro_TEC50 <- subset(microTinhraw,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = micro_TEC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Tmic <- rbind(EC50.Tmic,EC50.temp)}

write.table(EC50.Tmic, file = "./result/EC50_T_LL4.csv")


####################Amended#####################
##==Fluopyram==##
isolates<-c("19F042","CS11_10","CS11_1a","CS21_4","KSBJ11_1","KSBJ14_1","KSBJ17_1","KSBJ22","KSBJ24_1","KSBJ25")
Raw <- import("./data/amended_all.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="F")
    Tavr <- sapply(T0raw[6], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="F")
    Tinh <- Tinh[-c(3,4,5,7)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/am_F_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/am_F_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_am_F_LL4.csv")

##==Switch==##
isolates<-c("KSBJ11_1","CS12_1","CS13_1","CS21_1","F018","F022","KSBJ13_3","KSBJ29","CS13_8","PTN5")
Raw <- import("./data/amended_all.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="S")
        Tavr <- sapply(T0raw[6], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="S")
    Tinh <- Tinh[-c(3,4,5,7)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/am_S_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/am_S_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_am_S_LL4.csv")


##==Tebuconazole==##
isolates<-c("19F026","19F037","19F042","CS11_14","CS11_15","CS12_1","CS21_5","KSBJ10_1","KSBJ11_1","KSBJ31")
Raw <- import("./data/amended_all.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="T")
    Tavr <- sapply(T0raw[6], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="T")
    Tinh <- Tinh[-c(3,4,5,7)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/am_T_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/am_T_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_am_T_LL4.csv")

#################SPORE###################

##==Fluopyram==##
isolates<-c("19F042","CS11_10","CS11_1a","CS21_4","KSBJ11_1","KSBJ14_1","KSBJ17_1","KSBJ22","KSBJ24_1","KSBJ25")
Raw <- import("./data/spore.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="F")
    Tavr <- sapply(T0raw[4], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="F")
    Tinh <- Tinh[-c(3,5)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/sp_F_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/sp_F_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_sp_F_LL4.csv")

##==Switch==##
isolates<-c("KSBJ11_1","CS12_1","CS13_1","CS21_1","F018","F022","KSBJ13_3","KSBJ29","CS13_8","PTN5")

Raw <- import("./data/spore.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="S")
    Tavr <- sapply(T0raw[4], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="S")
    Tinh <- Tinh[-c(3,5)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/sp_S_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/sp_S_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_sp_S_LL4.csv")

##==Tebuconazole==##
isolates<-c("19F026","19F037","19F042","CS11_14","CS11_15","CS12_1","CS21_5","KSBJ10_1","KSBJ11_1","KSBJ31")

Raw <- import("./data/spore.csv", format="csv")
Tinhrall<- NULL
for (i in isolates){
  
  for(j in 1:3){
    T0raw <- subset(Raw, con==0&rep==j&iso==i &fun=="T")
    Tavr <- sapply(T0raw[4], mean)
    Tavr<-c(1,1,Tavr)
    Tinh <- subset(Raw, rep==j&iso==i&fun=="T")
    Tinh <- Tinh[-c(3,5)]
    Tinh <- data.matrix(Tinh)
    Tinh <- t(Tinh)
    Tinhr <- (Tinh/Tavr)
    S400 <- c(1,1,100)
    Tinhr <- (Tinhr*S400)
    Tinhr <- t(Tinhr)
    repiso<-rep(i,12)
    repiso<-data.matrix(repiso)
    Tinhr<-cbind(repiso,Tinhr)
    Tinhrall <-rbind(Tinhrall,Tinhr)
  }
}

colnames(Tinhrall)<-c("iso","rep","con","rgr")
write.table(Tinhrall,file = "./data/sp_T_inh.csv")

#Calculate EC50
Raw_inh<- import("./data/sp_T_inh.csv", format="csv")
EC50.Amic<-NULL
for(i in isolates){
  Raw_EC50 <- subset(Raw_inh,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = Raw_EC50, fct = LL.4())
  EC50.temp <- ED(drm.temp, c(50))
  EC50.Amic <- rbind(EC50.Amic,EC50.temp)}

write.table(EC50.Amic, file = "./result/EC50_sp_T_LL4.csv")

####################modelselection##############
mslect.Tmic<-NULL
for(i in isolate_micro){
  micro_TEC50 <- subset(microTinhraw,iso==i)
  drm.temp <- drm(rgr ~ con, iso ,data = micro_TEC50, fct = LL.４())
  modelselect <- mselect(drm.temp, list(LL.3(), LL.5(), CRS.4a(), W2.4(), BC.4()))
  repiso<-rep(i,6)
  repiso<-data.matrix(repiso)
  modelselect2<-cbind(repiso,modelselect)
  mslect.Tmic<-rbind(mslect.Tmic,modelselect2)
}
write.table(mslect.Tmic, file = "./result/mselect_S.csv")
