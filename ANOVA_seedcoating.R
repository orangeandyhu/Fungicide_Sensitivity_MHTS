#########Seedcoating############


#EC50
rawdata <- import("./data/Seedcoating.csv", format="csv")
#先算residue
aov = aov(DIT ~ treatment, data =rawdata) 
par(mfrow=c(1,2))
plot(aov,which=(c(1,2)))
rawdata.stdres = rstandard(aov)

### Shapiro-Wilk normality test
shapiro.test(rawdata.stdres)
#Equal varience-Bartlett test
bartlett.test(DIT~treatment, data=rawdata)
HSD <- HSD.test(aov,"treatment",consol=TRUE)
TukeyHSD(aov)



#########Mortality###########
wdata<-import("./data/Mortality.csv", format="csv")
aov<-aov(mortality~treatment, data=wdata)
summary(aov)
par(mfrow=c(1,2))
plot(aov,which=(c(1,2)))
rawdata.stdres = rstandard(aov)
shapiro.test(rawdata.stdres)
#Equal varience-Bartlett test
bartlett.test(mortality~treatment, data=wdata)

#Tukey
HSD.test(aov,"treatment",consol=TRUE)
#LSD
LSD.test(aov,"treatment",consol=TRUE)

#kruskal-wallis
kruskal.test(mortality~treatment, data=wdata)
#dunn's
dunn<-dunnTest(mortality~treatment, data=wdata)
dunn.datafram<- as.data.frame(dunn$res)
cldList(P.adj~Comparison, data=dunn.datafram,
        threshold  = 0.05)

pairwise.wilcox.test(wdata$mortality,wdata$treatment, p.adjust.method="none")
pairwise.wilcox.test(wdata$mortality,wdata$treatment, p.adjust.method="bonferroni")
pairwise.wilcox.test(wdata$mortality,wdata$treatment, p.adjust.method="holm")