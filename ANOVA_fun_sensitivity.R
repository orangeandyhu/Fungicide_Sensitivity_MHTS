library(agricolae)
library(rio)
library(lawstat)
library(FSA)
library(rcompanion)


###########Inisial selection############
rawdata <- import("./data/Fun_sensitivity.csv", format="csv")
#ANOVA
aov <- aov(ori ~ treatment, data=rawdata)
summary(aov)
#Q-Q and vareince plot
par(mfrow=c(1,2))
plot(aov,which=(c(1,2)))
#算residue
rawdata.stdres = rstandard(aov)
### Shapiro-Wilk normality test
shapiro.test(rawdata.stdres)
#Equal varience-Bartlett test
bartlett.test(ori~treatment, data=rawdata)
#posthoc
HSD <- HSD.test(aov,"treatment",consol=TRUE)
boxplot(ori~treatment,data=rawdata)

#kruskal-wallis
kruskal.test(ori ~ treatment, data=rawdata)
#dunn's
dunn<-dunnTest(ori ~ treatment, data=rawdata)
dunn.datafram<- as.data.frame(dunn$res)
cldList(P.adj~Comparison, data=dunn.datafram,
        threshold  = 0.05)