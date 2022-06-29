library(grid)
library(gridExtra)
library(survival)
library(ggplot2)
library("survminer")

TotalData<-data.frame(read.table("../../Data/KIF15_BRCAKaplerMeier.txt",header = TRUE,
                                 stringsAsFactors = FALSE))


Clinical.TotalData<-TotalData
Clinical.TotalData$GII<-factor(Clinical.TotalData$GII,levels=c("Low","High"))
Clinical.TotalData$TMB<-factor(Clinical.TotalData$TMB,levels=c("Low","High"))

sGene<-"KIF15"
sTitle<-"KIF15 (BRCA)"


Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="RLOH"),]$Group<-"rLOH"}




Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))
logtest <- survdiff(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=7)
Kapler.meier.data <- survfit(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OSdays, OS=="1")~Group+TumorStage+Age+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=4)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=6)

TCGA.TNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata", 
                     title=sTitle,legend.title="",xlab="Time in days",pval=paste(
                       "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                       "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,ncensor.plot.height = 0.15,risk.table.height = 0.30,
                     pval.coord = c(3500, 0.1))


TCGA.TNBCsurvp.SS$table<- ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
            y.text=F,ylab="",lab="",tables.theme = theme_cleantable())





TotalData<-data.frame(read.table("../../Data/KIF15_mBRCAKaplerMeier.txt",header = TRUE,
                                 stringsAsFactors = FALSE))






Clinical.TotalData<-TotalData
Clinical.TotalData$GII<-factor(Clinical.TotalData$GII,levels=c("Low","High"))
Clinical.TotalData$TMB<-factor(Clinical.TotalData$TMB,levels=c("Low","High"))



Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="RLOH"),]$Group<-"rLOH"}



Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))


logtest <- survdiff(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=4)
Kapler.meier.data <- survfit(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OS_Days, Expired=="1")~Group+HistologicGrade+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=4)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=4)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=4)

survp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata", 
                     title=sTitle,legend.title="",xlab="Time in days",pval=paste(
                       "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                       "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,ncensor.plot.height = 0.15,risk.table.height = 0.30,
                     pval.coord = c(800, 0.75))



survp.SS$table<- ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                             y.text=F,ylab="",lab="",tables.theme = theme_cleantable())



ssplots<-list()
ssplots[[1]]<-TCGA.TNBCsurvp.SS
ssplots[[2]]<-survp.SS


ss<-arrange_ggsurvplots(ssplots, print = TRUE,
                        ncol = 2, nrow = 1, risk.table.height = 0.2)

