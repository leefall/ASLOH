library(grid)
library(gridExtra)
library(survival)
library(ggplot2)
library("survminer")



TotalData<-data.frame(read.table("../../Data/TNBC_KaplerMeier.txt",header = TRUE,
                                 stringsAsFactors = FALSE))



Clinical.TotalData<-TotalData

Clinical.TotalData$GIIGroup<-factor(Clinical.TotalData$GII,levels=c("Low","High"))
Clinical.TotalData$TMB<-factor(Clinical.TotalData$TMB,levels=c("Low","High"))


Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="RLOH"),]$Group<-"rLOH"}


sGene<-"KIF15"
sTitle<-"KIF15 in TNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))


logtest <- survdiff(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OSdays, OS=="1")~Group+TumorStage+Age+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=4)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=4)
KIF15.TNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                              legend.labs=c("Non-LOH","fLOH","rLOH"),
                              title=sTitle,xlab="Time in days",pval=paste(
                                "fLOH vs Non-LOH Cox P=",round(Hetero.P.Cox,3),
                                "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(4000,0.25))



KIF15.TNBCsurvp.SS$table<-ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                     y.text=F,ylab="",lab="",tables.theme = theme_cleantable())





Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v17_10219113_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v17_10219113_C_T=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v17_10219113_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v17_10219113_C_T=="RLOH"),]$Group<-"rLOH"}

if ("LOH" %in% Clinical.TotalData$v17_10223714_T_C){
  Clinical.TotalData[which(Clinical.TotalData$v17_10223714_T_C=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v17_10223714_T_C){
  Clinical.TotalData[which(Clinical.TotalData$v17_10223714_T_C=="RLOH"),]$Group<-"rLOH"}


sGene<-"MYH13"
sTitle<-"MYH13 in TNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))



logtest <- survdiff(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OSdays, OS=="1")~Group+TumorStage+Age+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=3)

MYH13.TNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                              legend.labs=c("Non-LOH","fLOH","rLOH"),
                              title=sTitle,xlab="Time in days",pval=paste(
                                "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                                "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(4000, 0.25))




MYH13.TNBCsurvp.SS$table<- ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                      y.text=F,ylab="",lab="",tables.theme = theme_cleantable())






Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v9_130536717_G_A){
  Clinical.TotalData[which(Clinical.TotalData$v9_130536717_G_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v9_130536717_G_A){
  Clinical.TotalData[which(Clinical.TotalData$v9_130536717_G_A=="RLOH"),]$Group<-"rLOH"}

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))
sGene<-"SH2D3C"
sTitle<-"SH2D3C in TNBC"


#Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("HomoReference","GermlineHetero","GermlineHomo","Somatic","RLOH","LOH"))
logtest <- survdiff(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OSdays, OS=="1")~Group+TumorStage+Age+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=3)


SH2D3C.TNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                               legend.labs=c("Non-LOH","fLOH","rLOH"),
                               title=sTitle,xlab="Time in days",pval=paste(
                                 "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                                 "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(4000, 0.25))




SH2D3C.TNBCsurvp.SS$table<- ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                       y.text=F,ylab="",lab="",tables.theme = theme_cleantable())





Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v15_55790414_G_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55790414_G_T=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55790414_G_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55790414_G_T=="RLOH"),]$Group<-"rLOH"}

if ("LOH" %in% Clinical.TotalData$v15_55722882_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722882_C_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55722882_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722882_C_A=="RLOH"),]$Group<-"rLOH"}

if ("LOH" %in% Clinical.TotalData$v15_55789910_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55789910_C_T=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55789910_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55789910_C_T=="RLOH"),]$Group<-"rLOH"}

if ("LOH" %in% Clinical.TotalData$v15_55722872_G_C){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722872_G_C=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55722872_G_C){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722872_G_C=="RLOH"),]$Group<-"rLOH"}

sGene<-"DYX1C1"
sTitle<-"DYX1C1 in TNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))


logtest <- survdiff(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OSdays, OS=="1")~Group,data=Clinical.TotalData)
coxsummary<-summary(coxph(Surv(OSdays, OS=="1")~Group+TumorStage+Age+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=4)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=4)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=3)



DYX1C1.TNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                                legend.labs=c("Non-LOH","fLOH","rLOH"),
                                title=sTitle,xlab="Time in days",pval=paste(
                                  "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                                  "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(4700, 0.25))




DYX1C1.TNBCsurvp.SS$table<- ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                        y.text=F,ylab="",lab="",tables.theme = theme_cleantable())



ssplots<-list()
ssplots[[1]]<-KIF15.TNBCsurvp.SS
ssplots[[3]]<-MYH13.TNBCsurvp.SS
ssplots[[2]]<-SH2D3C.TNBCsurvp.SS
ssplots[[4]]<-DYX1C1.TNBCsurvp.SS

arrange_ggsurvplots(ssplots, print = TRUE,
                        ncol = 2, nrow = 2, risk.table.height = 0.2)

