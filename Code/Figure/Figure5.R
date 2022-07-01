library(grid)
library(gridExtra)
library(survival)
library(ggplot2)
library("survminer")


mTNBCData<-data.frame(read.table("../../Data/mTNBC_KaplerMeier.txt",header = TRUE,stringsAsFactors = FALSE))



Clinical.TotalData<-mTNBCData


Clinical.TotalData$GII<-factor(Clinical.TotalData$GII,levels=c("Low","High"))
Clinical.TotalData$TMB<-factor(Clinical.TotalData$TMB,levels=c("Low","High"))
Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="LOH"),]$Group<-"fLOH"}

if ("RLOH" %in% Clinical.TotalData$v3_44884647_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v3_44884647_C_A=="RLOH"),]$Group<-"rLOH"}


sGene<-"KIF15"
sTitle<-"KIF15 in mTNBC"




logtest <- survdiff(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)

summary(coxph(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData))

coxsummary<-summary(coxph(Surv(OS_Days, Expired=="1")~Group+HistologicGrade+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=3)
KIF15.mTNBC.survp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                     legend.labs=c("Non-LOH","fLOH","rLOH"),
                     title=sTitle,xlab="Time in days",pval=paste(
                       "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                       "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(400, 0.75))


KIF15.mTNBC.survp.SS$table<-ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
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
sTitle<-"MYH13 in mTNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","fLOH","rLOH"))

logtest <- survdiff(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)

summary(coxph(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData))


coxsummary<-summary(coxph(Surv(OS_Days, Expired=="1")~Group+HistologicGrade+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupfLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupfLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GrouprLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GrouprLOH"],digit=3)
MYH13.mTNBCsurvp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                          legend.labs=c("Non-LOH","fLOH","rLOH"),
                          title=sTitle,xlab="Time in days",pval=paste(
                            "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                            "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(550, 0.75))

MYH13.mTNBCsurvp.SS$table<-ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                 y.text=F,ylab="",lab="",tables.theme = theme_cleantable())


Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v9_130536717_G_A){
  Clinical.TotalData[which(Clinical.TotalData$v9_130536717_G_A=="LOH"),]$Group<-"FLOH"}

if ("RLOH" %in% Clinical.TotalData$v9_130536717_G_A){
  Clinical.TotalData[which(Clinical.TotalData$v9_130536717_G_A=="RLOH"),]$Group<-"RLOH"}

sGene<-"SH2D3C"
sTitle<-"SH2D3C in mTNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","FLOH","RLOH"))

logtest <- survdiff(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)

coxsummary<-summary(coxph(Surv(OS_Days, Expired=="1")~Group+HistologicGrade+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupFLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupFLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GroupRLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GroupRLOH"],digit=3)
SH2D3C.mTNBC.survp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                                 legend.labs=c("Non-LOH","fLOH","rLOH"),
                                 title=sTitle,xlab="Time in days",pval=paste(
                                   "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                                   "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(400, 0.75))


SH2D3C.mTNBC.survp.SS$table<-ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                        y.text=F,ylab="",lab="",tables.theme = theme_cleantable())



Clinical.TotalData$Group<-"NonLOH"

if ("LOH" %in% Clinical.TotalData$v15_55722882_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722882_C_A=="LOH"),]$Group<-"FLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55722882_C_A){
  Clinical.TotalData[which(Clinical.TotalData$v15_55722882_C_A=="RLOH"),]$Group<-"RLOH"}

if ("LOH" %in% Clinical.TotalData$v15_55789910_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55789910_C_T=="LOH"),]$Group<-"FLOH"}

if ("RLOH" %in% Clinical.TotalData$v15_55789910_C_T){
  Clinical.TotalData[which(Clinical.TotalData$v15_55789910_C_T=="RLOH"),]$Group<-"RLOH"}


sGene<-"DYX1C1"
sTitle<-"DYX1C1 in mTNBC"

Clinical.TotalData$Group<-factor(Clinical.TotalData$Group,levels=c("NonLOH","FLOH","RLOH"))

logtest <- survdiff(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)
real.log.p<-1-pchisq(logtest$chisq, length(unique(Clinical.TotalData$Group))-1)
real.log.p<-round(real.log.p,digit=3)
Kapler.meier.data <- survfit(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData)

summary(coxph(Surv(OS_Days, Expired=="1")~Group,data=Clinical.TotalData))


coxsummary<-summary(coxph(Surv(OS_Days, Expired=="1")~Group+HistologicGrade+GII+TMB,data=Clinical.TotalData))
Hetero.H.Cox<-round(coxsummary$coefficients[,2]["GroupFLOH"],digit=3)
Hetero.P.Cox<-round(coxsummary$coefficients[,5]["GroupFLOH"],digit=3)
HomoAlt.H.Cox<-round(coxsummary$coefficients[,2]["GroupRLOH"],digit=3)
HomoAlt.P.Cox<-round(coxsummary$coefficients[,5]["GroupRLOH"],digit=3)
DYX1C1.mTNBC.survp.SS<-ggsurvplot(Kapler.meier.data,data=Clinical.TotalData,size=1,risk.table=TRUE,risk.table.col="strata",
                                  legend.labs=c("Non-LOH","fLOH","rLOH"),
                                  title=sTitle,xlab="Time in days",pval=paste(
                                    "fLOH vs Non-LOH Cox P=",Hetero.P.Cox,
                                    "
rLOH vs Non-LOH Cox P=",HomoAlt.P.Cox,"
Log-rank P=",real.log.p,sep=""),ncensor.plot = FALSE,legend.title="",ncensor.plot.height = 0.15,risk.table.height = 0.30, pval.coord = c(10, 0.15))


DYX1C1.mTNBC.survp.SS$table<-ggrisktable(Kapler.meier.data,data=Clinical.TotalData,color="strata",
                                         y.text=F,ylab="",lab="",tables.theme = theme_cleantable())






ssplots<-list()
ssplots[[1]]<-KIF15.mTNBC.survp.SS
ssplots[[3]]<-MYH13.mTNBCsurvp.SS
ssplots[[2]]<-SH2D3C.mTNBC.survp.SS
ssplots[[4]]<-DYX1C1.mTNBC.survp.SS

arrange_ggsurvplots(ssplots, print = TRUE,
                        ncol = 2, nrow = 2, risk.table.height = 0.2)



