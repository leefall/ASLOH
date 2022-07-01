library(ggplot2)
library(ggsignif)

write.table(Total.Call.Algorithm,file="../../Data/BenchmarkResultwithoutFilter.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)

Total.Call.Algorithm<-data.frame(read.table("../../Data/BenchmarkResultwithoutFilter.txt",header = TRUE,
                                      stringsAsFactors = FALSE))



#ASLOH
no.ASLOH.mean<-mean(Total.Call.Algorithm$ASLOH,na.rm = TRUE)
no.ASLOH.sd<-sd(Total.Call.Algorithm$ASLOH,na.rm = TRUE)
#Varscan2
no.Varscan2.mean<-mean(Total.Call.Algorithm$Varscan2,na.rm = TRUE)
no.Varscan2.sd<-sd(Total.Call.Algorithm$Varscan2,na.rm = TRUE)
#SomaticSniper
no.SomaticSniper.mean<-mean(Total.Call.Algorithm$SomaticSniper,na.rm = TRUE)
no.SomaticSniper.sd<-sd(Total.Call.Algorithm$SomaticSniper,na.rm = TRUE)
#SENVsniffer
no.SNVsniffer.mean<-mean(Total.Call.Algorithm$SNVsniffer,na.rm = TRUE)
no.SNVsniffer.sd<-sd(Total.Call.Algorithm$SNVsniffer,na.rm = TRUE)
#Sequenza
no.Sequenza.mean<-mean(Total.Call.Algorithm$Sequenza,na.rm = TRUE)
no.Sequenza.sd<-sd(Total.Call.Algorithm$Sequenza,na.rm = TRUE)


#t.test(Total.Call.Algorithm$ASLOH, Total.Call.Algorithm$Varscan2, paired = TRUE, alternative = "two.sided")
no.t.Varscan2<-t.test(Total.Call.Algorithm$ASLOH, Total.Call.Algorithm$Varscan2, paired = FALSE, alternative = "two.sided")
no.t.SomaticSniper<-t.test(Total.Call.Algorithm$ASLOH, Total.Call.Algorithm$SomaticSniper, paired = FALSE, alternative = "two.sided")
no.t.SNVsniffer<-t.test(Total.Call.Algorithm$ASLOH, Total.Call.Algorithm$SNVsniffer, paired = FALSE, alternative = "two.sided")
no.t.Sequenza<-t.test(Total.Call.Algorithm$ASLOH, Total.Call.Algorithm$Sequenza, paired = FALSE, alternative = "two.sided")

#Deleterious



Total.Del.Algorithm<-data.frame(read.table("../../Data/BenchmarkResultwithoutFilter.txt",header = TRUE,
                                      stringsAsFactors = FALSE))






#ASLOH
Del.ASLOH.mean<-mean(Total.Del.Algorithm$ASLOH,na.rm = TRUE)
Del.ASLOH.sd<-sd(Total.Del.Algorithm$ASLOH,na.rm = TRUE)
#Varscan2
Del.Varscan2.mean<-mean(Total.Del.Algorithm$Varscan2,na.rm = TRUE)
Del.Varscan2.sd<-sd(Total.Del.Algorithm$Varscan2,na.rm = TRUE)
#SomaticSniper
Del.SomaticSniper.mean<-mean(Total.Del.Algorithm$SomaticSniper,na.rm = TRUE)
Del.SomaticSniper.sd<-sd(Total.Del.Algorithm$SomaticSniper,na.rm = TRUE)
#SENVsniffer
Del.SNVsniffer.mean<-mean(Total.Del.Algorithm$SNVsniffer,na.rm = TRUE)
Del.SNVsniffer.sd<-sd(Total.Del.Algorithm$SNVsniffer,na.rm = TRUE)
#Sequenza
Del.Sequenza.mean<-mean(Total.Del.Algorithm$Sequenza,na.rm = TRUE)
Del.Sequenza.sd<-sd(Total.Del.Algorithm$Sequenza,na.rm = TRUE)


Del.t.Varscan2<-t.test(Total.Del.Call$ASLOH, Total.Del.Call$Varscan2, paired = FALSE, alternative = "two.sided")
Del.t.SomaticSniper<-t.test(Total.Del.Call$ASLOH, Total.Del.Call$SomaticSniper, paired = FALSE, alternative = "two.sided")
Del.t.SNVsniffer<-t.test(Total.Del.Call$ASLOH, Total.Del.Call$SNVsniffer, paired = FALSE, alternative = "two.sided")
Del.t.Sequenza<-t.test(Total.Del.Call$ASLOH, Total.Del.Call$Sequenza, paired = FALSE, alternative = "two.sided")




#GTEX Filtering

Total.GTEXCall.Algorithm<-data.frame(read.table("../../Data/BenchmarkResultinSilicoeQTL.txt",header = TRUE,
                                           stringsAsFactors = FALSE))



GTEX.t.Varscan2<-t.test(Total.GTEX.Call$ASLOH, Total.GTEX.Call$Varscan2, paired = FALSE, alternative = "two.sided")
GTEX.t.SomaticSniper<-t.test(Total.GTEX.Call$ASLOH, Total.GTEX.Call$SomaticSniper, paired = FALSE, alternative = "two.sided")
GTEX.t.SNVsniffer<-t.test(Total.GTEX.Call$ASLOH, Total.GTEX.Call$SNVsniffer, paired = FALSE, alternative = "two.sided")
GTEX.t.Sequenza<-t.test(Total.GTEX.Call$ASLOH, Total.GTEX.Call$Sequenza, paired = FALSE, alternative = "two.sided")



#ASLOH
GTEX.ASLOH.mean<-mean(Total.GTEXCall.Algorithm$ASLOH,na.rm = TRUE)
GTEX.ASLOH.sd<-sd(Total.GTEXCall.Algorithm$ASLOH,na.rm = TRUE)
#Varscan2
GTEX.Varscan2.mean<-mean(Total.GTEXCall.Algorithm$Varscan2,na.rm = TRUE)
GTEX.Varscan2.sd<-sd(Total.GTEXCall.Algorithm$Varscan2,na.rm = TRUE)
#SomaticSniper
GTEX.SomaticSniper.mean<-mean(Total.GTEXCall.Algorithm$SomaticSniper,na.rm = TRUE)
GTEX.SomaticSniper.sd<-sd(Total.GTEXCall.Algorithm$SomaticSniper,na.rm = TRUE)
#SENVsniffer
GTEX.SNVsniffer.mean<-mean(Total.GTEXCall.Algorithm$SNVsniffer,na.rm = TRUE)
GTEX.SNVsniffer.sd<-sd(Total.GTEXCall.Algorithm$SNVsniffer,na.rm = TRUE)
#Sequenza
GTEX.Sequenza.mean<-mean(Total.GTEXCall.Algorithm$Sequenza,na.rm = TRUE)
GTEX.Sequenza.sd<-sd(Total.GTEXCall.Algorithm$Sequenza,na.rm = TRUE)



mydata <- data.frame(Filtering=c("Non-Filtered","Non-Filtered","Non-Filtered","Non-Filtered","Non-Filtered",
                                 "In silico","In silico","In silico","In silico","In silico",
                                 "eQTL","eQTL","eQTL","eQTL","eQTL"),
                     SD = c(no.ASLOH.sd, no.Varscan2.sd, no.SomaticSniper.sd, no.SNVsniffer.sd, no.Sequenza.sd,
                            Del.ASLOH.sd, Del.Varscan2.sd, Del.SomaticSniper.sd, Del.SNVsniffer.sd, Del.Sequenza.sd,
                            GTEX.ASLOH.sd, GTEX.Varscan2.sd, GTEX.SomaticSniper.sd, GTEX.SNVsniffer.sd, GTEX.Sequenza.sd), 
                     Algorithm = c("ASLOH", "Varscan2", "SomaticSniper", "SNVsniffer", "Sequenza",
                                    "ASLOH", "Varscan2", "SomaticSniper", "SNVsniffer", "Sequenza",
                                    "ASLOH", "Varscan2", "SomaticSniper", "SNVsniffer", "Sequenza"), 
                     Mean = c(no.ASLOH.mean,no.Varscan2.mean,no.SomaticSniper.mean,no.SNVsniffer.mean, no.Sequenza.mean,
                              Del.ASLOH.mean,Del.Varscan2.mean,Del.SomaticSniper.mean,Del.SNVsniffer.mean, Del.Sequenza.mean,
                              GTEX.ASLOH.mean,GTEX.Varscan2.mean,GTEX.SomaticSniper.mean,GTEX.SNVsniffer.mean, GTEX.Sequenza.mean))


mydata$Algorithm<-factor(mydata$Algorithm,levels = c("ASLOH","Varscan2","SomaticSniper","SNVsniffer","Sequenza"))
mydata$Filtering<-factor(mydata$Filtering,levels = c("Non-Filtered","In silico","eQTL"))

ggplot(mydata, aes(x=Filtering, y=Mean, fill=Algorithm)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9)) +ylab("Average Concordance rates (%)")+ 
  scale_y_continuous(breaks=c(0,25,50,75,100),limits = c(-6.52,115))+
  geom_signif(stat="identity",
              data=data.frame(x=c(0.65, 0.65,0.65,0.65,
                                  1.65,1.65,1.65,1.65,
                                  2.65,2.65,2.65,2.65), 
                              xend=c(0.82, 1,1.18,1.36,
                                     1.82,2,2.18,2.36,
                                     2.82,3,3.18,3.36),
                              y=c(114, 113,112,111,
                                  114,113,112,111,
                                  114,113,112,107), 
                              annotation=c("**", "**","**","**",
                                           "** ","** ","** ","** ",
                                           "**  ","**  ","**  ","*"),
                              Algorithm=c("ASLOH","ASLOH","ASLOH","ASLOH",
                                           "ASLOH","ASLOH","ASLOH","ASLOH",
                                           "ASLOH","ASLOH","ASLOH","ASLOH")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))+theme(axis.line = element_line(colour = "black"),
                                                                            panel.grid.major = element_blank(),
                                                                            panel.grid.minor = element_blank(),
                                                                            panel.border = element_blank(),
                                                                            panel.background = element_blank(),
                                                                            text = element_text(size=15))




  
