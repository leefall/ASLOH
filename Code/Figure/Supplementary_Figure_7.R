library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)


TCGA.Result.Data<-data.frame(read.table("../../Data/TNBC_Survival_Result.txt",stringsAsFactors = FALSE,header = TRUE))



colnames(TCGA.Result.Data)<-c("Gene_LOH","Pvalue","HR","Qvalue")

FilterTCGA.Result.Data<-TCGA.Result.Data[which(TCGA.Result.Data$N>3),]


ggplot(FilterTCGA.Result.Data, aes(x=HR,y=-log10(Pvalue),label=Gene_LOH))+
  geom_point(size=2,color=dplyr::case_when(FilterTCGA.Result.Data$Qvalue<0.05~"red",FilterTCGA.Result.Data$Qvalue>0.05 ~"#CCCCCC"))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",color="red")+
  geom_vline(xintercept = 1,linetype="dashed",color="black")+
  
  geom_label_repel(aes(label=ifelse((Qvalue<0.05),as.character(Gene_LOH),'')),box.padding = unit(5, 'lines'),
                   point.padding=unit(0.2,"lines"),segment.color='grey50',max.overlaps = 1000,max.iter=1e4)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=15),legend.position = "none")+
  labs(y="-log10(P-value)",x="Hazard Ratio (HR)")






