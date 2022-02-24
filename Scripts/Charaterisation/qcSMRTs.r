## perform SMRT cell level summary of data contribution


library(Rmisc)
library(ggplot2)
library(ggpubr)

sampleSheet <- read.csv("RawData/PacBio/"
smrtCellStats<-read.table("MasterTranscriptome/TALON/pfc_merge_smrt_talon_summary.tsv", stringsAsFactors = FALSE, header = TRUE)

smrtCellStats$genes<-rowSums(smrtCellStats[,c("known_genes", "antisense_genes", "other_novel_genes")])
smrtCellStats$transcripts<-rowSums(smrtCellStats[,c("known_transcripts", "novel_transcripts")])
smrtCellStats$novel_genes<-rowSums(smrtCellStats[,c("antisense_genes", "other_novel_genes")])


## summarise data from each smrt cell

figa1<-ggplot(smrtCellStats, aes(x = "", y=reads_annotated)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
    ylab("Reads annotated") +
    xlab("") +
    theme_bw() 
	
figa2<-ggplot(smrtCellStats, aes(x = "", y=genes)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
    ylab("Number of genes")+
    xlab("") +
    theme_bw() 

figa3<-ggplot(smrtCellStats, aes(x = "", y=other_novel_genes)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
    ylab("Number of novel genes")+
    xlab("") +
    theme_bw() 
	
figa4<-ggplot(smrtCellStats, aes(x = "", y=transcripts)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
    ylab("Number of transcripts")+
    xlab("") +
    theme_bw() 

figa5<-ggplot(smrtCellStats, aes(x = "", y=novel_transcripts)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
    ylab("Number of novel transcripts")+
    xlab("") +
    theme_bw() 


  
figa<-ggarrange(figa1, figa2, figa3,figa4,figa5,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 5, nrow = 1)

pdf("RawData/PacBio/QCOutput/ViolinPlotNumberofFeaturesBySMRTcell.pdf", width = 15, height = 7)
print(figa)
dev.off()



## concerned about samples that have more novel genes/transcripts than expected

figb1<-ggplot(smrtCellStats, aes(x=known_genes, y=novel_genes)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Known genes") +
    ylab("Novel & Antisense genes") +
    theme_bw() 
	
figb2<-ggplot(smrtCellStats, aes(x=known_transcripts, y=novel_transcripts)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Known transcripts") +
    ylab("Novel transcripts") +
    theme_bw() 

figb<-ggarrange(figb1, figb2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

pdf("RawData/PacBio/QCOutput/ScatterplotNumberofFeaturesBySMRTcell.pdf", width = 6, height = 7)
print(figb)
dev.off()

figc1<-ggplot(smrtCellStats, aes(x=novel_transcripts, y=ISMs)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Novel transcripts") +
    ylab("Incomplete splice match") +
    theme_bw() 
	
figc2<-ggplot(smrtCellStats, aes(x=novel_transcripts, y=prefix_ISMs)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Novel transcripts") +
    ylab("Prefix incomplete splice match") +
    theme_bw() 

figc3<-ggplot(smrtCellStats, aes(x=novel_transcripts, y=suffix_ISMs)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Novel transcripts") +
    ylab("Suffix incomplete splice match") +
    theme_bw() 
	
figc4<-ggplot(smrtCellStats, aes(x=novel_transcripts, y=NICs)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Novel transcripts") +
    ylab("Novel in catelogue") +
    theme_bw() 
	
figc5<-ggplot(smrtCellStats, aes(x=novel_transcripts, y=NNCs)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm', se = TRUE, fill="#69b3a2") +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Novel transcripts") +
    ylab("Novel not in catelogue") +
    theme_bw() 
	
figc<-ggarrange(figc1, figc2, figc3,figc4,figc5,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

pdf("RawData/PacBio/QCOutput/ScatterplotNumberofNovelTranscriptsBySMRTcell.pdf", width = 15, height = 10)
print(figc)
dev.off()