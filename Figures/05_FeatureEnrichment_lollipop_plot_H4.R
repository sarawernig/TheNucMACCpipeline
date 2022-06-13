#!/usr/bin/env Rscript

#LogFC lollipop plot
setwd("path/to/plots")

library(dplyr)
library(ggplot2)

input <- read.csv2("FeatureDistribution_LogFC_enrichment_H4.csv", encoding="UTF-8",sep="\t")

############################################################################
########### Hypogeometric test of enrichment/depletion #####################
############################################################################

#phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

overlap.group1 <- as.numeric(input$Overlap)
#no.overlap.group1 <- as.numeric(input$NotOverlap) #for Fishers only
group1 <- as.numeric(input$Total)
group2 <- as.numeric(input$totalNucs)
overlap.group2 <- as.numeric(input$Overlap_subNucs)
  
probabilities <- phyper(overlap.group1-1, group1, group2, overlap.group2, lower.tail=TRUE)
input$dhyper <- as.numeric(probabilities)

#############################################################################
######################## Colored by significance ############################
#############################################################################

plot1 <- ggplot(input[1:6,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color= as.numeric(-log10(dhyper)))) +
  geom_point(aes(size = as.numeric(Size), color= as.numeric(-log10(dhyper))), alpha = 0.5) +
  scale_color_gradient2(low="darkblue",mid="purple",high="darkred",midpoint=1.3) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.y=element_text(angle=50, size=12, vjust=0.5,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
    ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,15)) +
  coord_cartesian(xlim = c(-0.6,0.4))

plot2 <- ggplot(input[7:12,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color= as.numeric(-log10(dhyper)))) +
  geom_point(aes(size = as.numeric(Size), color= as.numeric(-log10(dhyper))), alpha = 0.5) +
  scale_color_gradient2(low="darkblue",mid="purple",high="darkred",midpoint=1.3) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.y=element_text(angle=50, size=12, vjust=0.5,color="black"),
    axis.text.x=element_text(size=12, color="black")) +
    #plot.margin = unit(c(3, 3, 4, 3), "cm")) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,15)) +
  coord_cartesian(xlim = c(-0.6,0.4))


plot3 <- ggplot(input[13:18,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color= as.numeric(-log10(dhyper)))) +
  geom_point(aes(size = as.numeric(Size), color= as.numeric(-log10(dhyper))), alpha = 0.5) +
  scale_color_gradient2(low="darkblue",mid="purple",high="darkred",midpoint=1.3) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.y=element_text(angle=50, size=12, vjust=0.5,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,15)) +
  coord_cartesian(xlim = c(-0.6,0.4))


plot4 <- ggplot(input[19:24,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color= as.numeric(-log10(dhyper)))) +
  geom_point(aes(size = as.numeric(Size), color= as.numeric(-log10(dhyper))), alpha = 0.5) +
  scale_color_gradient2(low="darkblue",mid="purple",high="darkred",midpoint=1.3) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.y=element_text(angle=50, size=12, vjust=0.5,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,15)) +
  coord_cartesian(xlim = c(-0.6,0.4))


pdf("lollipop_plot_p-values.pdf")
par(mfrow=c(2,2))
plot1
plot2
plot3 
plot4
dev.off()

###########################################
########### Colored by feature ############
###########################################

plot1 <- ggplot(input[1:6,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color = Feature)) +
  geom_point(aes(size = as.numeric(Size), color = Feature), alpha = 0.5) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,20)) +
  coord_cartesian(xlim = c(-0.6,0.4))

plot2 <- ggplot(input[7:12,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color = Feature)) +
  geom_point(aes(size = as.numeric(Size), color = Feature), alpha = 0.5) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,20)) +
  coord_cartesian(xlim = c(-0.6,0.4))


plot3 <- ggplot(input[13:18,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color = Feature)) +
  geom_point(aes(size = as.numeric(Size), color = Feature), alpha = 0.5) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,20)) +
  coord_cartesian(xlim = c(-0.6,0.4))


plot4 <- ggplot(input[19:24,], aes(y=Feature, x=as.numeric(logFC)),size=as.numeric(Size)) +
  geom_segment( aes(y= Feature, yend= Feature, x=0, xend=as.numeric(logFC), color = Feature)) +
  geom_point(aes(size = as.numeric(Size), color = Feature), alpha = 0.5) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
    #plot.margin = unit(c(3, 3, 4, 3), "cm")      #top, right, bottom, left) 
  ) +  
  ylab("Genomic feature") +
  xlab("Fold enrichment") +
  geom_vline(aes(xintercept = 0), color = "black", size = 1) +
  scale_size(range = c(5,20)) +
  coord_cartesian(xlim = c(-0.6,0.4))


pdf("lollipop_plot.pdf")
par(mfrow=c(2,2))
plot1
plot2
plot3 
plot4
dev.off()

write.table(input,file="FeatureDistribution_LogFC_enrichment_phyper_H4.csv", quote=FALSE,row.names = FALSE,sep = "\t")


#############################################################################
########################### FINAL PLOTS #####################################
#############################################################################

input.list <- split(input, f = input$Category)

#convert to numeric
for(i in names(input.list)){
    input.list[[i]]$logFC<-as.numeric(as.character(input.list[[i]]$logFC))
}

#get maximum values
max.X<-max(round(sapply(input.list, function(x) max(abs(x$logFC))),1))

    
for(i in names(input.list)){
     g<-ggplot(input.list[[i]], aes(y = logFC, x=Feature ,size=as.numeric(Size))   )+
         geom_hline(yintercept = 0, color = "black", size = 1)+
        geom_bar(aes(Feature, logFC, fill=Feature), stat="identity", width = 0.05 )+
        geom_point(aes(color=Feature))+
        coord_flip()+ylim(-max.X,max.X)+
        scale_size(limits = c(0.001,0.6),range=c(0.5,15))+
         guides(fill=FALSE,color=F)+
        theme_bw()
     pdf(file = paste0(i,".pdf"), width = 5, height = 2.5)
        print(g)
      dev.off()
}


#The End
