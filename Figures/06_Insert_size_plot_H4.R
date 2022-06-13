#!/usr/bin/env Rscript

outpath <- "path/to/plots"
colors = cbind("#00d455", "#ff00cc", "#8800aa","#00ccff")

################################################
######## Each samples in it's own plot #########
################################################

dirs <- list.dirs("path/to/RUN/03_QUALIMAP", recursive = FALSE)
names <- c("1.5U","100U","25U","6.25U")

counter = 1

par(mfrow = c(2, 2))
for (bamqc in dirs){
  file <- file.path(bamqc, "/raw_data_qualimapReport/insert_size_histogram.txt")
  hist <- read.table(file)
  hist["V2"] = lapply(hist["V2"], function(x) x/sum(hist["V2"]))
  plot(hist, type = "l", ylim = c(0, 0.02), col = colors[counter], lwd = 2, main = names[counter], xlab = "fragment size", ylab = "frequency")
  counter = counter + 1
}

######################################
###### All samples in one plot #######
######################################

path <- "path/to/RUN/03_QUALIMAP"
file_list <- list.files(path, "insert_size_histogram.txt$",recursive = TRUE,full.names=TRUE)

multmerge = function(file_list){
  datalist = lapply(file_list, function(x){read.table(file=x,header=F)})
  Reduce(function(x,y) merge(x,y,by="V1",all=TRUE), datalist)
}

all <- multmerge(file_list)
names <- c("Length","1.5U","100U","25U","6.25U")
colnames(all) <- names

all[is.na(all)] <- 0

#Normalize read counts to sequencing depth
all$`1.5U` <- as.numeric(all$`1.5U`) / sum(all$`1.5U`)
all$`6.25U` <- as.numeric(all$`6.25U`) / sum(all$`6.25U`)
all$`25U` <- as.numeric(all$`25U`) / sum(all$`25U`)
all$`100U` <- as.numeric(all$`100U`) / sum(all$`100U`)

head(all)

library(reshape2)
library(ggplot2)

setwd(outpath)
theme_set(
  theme_classic()
)

insert_sizes <- melt(all,id.vars = "Length")

p <- ggplot(insert_sizes, aes(y=value,x=Length,color=variable)) + 
  geom_line(lwd=1.5) + scale_color_manual(values=colors)

pdf("Insert_size_plot.pdf")
p + xlab("Insert size [nt]") + 
  theme(axis.text.x = element_text(colour = 'black', size = 14),
          axis.title.x = element_text(colour = 'black',size = 14)) + 
  ylab("Proportion of reads") +
  theme(axis.text.y = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black',size = 14)) +
  labs(color = "MNase titrations") +
  theme(legend.title = element_text(color = "black", size = 14), 
        legend.text = element_text(color = "black", size = 14))
dev.off()
