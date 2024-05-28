library(ggplot2)
library(dplyr)
library(data.table)


plotAHMM<-function(sample){
  AFP<-read.csv(paste0("filter_0.95/",sample,".posterior.filter0.95"),sep = "\t")
  AFP$prob<-round((AFP$X2.0),2)+(1/2)*round((AFP$X1.1),2)
  
  # Round to closest 0,.5 or 1
  AFP$prob <- round(AFP$prob * 2) / 2
  FLprop<-mean(AFP$prob)
  TXprop<-1-FLprop
  AFP<-AFP[AFP$chrom%in%scaffolds100Kb$chrom,]
  
  # get chromosome length for each chromosome
  AFP<-AFP %>% group_by(chrom) %>% mutate(chr_len=n())
  
  # change the chromosome name (string) to factor
  AFP$chrom <- factor(AFP$chrom, levels=unique(AFP$chrom))
  
  #Create chunks of ancetry by joining adjacent SNPs with the same ancestry
  result <- AFP %>%
    mutate(group_id = rleid(prob)) %>%
    group_by(chrom, group_id) %>%
    summarise(start = min(position), end = max(position), prob = first(prob)) %>%
    ungroup()

  #Chunk length
  result$chunk<-result$end-result$start
  #Modify if you want to exclude chunks smaller than certain number of bp
  smaller<-result[result$chunk>0,]
  subset<-smaller
  
  AFP<-subset%>% group_by(chrom) %>%mutate(chr_len = max(end))

  # line up positions on each chromosome
  AFP<-AFP %>% group_by(chrom)%>% arrange(start, .by_group = TRUE)
  chroms<-data.frame(unique(AFP[c("chrom","chr_len")]))
  chroms$cumlen <-cumsum(as.numeric(chroms$chr_len))-chroms$chr_len
  
  # set up how to line data up along the plotting X axis
  AFP <- merge(AFP, chroms, by = c("chrom","chr_len"), all.x = TRUE)
  AFP$startA<-AFP$start+AFP$cumlen
  AFP$endA<-AFP$end+AFP$cumlen
  
  #AFP<-AFP%>% group_by(chrom) %>% arrange(desc(chr_len), .by_group = FALSE)
  AFP$chrom <- factor(AFP$chrom, levels=chroms$chrom)
  
  # Set up colors for each chromosome
  unique_chroms <- unique(AFP$chrom)
  chrom_colors <- rep(c("#62d0f5", "#1d7052"), length.out = length(unique_chroms))
  
  # Plot the background rectangles
  bg_rectangles <- data.frame(
    xmin = rep(-Inf, length(unique_chroms)),
    xmax = rep(Inf, length(unique_chroms)),
    ymin = rep(-Inf, length(unique_chroms)),
    ymax = rep(Inf, length(unique_chroms)),
    chrom = unique_chroms
  )
   
  dkcol="#8bae24"
  ltcol="#bbdc58"

  #Calculate proportion of ancestry based on chunks
  proportions<-AFP%>% 
    group_by(prob) %>% 
    summarize(total = sum(chunk)) 
  proportions$prop<-proportions$total/sum(proportions$total)
  propFL<-(proportions[proportions$prob==0.5,]$prop)/2+proportions[proportions$prob==1,]$prop
  propTX<-(proportions[proportions$prob==0.5,]$prop)/2+proportions[proportions$prob==0,]$prop
  propFLFL<-proportions[proportions$prob==1,]$prop
  propFLTX<-proportions[proportions$prob==.5,]$prop
  propTXTX<-proportions[proportions$prob==0,]$prop
  
  dfNew<-data.frame(Name=sample, 
                          FL=round(propFL, 4)*100, 
                          TX=round(propTX, 4)*100,
                          FLFL=round(propFLFL, 4)*100,
                          FLTX=round(propFLTX, 4)*100,
                          TXTX=round(propTXTX, 4)*100)
  
  #Modify if you want to only look at first scaffold. Example 200 scaffolds
  extractChrom<-unique(AFP$chrom)[1:200]
  
  AFP<-AFP[AFP$chrom%in%extractChrom,]
  axis_set <- AFP %>% 
    group_by(chrom) %>% 
    summarize(center = mean(startA)) %>% arrange(center, .by_group = FALSE)
  # Assuming axis_set is a data frame containing information about the axis
  # You may need to replace the following code based on your actual data structure
  axis_set$chrom <- factor(axis_set$chrom, levels = levels(AFP$chrom))
  axis_set$start <- sapply(axis_set$chrom, function(chrom) min(AFP$startA[AFP$chrom == chrom]))
  axis_set$end <- sapply(axis_set$chrom, function(chrom) max(AFP$endA[AFP$chrom == chrom]))
  axis_set<-droplevels(axis_set)
  # Calculate adjusted breaks
  breaks <- axis_set$start  # Using start positions as breaks
  
  ggplot(AFP, aes(x = startA, xend = endA, y = prob, yend = prob, color = chrom)) +
    geom_segment(size = 5) +
    theme_classic(base_size = 10) +
    scale_x_continuous(position = "top",
                       breaks = breaks,
                       labels = paste0("         ",1:length(axis_set$chrom)),
                       expand = c(0, 0)  # No expansion of axis limits
    ) +
    scale_y_continuous(
      breaks = c(0,.5,1),
      labels = c("TX/TX", "FL/TX","FL/FL"),
      expand = c(.1,.1)
    ) +
    scale_color_manual(values = rep(c(dkcol, ltcol), unique(length(AFP$chrom)))) +
    labs(caption = paste0(sample, " FL : ", round(propFL, 4) * 100, "%, TX : ", round(propTX, 4) * 100, "%"),
         y = "Ancestry")+
    theme(legend.position = "none",
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 8,color="black", 
                                     vjust = c(c(2.2,1.5,1.5,1.5,
                                                 1.5,1.5,1.5,1.25,1.25), rep(1,39 )), 
                                     hjust = 0),
          axis.title.x = element_blank(),
          plot.caption= element_text(hjust = 0.5, size = 10,vjust=1) )
  ggsave(paste0("aHMM_",sample,".png"), width=12, height=1.3)
  return(dfNew)
}
