getwd()
#wd="C:/Users/1/Documents/cov_cal_outgroup", under this dir there are 10 .csv files
library(ggplot2)
library(reshape2)

setwd("C:/Users/1/Documents/cov_cal_outgroup")
male1<-read.csv("RNDAG02-C_FP245.yag.bam_window.csv", header=FALSE)
male2<-read.csv("RNDAG02-D_UCFP160.yag.bam_window.csv", header=FALSE)
male3<-read.csv("RNDAG02-F_UCFP191.yag.bam_window.csv", header=FALSE)
male4<-read.csv("RNDAG02-N_UCFP244.yag.bam_window.csv", header=FALSE)
male5<-read.csv("RNDAG02-X_UCFP287.yag.bam_window.csv", header=FALSE)
female1<-read.csv("RNDAG02-AB_UCFP304.yag.bam_window.csv", header=FALSE)
female2<-read.csv("RNDAG02-H_UCFP217.yag.bam_window.csv", header=FALSE)
female3<-read.csv("RNDAG02-L_UCFP229.yag.bam_window.csv", header=FALSE)
female4<-read.csv("RNDAG02-T_UCFP279.yag.bam_window.csv", header=FALSE)
female5<-read.csv("RNDAG02-Y_UCFP289.yag.bam_window.csv", header=FALSE)

# position end - position start = size of window
process<-function(female){
  female$sizewin<-female$V3-female$V2
  female<-female[female$sizewin>10000,]
  # per bp coverage = coverage of the whole window / size of window
  female$coverage_per_site<-female$V4/female$sizewin
  # normalize the per bp coverage by the overall coverage of that individual
  female$normalized_cov<-female$coverage_per_site/mean(female$coverage_per_site)
  return(female)
}

f1<-process(female1)
f2<-process(female2)
f3<-process(female3)
f4<-process(female4)
f5<-process(female5)

m1<-process(male1)
m2<-process(male2)
m3<-process(male3)
m4<-process(male4)
m5<-process(male5)


library(dplyr)
# merge 5 males together; merge 5 females together
male_table<-merge(m1, m2, by=c("V1", "V2", "V3","sizewin"),suffixes = c("_m1","_m2"))
male_table<-merge(male_table, m3, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_m3"))
male_table<-merge(male_table, m4, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_m4"))
male_table<-merge(male_table, m5, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_m5"))

#calculate the mean coverage per bp for each window for 5 males
male_table<-male_table %>% mutate(male_cov=(normalized_cov_m1+normalized_cov_m2+normalized_cov+normalized_cov_m4+normalized_cov_m5)/5)



female_table<-merge(f1, f2, by=c("V1", "V2", "V3","sizewin"),suffixes = c("_f1","_f2"))
female_table<-merge(female_table, f3, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_f3"))
female_table<-merge(female_table, f4, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_f4"))
female_table<-merge(female_table, f5, by=c("V1", "V2", "V3","sizewin"),suffixes =c("","_f5"))

#calculate the mean coverage per bp for each window for 5 males
female_table<-female_table %>% mutate(female_cov=(normalized_cov_f1+normalized_cov_f2+normalized_cov+normalized_cov_f4+normalized_cov_f5)/5)

# subset the table
male_table<-male_table[,c('V1','V2', 'V3','male_cov')]
female_table<-female_table[,c('V1','V2','V3', 'female_cov')]

# group by chromosome, line up by starting position
male_table<-male_table %>% group_by(V1) %>% arrange(V2, .by_group = TRUE)
female_table<-female_table %>% group_by(V1) %>% arrange(V2, .by_group = TRUE)
# merge male and female data
whole_table<-merge(male_table,female_table,by=c("V1", "V2","V3"))
# get chromosome length for each chromosome
whole_table<-whole_table %>% group_by(V1) %>% mutate(chr_len=max(V3))
# line up based on chromosome size (big to small)
whole_table<-whole_table %>% arrange(desc(chr_len), .by_group = FALSE) %>% group_by(V1)
#whole_table$male_cov<-log(whole_table$male_cov)
#whole_table$female_cov<-log(whole_table$female_cov)

whole_table$SNPs<-1:nrow(whole_table)

axis_set <- whole_table %>% 
  group_by(V1) %>% 
  summarize(center = mean(SNPs))

whole_table$ratio<-whole_table$male_cov/whole_table$female_cov

#hola<-whole_table%>% 
#  group_by(V1) %>% summarize(meanratio=median(ratio))

#sexX<-hola[hola$meanratio<=.6,]
#sexY<-hola[hola$meanratio>=2,]





whole_table$V1 <- factor(whole_table$V1, levels=unique(whole_table$V1))

coverage_plot2<-ggplot(whole_table, aes(whole_table,color=as.factor(V1)))+
  geom_point(aes(x=SNPs, y=ratio), size=1,alpha=.8)+theme_classic(base_size=14)+
  ylim(0, 3)+geom_hline(yintercept = 1.15)+geom_hline(yintercept = 0.85)+
  scale_x_continuous(label = axis_set$V1, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("red", "gray"),unique(length(whole_table$V1))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  xlab("Chromosome")+ylab("coverage_ratio_male/female")  

coverage_plot2
#ggsave("coverage.png", width=10, height=5)


sum(is.na(whole_table$ratio))
# make sure the table doesn't contain NA value
windows2remove<-whole_table[whole_table$ratio<0.85 |whole_table$ratio>1.15, ]
windows2remove<-windows2remove[, c('V1','V2','V3')]
#write.csv(windows2remove, "win2remove_yagAsRef.csv", row.names = FALSE,col.names=FALSE)
# the above part is to paint the overview and output the windows that have strange ratio (potential sex chromosome fragments)


# the below part is to paint the putative X and Y
chr_putative_table<-whole_table%>% group_by(V1) %>% count(V1, name='total')
X_putative<-whole_table[whole_table$ratio<0.85,]
Y_putative<-whole_table[whole_table$ratio>1.15,]
X_putative<-X_putative%>% group_by(V1) %>% count(V1, name='X')
Y_putative<-Y_putative%>% group_by(V1) %>% count(V1, name='Y')
#chr_putative_table<-merge(chr_putative_table, X_putative, by="V1", all=TRUE)
#chr_putative_table<-merge(chr_putative_table, Y_putative, by="V1", all=TRUE)
X_putative<-merge(chr_putative_table, X_putative, by="V1")
Y_putative<-merge(chr_putative_table, Y_putative, by="V1")
#chr_putative_table$X<-chr_putative_table$X/chr_putative_table$total
#chr_putative_table$Y<-chr_putative_table$Y/chr_putative_table$total
X_putative$ratio<-X_putative$X/X_putative$total
Y_putative$ratio<-Y_putative$Y/Y_putative$total
X_putative<-X_putative[X_putative$ratio>=0.8, ]
Y_putative<-Y_putative[Y_putative$ratio>=0.8, ]

#knownsex<-c("NW_020340090.1","NW_020340091.1","NW_020340092.1",
#            "NW_020339946.1","NW_020338877.1")
whole_table$type<-"autosome"
whole_table$type[whole_table$V1%in%X_putative$V1]<-"X?"
whole_table$type[whole_table$V1%in%Y_putative$V1]<-"Y?"
#whole_table$type[whole_table$V1%in%knownsex]<-"X_confirmed"
#whole_table$type[whole_table$V1%in%sexY$V1]<-"Y?"

# this is for highlighting putative X/Y on the plot; autosome will be painted grey, XY will be blue/pink
coverage_plot_XYhighlight<-ggplot(whole_table, aes(whole_table,color=type,alpha=as.factor(V1)))+
  geom_point(aes(x=SNPs, y=ratio), size=1)+theme_classic(base_size=14)+
  ylim(0, 3)+
    scale_x_continuous(label = axis_set$V1, breaks = axis_set$center) +
  scale_color_manual(values = c("gray","blue","hotpink"))+
  scale_alpha_manual(values = rep(c(.6,.8),unique(length(whole_table$V1))))+
  guides(alpha = "none")+
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  xlab("Scaffold")+ylab("coverage_ratio_male/female")  
coverage_plot_XYhighlight
#ggsave("coverage_XYhighlight.png", width=10, height=5)


# in our data analysis, we are going to filter out both
# 1) windows on every scaffold that have ratio >1.15 or <0.85, or
# 2) all windows on all putative X/Y scaffolds; 
  # the definition for putative X(Y) is: over 80% of windows on that scaffold has ratio <0.85 (>1.15)

# for windows having ratio >1.15/<0.85, already stored in the table "windows2remove"

# for putative scaffolds, extract from whole table according to "type"
scaffold_table<-whole_table[,c('V1','V2', 'V3','type')]
scaffold_table<-scaffold_table[scaffold_table$type=="X?"|scaffold_table$type=="Y?",]
scaffold_table<-scaffold_table[,c('V1','V2', 'V3')]

windows2remove_both<-rbind(scaffold_table,windows2remove)
windows2remove_both<-distinct(windows2remove_both)
scaffold_unique_list<-scaffold_table[,'V1']
scaffold_unique_list<-unique(scaffold_unique_list)

write.csv(windows2remove_both,"windows2remove_both_Feb10.csv")
write.csv(scaffold_unique_list, "sex_scaffold_list_Feb10.csv")


#coverage_plot<-ggplot(whole_table, aes(whole_table))+
#  geom_point(aes(x=SNPs, y=male_cov),color="blue",alpha=.6)+
#  geom_point(aes(x=SNPs, y=female_cov),color="pink",alpha=.6, size=1)+theme_classic(base_size=14)+
#  ylim(0, 3)+
#  #  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
#  guides(alpha = "none")+
#  theme(panel.border = element_blank(),
#        panel.grid.major.x = element_blank(),        panel.grid.minor.x = element_blank(),
#        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
#  xlab("Chromosome")+ylab("coverage_ratio_male/female")  
#coverage_plot
#ggsave("coverage_overlap.png", width=10, height=5)



#Yputative_table<-whole_table[whole_table$type=="Y?", ]

#Xputative_table<-whole_table[whole_table$type=="X?", ]
#coverage_plot_XYhighlight<-ggplot(Xputative_table, aes(Xputative_table))+
#  geom_point(aes(x=SNPs, y=male_cov),color="blue",alpha=.6)+
#  geom_point(aes(x=SNPs, y=female_cov),color="pink",alpha=.6, size=1)+theme_classic(base_size=14)+
#  ylim(0, 3)+
#    scale_x_continuous(label = axis_set$V1, breaks = axis_set$center) +
#  guides(alpha = "none")+
#  theme(panel.border = element_blank(),
#        panel.grid.major.x = element_blank(),
#        panel.grid.minor.x = element_blank(),
#        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
#  xlab("Chromosome")+ylab("coverage_ratio_male/female")  
#coverage_plot_XYhighlight


#dat.m <- melt(Xputative,id.vars=c('V1','V2','type'), measure.vars=c('male_cov','female_cov'))

#coverage_boxplot<-ggplot(dat.m)+
#  geom_boxplot(aes(x=V1, y=value,color=variable,fill=type))+
#  scale_color_manual(values=c("blue","pink"))+
#  scale_fill_manual(values=c("black","white"))+
#  theme_classic(base_size=14)+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8))+
#  ylim(0, 3)+
#  xlab("Chromosome")+ylab("Normalized coverage")  
#coverage_boxplot
#ggsave("boxplot_Xcrom.png",width=10,height=5)
