###ACHI VIS data analysis
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(plyr)
library(data.table)
library(stringr)
library(reshape2)
#library(svMisc)
#library(clipr)
library(qvalue)
library(scales)
library(ggiraph)
library(ggiraphExtra)
library(plotrix)
library(AICcmodavg)

library(cowplot)



path="/Users/Nicolas/WORK/Achi Vis"


DF <- read.delim("/Users/Nicolas/WORK/Achi Vis/achi_vis_testes_strand.6.36.tpm.txt",
                 header = TRUE,
                 #colClasses = "character", # force a col_type
                 sep = "\t")



TFAav_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_regulated_tfgenes_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")

TFAavo_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_not_regulated_tf_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")

DNGav_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_regulated_denovo_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")

DNGavo_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_not_regulated_denovo_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")

TSPav_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_regulated_testis_biased_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")

TSPavo_info <- read.delim("/Users/Nicolas/WORK/Achi Vis/visachi_not_regulated_testis_biased_fbid.txt",
                   header = FALSE,
                   #colClasses = "character", # force a col_type
                   sep = "\t")


#OTher
DNG_ID_index_i <- read.delim("/Users/Nicolas/WORK/Achi Vis/DNG_current_symbol.txt",
                             header = TRUE,
                             #colClasses = "character", # force a col_type
                             sep = "\t")
DNG_ID_index <- subset(DNG_ID_index_i, select=c(Gene_ID,current_symbol))

DNG_age_info_i <- read.delim("/Users/Nicolas/WORK/Achi Vis/supplementary_File_S2_denovo_property_with_otherISDs.txt",
                             header = TRUE,
                             #colClasses = "character", # force a col_type
                             sep = "\t")

DNG_age_info <- subset(DNG_age_info_i, select=c(Gene_ID,lineage_age))

#Data formatting

DF1 <- DF[1:16]
#DF1[,2:16] <- as.numeric(DF1[,2:16])

setnames(DF1, old=c(names(DF1)), new=c("Gene_ID","A2V1_1","A2V1_2","A2V1_3",
                                                 "A1V2_1","A1V2_2","A1V2_3",
                                                 "A2V2_1","A2V2_2","A2V2_3",
                                                 "A3V3_1","A3V3_2","A3V3_3",
                                                 "A1V1_1","A1V1_2","A1V1_3"))


DF1$Mean_TPM <- apply(DF1[,2:16], 1, mean, na.rm=TRUE)
DF1$Median_TPM <- apply(DF1[,2:16], 1, median, na.rm=TRUE)

TF_list_av <- c(dplyr::pull(TFAav_info, V1),dplyr::pull(TFAavo_info, V1))
DNG_list_av <- c(dplyr::pull(DNGav_info, V1),dplyr::pull(DNGavo_info, V1))
TSP_list_av <- c(dplyr::pull(TSPav_info, V1),dplyr::pull(TSPavo_info, V1))
REG_list_av <- c(dplyr::pull(DNGav_info, V1),dplyr::pull(TFAav_info, V1),dplyr::pull(TSPav_info, V1))

DF1 <- within(DF1, DNG_av <- ifelse(Gene_ID%in%DNG_list_av ,TRUE,FALSE))
DF1 <- within(DF1, TFA_av <- ifelse(Gene_ID%in%TF_list_av ,TRUE,FALSE))
DF1 <- within(DF1, TSP_av <- ifelse(Gene_ID%in%TSP_list_av ,TRUE,FALSE))
DF1 <- within(DF1, REG_av <- ifelse(Gene_ID%in%REG_list_av ,TRUE,FALSE))




VIS_exp <- as.numeric(subset(DF1, (Gene_ID %in% c( "FBgn0033748")))[,2:16])
ACH_exp <- as.numeric(subset(DF1, (Gene_ID %in% c( "FBgn0033749")))[,2:16])
Gen_grp <- as.character(c("a","a","a","b","b","b","c","c","c","d","d","d","e","e","e"))

########NO TRANSFORMATION


DF2 <- DF1

DF2$Zero <- "N"

DF2$Slope_V <-"N"
DF2$p_V <- "N"

DF2$Slope_A <-"N"
DF2$p_A <- "N"

DF2$Slope_I <-"N"
DF2$p_I <- "N"

Rsq <- "N"
Padj <- "N"


#count the zeros
for(i in 1:nrow(DF2)){
  if(DF1$Mean_TPM[i]==0) {
    DF2$Zero[i] <- 15
    DF2[i,2:16] <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
  
  }
  else{
  
    DF2$Zero[i] <- sum(DF1[i,2:16]==0)
    DF2[i,2:16] <- gsub("Inf", NA, (DF1[i,2:16])) #apply transformation here if necessary
  }
  
}

Max_zero <- 3#Max count of 0 accepted

#Define a function to extract adjpval from summary
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#apply linear model
for(i in 1:nrow(DF2)){
  if(as.numeric(DF2$Zero[i])> Max_zero ) {
    
    DF2$Slope_V[i] <-"NA"
    DF2$p_V[i] <- "NA"
    
    DF2$Slope_A[i] <-"NA"
    DF2$p_A[i] <- "NA"
    
    DF2$Slope_I[i] <-"NA"
    DF2$p_I[i] <- "NA"
    
    DF2$Rsq[i] <- "NA"
    DF2$Padj[i] <- "NA"
    
    
  }
  else{
    Yg <- as.numeric(c(DF2[i,2:16]))
    TT <- data.frame("X"=VIS_exp,"Z"=ACH_exp, "Y"=Yg, "Gen_grp"=Gen_grp)
    TT[TT==0] <- NA
    #1)full
    Mo_full <- lm(Y~X*Z,TT,na.action=na.omit) #linear model
    SMo_full <- summary(Mo_full)
   
    DF2$Slope_V_full[i] <-SMo_full$coefficients[2,1] #slope for Vis factor
    DF2$p_V_full[i] <- SMo_full$coefficients[2,4] #pvalue
    
    DF2$Slope_A_full[i] <-SMo_full$coefficients[3,1] #slope for Achi factor
    DF2$p_A_full[i] <- SMo_full$coefficients[3,4]
    
    DF2$Slope_I_full[i] <-SMo_full$coefficients[4,1]#slope for VxA interaction factor
    DF2$p_I_full[i] <- SMo_full$coefficients[4,4]
    
    
    DF2$Rsq_full[i] <- SMo_full$adj.r.squared
    DF2$Padj_full[i] <- overall_p(Mo_full)# adjusted pval
   
    
    #2)vis + achi
    Mo_av <- lm(Y~X+Z,TT,na.action=na.omit) #linear model
    SMo_av <- summary(Mo_av)
    
    DF2$Slope_V_av[i] <-SMo_av$coefficients[2,1] #slope for Vis factor
    DF2$p_V_av[i] <- SMo_av$coefficients[2,4] #pvalue
    
    DF2$Slope_A_av[i] <-SMo_av$coefficients[3,1] #slope for Achi factor
    DF2$p_A_av[i] <- SMo_av$coefficients[3,4]
    
    
    DF2$Rsq_av[i] <- SMo_av$adj.r.squared
    DF2$Padj_av[i] <- overall_p(Mo_av)
    
   
    #3)vis
    Mo_v <- lm(Y~as.numeric(X),TT,na.action=na.omit) #linear model
    SMo_v <- summary(Mo_v)
    
    DF2$Slope_V_v[i] <-SMo_v$coefficients[2,1] #slope for Vis factor
    DF2$p_V_v[i] <- SMo_v$coefficients[2,4] #pvalue
    
    
    DF2$Rsq_v[i] <- SMo_v$r.squared
    DF2$Rsq_adj_v[i] <- SMo_v$adj.r.squared
    DF2$Padj_v[i] <- overall_p(Mo_v)
    
    
    #4) achi 
    Mo_a <- lm(Y~as.numeric(Z),TT,na.action=na.omit) #linear model
    SMo_a <- summary(Mo_a)
    
    DF2$Slope_A_a[i] <-SMo_a$coefficients[2,1] #slope for Vis factor
    DF2$p_A_a[i] <- SMo_a$coefficients[2,4] #pvalue
    
    DF2$Rsq_a[i] <- SMo_a$adj.r.squared
    DF2$Padj_a[i] <- overall_p(Mo_a)
    
  
    
    models <- list(Mo_full, Mo_av, Mo_v, Mo_a) 
    mod.names <- c('V.A.I', 'V.A', 'V', 'A') 
    AIC_scores <- aictab(cand.set = models, modnames = mod.names)
    DF2$AIC_1[i] <- AIC_scores$Modnames[1]
    DF2$AIC_2[i] <- AIC_scores$Modnames[2]
    DF2$AIC_3[i] <- AIC_scores$Modnames[3]
    DF2$AIC_4[i] <- AIC_scores$Modnames[4]
    
    DF2$AIC_D2[i] <- AIC_scores$Delta_AICc[2]
    DF2$AIC_D3[i] <- AIC_scores$Delta_AICc[3]
    DF2$AIC_D4[i] <- AIC_scores$Delta_AICc[4]
  }
  
}




#Count where T~V is best model

DF2$t.v <- ifelse(DF2$AIC_1=="V" , "TRUE", 
                  ifelse(DF2$AIC_2=="V" & DF2$AIC_D2<2, TRUE, 
                                                  ifelse(DF2$AIC_3=="V" & DF2$AIC_D3<2, TRUE, 
                                                         ifelse(DF2$AIC_4=="V"& DF2$AIC_D4<2, TRUE, FALSE))))


#Best model
ti="rawT=V" # title
# define a cutoff for minimal Mean_TPM
cu <- 0 

# use the cutoff; remove Vis and Achi from the table
DF4 <- subset(DF2, Mean_TPM >cu & Median_TPM >cu & !(Gene_ID %in% c( "FBgn0033748","FBgn0033749")))





#VISACHI
St <- 0.15
DF4tsp <- subset(DF4, TSP_av==TRUE & abs(Slope_V_v) >St, select=c(Gene_ID, Rsq_v, REG_av))
DF4tsp$set <- "1tsp"
DF4dng <- subset(DF4, DNG_av==TRUE & abs(Slope_V_v) >St, select=c(Gene_ID, Rsq_v, REG_av))
DF4dng$set <- "2dng"
DF4tfa <- subset(DF4, TFA_av==TRUE  & abs(Slope_V_v) >St, select=c(Gene_ID, Rsq_v, REG_av))
DF4tfa$set <- "3tfa"
DF4fig_i <- rbind(DF4tsp, DF4dng, DF4tfa )
DF4fig <- arrange(transform(DF4fig_i,set=factor(set,levels=c("1tsp","2dng","3tfa"))),set)

figRsq <- ggplot( DF4fig ,aes(x=REG_av, y=(as.numeric(Rsq_v)), fill=REG_av)) +
  geom_violin(position=position_dodge(0.8)) +
  scale_fill_manual(values=c( "#BEBEBE", "#52a3c2" )) +
  geom_boxplot(width=0.085,size=0.5,position=position_dodge(0.8),outlier.shape = NA,fill="white") +
  stat_compare_means(label = "p.format", comparisons = list(c("TRUE","FALSE")), size=3) + 
  ylim(-0.1,1.1)+
  scale_x_discrete ( name="", labels = c("TRUE"="Vis regulon", "FALSE"="Not Vis regulon")) +
  ylab("R-square between\nVis and target gene ") + 
  facet_grid ( .~ set, labeller = labeller(set = c("1tsp"="Testis specific genes","2dng"="De novo genes","3tfa"="Transcription factors")))+
  theme(strip.background = element_rect( fill="white"),
        strip.placement = "outside",
        panel.border = element_blank(), 
        panel.background = element_rect (fill = "white", colour= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Arial"),
        axis.text=element_text(size=8), 
        axis.title=element_text(size=10),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line.y = element_line (linewidth = 0.5, linetype= "solid", colour= "black"),
        axis.line.x = element_line (linewidth = 0.5, linetype= "solid", colour= "white"),
        legend.position='none',
        
        plot.title = element_text(size=10))


figRsq
ggsave(filename=paste("Rsq_Vis_to_target_AIC","ALL","VisONLY","2024_11_25.fig.eps",sep="_"), 
       plot= figRsq   , 
       path= path, 
       width = 17, 
       height = 8, 
       units = "cm")
P1L <- figRsq

aggregate(Rsq_v ~ REG_av, data = subset(DF4fig, set=="1tsp"), 
          FUN = function(x) c( N= length(x) , mean = mean(x), se = std.error(x)))

wilcox.test(Rsq_v ~ REG_av, data=subset(DF4fig, set=="1tsp"), alternative = "two.sided") #alternative: "two.sided";"greater";"less"

aggregate(Rsq_v ~ REG_av, data = subset(DF4fig, set=="2dng"), 
          FUN = function(x) c(N=length(x),mean = mean(x), se = std.error(x)))
wilcox.test(Rsq_v ~ REG_av, data=subset(DF4fig, set=="2dng"), alternative = "two.sided") #alternative: "two.sided";"greater";"less"


aggregate(Rsq_v ~ REG_av, data = subset(DF4fig, set=="3tfa"), 
          FUN = function(x) c(N=length(x),mean = mean(x), se = std.error(x)))
wilcox.test(Rsq_v ~ REG_av, data=subset(DF4fig, set=="3tfa"), alternative = "two.sided") #alternative: "two.sided";"greater";"less"





##### FDR cutoff 
DF4z <- subset(DF4, !is.na(as.numeric(Padj_v)))
qobj <- qvalue(p = as.numeric(DF4z$Padj_v), fdr.level = 0.05) #AT 5% FDR
DF4z$qval005 <- qobj$qvalues
DF4z$SIG005 <- qobj$significant
qobj <- qvalue(p = as.numeric(DF4z$Padj_v), fdr.level = 0.10) #AT 5% FDR
DF4z$qval010 <- qobj$qvalues
DF4z$SIG010 <- qobj$significant


DF5 <- DF4z



# %age genes with positive slope
nrow(subset(DF5, SIG005==TRUE &  Slope_V_v>0))
nrow(subset(DF5, SIG005==TRUE &  Slope_V_v>0))/nrow(subset(DF5, SIG005==TRUE))

#CANDIDATE PLOTS

DF_allSigDNGii <- subset(DF5, DNG_av==TRUE & SIG005==TRUE)
DF_allSigDNGi <- left_join(DF_allSigDNGii,DNG_ID_index, by="Gene_ID" , copy=FALSE)
DF_allSigDNG <- left_join(DF_allSigDNGi,DNG_age_info, by="Gene_ID" , copy=FALSE)



#candidate plotting 1
Tg <- "FBgn0265423"  #FBgn0036443
rw <- which(DF5$Gene_ID==Tg)
Tgsy <- subset(DF_allSigDNG, Gene_ID==Tg, select=c(current_symbol))[1,1]
YTg <- (as.numeric(c(DF5[rw,2:16])))
TTg <- data.frame("X"=VIS_exp,"Z"=ACH_exp, "Y"=YTg, "Gen_grp"=Gen_grp)
Label1 <- paste(paste("R2=",round(DF5$Rsq_v[rw],2),sep=""),
                 #paste("a=",round(DF_allSigDNG$Slope_V_v[j],2),sep=""),
                 paste("p=",formatC(DF5$Padj_v[rw],format="e",digits=1),sep=""),
                 # paste("Reg=",DF_allSigDNG$REG_av[j],sep=""),
                 #paste(DF_allSigDNG$lineage_age[j],sep=""),
                 sep="\n")
MoTg <- lm(Y~X,TTg) #linear model
SMoTg <- summary(MoTg)
SMoTg

p4 <-  ggplot(TTg,aes(y=Y,x=X))+
  geom_point(color="#BEBEBE")+
  stat_smooth(method="lm",se=FALSE,col="#52a3c2")+
  ggtitle(paste("", Tgsy, sep=" ")) +
  #geom_text(x=50, y=55,  size=3, family = "Arial", label=Label1)+
  scale_x_continuous(breaks=c(25,75,125)) +
  scale_y_continuous(breaks=c(15,30,45,60)) +
  ylab(paste(Tgsy, "TPM", sep=" ")) +
  xlab("Vis TPM") 
p4

#candidate plotting 2
Tg <- "FBgn0032658"  #FBgn0036443
rw <- which(DF5$Gene_ID==Tg)
Tgsy <- subset(DF_allSigDNG, Gene_ID==Tg, select=c(current_symbol))[1,1]
YTg <- (as.numeric(c(DF5[rw,2:16])))
TTg <- data.frame("X"=VIS_exp,"Z"=ACH_exp, "Y"=YTg, "Gen_grp"=Gen_grp)
Label1 <- paste(paste("R2=",round(DF5$Rsq_v[rw],2),sep=""),
                #paste("a=",round(DF_allSigDNG$Slope_V_v[j],2),sep=""),
                paste("p=",formatC(DF5$Padj_v[rw],format="e",digits=1),sep=""),
                # paste("Reg=",DF_allSigDNG$REG_av[j],sep=""),
                #paste(DF_allSigDNG$lineage_age[j],sep=""),
                sep="\n")
MoTg <- lm(Y~X,TTg) #linear model
SMoTg <- summary(MoTg)
SMoTg

p5 <-  ggplot(TTg,aes(y=Y,x=X))+
  geom_point(color="#BEBEBE")+
  stat_smooth(method="lm",se=FALSE,col="#52a3c2")+
  ggtitle(paste("", Tgsy, sep=" ")) +
  #geom_text(x=110, y=60, size=3, family = "Arial", label=Label1)+
  scale_x_continuous(breaks=c(25,75,125)) +
  scale_y_continuous(breaks=c(40,60,80,100)) +
  ylab(paste(Tgsy, "TPM", sep=" ")) +
  xlab("Vis TPM") 
p5 
   
#candidate plotting 3

Tg <- "FBgn0028943"  #FBgn0036443
rw <- which(DF5$Gene_ID==Tg)
Tgsy <- subset(DF_allSigDNG, Gene_ID==Tg, select=c(current_symbol))[1,1]
YTg <- (as.numeric(c(DF5[rw,2:16])))
TTg <- data.frame("X"=VIS_exp,"Z"=ACH_exp, "Y"=YTg, "Gen_grp"=Gen_grp)
Label1 <- paste(paste("R2=",round(DF5$Rsq_v[rw],2),sep=""),
                #paste("a=",round(DF_allSigDNG$Slope_V_v[j],2),sep=""),
                paste("p=",formatC(DF5$Padj_v[rw],format="e",digits=1),sep=""),
                # paste("Reg=",DF_allSigDNG$REG_av[j],sep=""),
                #paste(DF_allSigDNG$lineage_age[j],sep=""),
                sep="\n")
MoTg <- lm(Y~X,TTg) #linear model
SMoTg <- summary(MoTg)
SMoTg

p6 <-  ggplot(TTg,aes(y=Y,x=X))+
  geom_point(color="#BEBEBE")+
  stat_smooth(method="lm",se=FALSE,col="#52a3c2")+
  ggtitle(paste("", Tgsy, sep=" ")) +
  #geom_text(x=100, y=360,size=3, family = "Arial", label=Label1)+
  scale_x_continuous(breaks=c(25,75,125)) +
  scale_y_continuous(breaks=c(100,200,300,400)) +
  ylab(paste(Tgsy, "TPM", sep=" ")) +
  xlab("Vis TPM") 
p6

  Theme_layout_Regs <- theme(strip.background = element_rect( fill="white"),
                        strip.placement = "outside",
                        panel.border = element_blank(), 
                        panel.background = element_rect (fill = "white", colour= "white"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        text = element_text(family = "Arial"),
                        axis.text=element_text(size=8), 
                        axis.title=element_text(size=8),
                        axis.ticks.y = element_line(color = "black"),
                        axis.ticks.x = element_line(color = "black"),
                        axis.line.y = element_line (linewidth = 0.5, linetype= "solid", colour= "black"),
                        axis.line.x = element_line (linewidth = 0.5, linetype= "solid", colour= "white"),
                        legend.text=element_text(size=8),
                        legend.title=element_text(size=8),
                        plot.title = element_text(size=8, hjust = 0.5)) 
P4 <-p4+ Theme_layout_Regs
P5 <-p5+ Theme_layout_Regs
P6 <-p6+ Theme_layout_Regs  

########  FIGURE PLOTTING


FIG4 <- ggdraw()+
  draw_plot(P1L, x = 0    , y = .5, width = 0.95 , height = .49) +
  #draw_plot(P2L, x = 0 , y = 0.5, width = 0.9 , height = .25) +
  #draw_plot(P3L, x = 0    , y = 0.25  , width = 0.9 , height = .25) +
  #draw_plot(Legend, x = 0.91    , y = 0.5  , width = 0.09 , height = .25) +
  draw_plot(P4, x = 0 , y = 0  , width = .32 , height = .49) +
  draw_plot(P5, x = 0.33 , y = 0  , width = .32 , height = .49) +
  draw_plot(P6, x = 0.66 , y = 0  , width = .32 , height = .49) +
  #draw_plot(p5, x = 0.65, y = 0, width = 0.32, height = .49) +
  draw_plot_label(label = c("A", "B", "C", "D"), #, "FIGURE 4"
                  size = 12,
                  x = c(0,   0 , 0.33, 0.66), #, 0.87
                  y = c(1, 0.5 , 0.5, 0.5)) #, 0.05

FIG4
#ggsave(paste(path,"_Figure_RNAseq2_2024_12_10.eps",sep=""), plot= FIG1 , width = 16, height = 11, units = "cm")


##############
#Supplementary Figure S8

p1 <- ggplot(subset(DF5, TSP_av==TRUE), aes(x = as.numeric(Rsq_v), y = (abs(Slope_V_v)))) + 
  geom_density_2d_filled(contour_var = "ndensity") + 
  ggtitle("De novo genes") +
  scale_x_continuous(name = "Regression R-square",  limits=c(0,0.9) ,breaks=c(0, 0.15, 0.30, 0.45, 0.60 ,0.75, 0.9)) +
  scale_y_continuous(name = "Regression slope", trans = "log2", limits=c(5E-5,5E1) ,breaks=c(5E-5, 5E-4, 5E-3, 5E-2,5E-1,5E0, 5E1)) +  
  facet_grid ( .~ REG_av, labeller = labeller(REG_av = c("TRUE"="Vis regulon","FALSE"="Not Vis regulon")))


p2 <- ggplot(subset(DF5, DNG_av==TRUE), aes(x = as.numeric(Rsq_v), y = (abs(Slope_V_v)))) + 
  geom_density_2d_filled(contour_var = "ndensity") + 
  ggtitle("De novo genes") +
  scale_x_continuous(name = "Regression R-square",  limits=c(0,0.9) ,breaks=c(0, 0.15, 0.30, 0.45, 0.60 ,0.75, 0.9)) +
  scale_y_continuous(name = "Regression slope", trans = "log2", limits=c(5E-5,5E1) ,breaks=c(5E-5, 5E-4, 5E-3, 5E-2,5E-1,5E0, 5E1)) +  
  facet_grid ( .~ REG_av, labeller = labeller(REG_av = c("TRUE"="Vis regulon","FALSE"="Not Vis regulon")))

p3 <- ggplot(subset(DF5, TFA_av==TRUE), aes(x = as.numeric(Rsq_v), y = (abs(Slope_V_v)))) + 
  geom_density_2d_filled(contour_var = "ndensity") + 
  ggtitle("Transcription Factors") +
  scale_x_continuous(name = "Regression R-square",  limits=c(0,0.9) ,breaks=c(0, 0.15, 0.30, 0.45, 0.60 ,0.75, 0.9)) +
  scale_y_continuous(name = "Regression slope", trans = "log2", limits=c(5E-5,5E1) ,breaks=c(5E-5, 5E-4, 5E-3, 5E-2,5E-1,5E0, 5E1)) +  
  facet_grid ( .~ REG_av, labeller = labeller(REG_av = c("TRUE"="Vis regulon","FALSE"="Not Vis regulon")))

Theme_layout_DPs <- theme(strip.background = element_rect( fill="white"),
                          strip.placement = "outside",
                          panel.border = element_blank(), 
                          panel.background = element_rect (fill = "white", colour= "white"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          text = element_text(family = "Arial"),
                          axis.text=element_text(size=8), 
                          axis.title=element_text(size=8),
                          axis.ticks.y = element_line(color = "black"),
                          axis.ticks.x = element_line(color = "black"),
                          axis.line.y = element_line (linewidth = 0.5, linetype= "solid", colour= "black"),
                          axis.line.x = element_line (linewidth = 0.5, linetype= "solid", colour= "white"),
                          legend.text=element_text(size=6),
                          legend.title=element_text(size=8),
                          plot.title = element_text(size=10, hjust = 0.5))   

P1 <- p1 + Theme_layout_DPs
P2 <- p2 + Theme_layout_DPs
P3 <- p3 + Theme_layout_DPs

Legend <-get_legend(P2)
P1L <- P1 + theme(legend.position='none')
P2L <- P2 + theme(legend.position='none')
P3L <- P3 + theme(legend.position='none')


FIG_S8 <- ggdraw()+
  draw_plot(P1L, x = 0    , y = .66, width = 0.89 , height = .32) +
  draw_plot(P2L, x = 0 , y = 0.33, width = 0.89 , height = .32) +
  draw_plot(P3L, x = 0    , y = 0  , width = 0.89 , height = .32) +
  draw_plot(Legend, x = 0.91    , y = 0.33  , width = 0.08 , height = .25) +
  draw_plot_label(label = c("A", "B", "C"),
                  size = 12,
                  x = c(0, 0,  0 ),
                  y = c(1, 0.66, 0.33 ))


FIG_S8
ggsave(paste(path,"AchiVis_Figure_S8_2024_12_12.eps",sep=""), plot= FIG_S8 , width = 16, height = 14, units = "cm")
