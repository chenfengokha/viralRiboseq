library(Biostrings)
library(plyr)
library(dplyr)
library(parallel)
library(ggplot2)
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAsupply/xijofD.hum.Rdata")
##############D
viruscds <- readDNAStringSet("/home/chenfeng/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/1.reference/genomeold/sce.norm.cds.fa")
virus.DP.exp.tai <- mclapply(mc.cores = 10,1:length(viruscds),function(x){
  tmpseq <- as.character(viruscds[[x]])
  mydf0 <- substring(as.character(tmpseq),seq(1,(nchar(tmpseq)-2),by=3),seq(3,nchar(tmpseq),by=3))
  mydf <- mydf0 %>% table() %>% as.data.frame()
  names(mydf)[1] <- c("codon")
  mydf$codon <- as.vector(mydf$codon)
  xij <- xijofD
  xij$yfre <- mydf$Freq[match(xij$codon,mydf$codon)]
  xij$yfre[which(is.na(xij$yfre))] <- 0
  xij %>% group_by(amino) %>% dplyr::mutate(ny=sum(yfre)) %>% group_by(amino,codon) %>% dplyr::mutate(yij=yfre/ny) -> xij
  xij$yij[which(xij$yij=="NaN")] <- 0
  xij %>%
    group_by(amino) %>% 
    dplyr::summarize(Di.exp=(sum((yij-xij)^2)^0.5),Di.tai=(sum((yij-xijoftai)^2)^0.5)) %>% 
    dplyr::summarize(DP.exp=round(prod(Di.exp)^(1/length(Di.exp)),3),DP.tai=round(prod(Di.tai)^(1/length(Di.tai)),3)) %>% 
    cbind(data.frame(stringsAsFactors = F,genename=names(viruscds)[x]),type="viral",myseq=tmpseq)
  
}) %>% rbind.fill()

humcds <- readDNAStringSet("/home/chenfeng/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/1.reference/genomeold/sce.human.norm.cds.fa")
mydata <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sce_cds_tpm_quant.txt",sep = "\t",stringsAsFactors = F)
hum.DP.exp.tai <- mclapply(mc.cores = 10,1:length(humcds),function(x){
  if(names(humcds)[x] %in% mydata$name){
    tmpseq <- as.character(humcds[[x]])
    mydf0 <- substring(as.character(tmpseq),seq(1,(nchar(tmpseq)-2),by=3),seq(3,nchar(tmpseq),by=3))
    mydf <- mydf0 %>% table() %>% as.data.frame()
    names(mydf)[1] <- c("codon")
    mydf$codon <- as.vector(mydf$codon)
    xij <- xijofD
    xij$yfre <- mydf$Freq[match(xij$codon,mydf$codon)]
    xij$yfre[which(is.na(xij$yfre))] <- 0
    xij %>% group_by(amino) %>% dplyr::mutate(ny=sum(yfre)) %>% group_by(amino,codon) %>% dplyr::mutate(yij=yfre/ny) -> xij
    xij$yij[which(xij$yij=="NaN")] <- 0
    res <- xij %>%
      group_by(amino) %>% 
      dplyr::summarize(Di.exp=(sum((yij-xij)^2)^0.5),Di.tai=(sum((yij-xijoftai)^2)^0.5)) %>% 
      dplyr::summarize(DP.exp=round(prod(Di.exp)^(1/length(Di.exp)),3),DP.tai=round(prod(Di.tai)^(1/length(Di.tai)),3)) %>% 
      cbind(data.frame(stringsAsFactors = F,genename=names(humcds)[x]),type="hum",myseq=tmpseq)
  } else {
    res <- data.frame(stringsAsFactors = F,DP.exp=99, DP.tai=99,genename=99,type=99,myseq=99)
  }
  res
}) %>% rbind.fill() %>% dplyr::filter(DP.exp!=99)


alldp <- rbind(virus.DP.exp.tai,hum.DP.exp.tai)

save(alldp,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.Rdata")

#################################################
#####
mydata <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sce_cds_tpm_quant.txt",sep = "\t",stringsAsFactors = F)
sample <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sample.inf.csv",header = F,stringsAsFactors = F)
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)

datares <- mclapply(mc.cores = 10,2:20,function(x){
  mysample <- colnames(mydata)[x]
  samptmp <- sample$V2[which(sample$V1 == mysample)]
  tmpdata <- mydata[,c(1,x)]
  
  names(tmpdata)[2] <- c("TPM") 
  merge(tmpdata,alldp[,-5],by.x="name",by.y="genename") %>% cbind(data.frame(stringsAsFactors = F,mysample,samptmp))
  #tt %>% group_by(type) %>% dplyr::summarize(rho=cor.test(DP.tai,TPM,method="s")$estimate)
  
}) %>% rbind.fill()

save(datares,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.cis.Rdata")

###plot fig3a-c
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.cis.Rdata")
viral <- datares %>% dplyr::filter(mysample!="GSM2097523_trimmed_cds_tpm" & samptmp != "HEK293T")
viralD <- viral[,c(1,3,4,5)] %>% unique()
source("~/Rfunction/style.print.R")
viral1 <- viralD %>% dplyr::filter(type=="viral")
end1 <- viralD %>% dplyr::filter(type!="viral")
summary(viral1$DP.exp)
summary(end1$DP.exp)
wilcox.test(viral1$DP.exp,end1$DP.exp)$p.value

summary(viral1$DP.tai)
summary(end1$DP.tai)
wilcox.test(viral1$DP.tai,end1$DP.tai)$p.value

viralD %>% ggplot(aes(DP.exp,fill = type))+geom_histogram()+style.print()
viralD %>% ggplot(aes(DP.tai,fill = type))+geom_histogram()+style.print()
viral %>% ggplot(aes(DP.exp,TPM))+geom_point()+
  facet_grid(type~mysample,scale="free")+
  geom_smooth(method = "lm",se=F)
###fig 3c
tt <- viral %>% group_by(type,mysample,samptmp) %>% dplyr::summarize(rho.exp=cor.test(DP.exp,TPM,method="s")$estimate,
                                                               p.exp=cor.test(DP.exp,TPM,method="s")$p.value,
                                                               rho.tai=cor.test(DP.tai,TPM,method="s")$estimate,
                                                               p.tai=cor.test(DP.tai,TPM,method="s")$p.value) %>% as.data.frame() %>% arrange(type)
library(tidyr)
df_long <- tt %>%
  pivot_longer(cols = c(rho.exp, rho.tai, p.exp, p.tai),
               names_to = c(".value", "metric"),
               names_pattern = "(rho|p)\\.(.*)")

df_long$metric_label <- ifelse(df_long$metric == "exp", "Expression", "TAI")

df_long$mysample <- factor(df_long$mysample, levels = unique(df_long$mysample))
ggplot(df_long, aes(x = mysample, y = metric_label, fill = rho)) +
  geom_tile(color = "white") +                              
  geom_text(aes(label = ifelse(p < 0.001, "***", "")),     
            size = 3.5, na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Rho") +        
  labs(x = "", y = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),                    
        legend.position = "right")

#################################################
####RPF per gene
mydata <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sce_cds_tpm_quant.txt",sep = "\t",stringsAsFactors = F)
sample <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sample.inf.csv",header = F,stringsAsFactors = F)
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)

datares <- mclapply(mc.cores = 10,4:20,function(x){
  mysample <- colnames(mydata)[x]
  samptmp <- sample$V2[which(sample$V1 == mysample)]
  if(samptmp == "HEK293T"){
    tmpdata <- mydata[,c(1,3,x)]
  } else {tmpdata <- mydata[,c(1,2,x)]}
  names(tmpdata)[2:3] <- c("control","exp") 
  tmpdata1 <- merge(tmpdata,alldp,by.x="name",by.y="genename") %>% dplyr::mutate(deltaTPM=exp-control,FCTPM=exp/control)
  ##1)virus yij
  virustmpdata <- tmpdata1 %>% dplyr::filter(type=="viral")
  virusyij <- mclapply(mc.cores=5,1:nrow(virustmpdata),function(i){
    seq <- substring(as.character(virustmpdata$myseq[i]),seq(1,(nchar(virustmpdata$myseq[i])-2),by=3),seq(3,nchar(virustmpdata$myseq[i]),by=3));
    a <- as.data.frame(table(seq));
    names(a) <- c("codon","fre");
    a$codon <- as.vector(a$codon)
    a$gene <- virustmpdata$name[i]
    a$amino <- codon$aa[match(a$codon,codon$codon)]
    a$TPM <- virustmpdata$exp[i]
    a
  }) %>% rbind.fill() %>%
    group_by(codon,amino) %>% dplyr::summarise(freq=sum(fre*TPM)) %>% as.data.frame() %>% 
    dplyr::filter(amino != "*" & amino != "W" & amino != "M") %>%
    group_by(amino) %>% dplyr::mutate(nn=sum(freq)) %>% dplyr::mutate(yijvirus=freq/nn)  
  ##hum yij and virus yij oushijuli 
  humtmpdata <- tmpdata1 %>% dplyr::filter(type=="hum") 
  sijofhumgene <- mclapply(mc.cores=5,1:nrow(humtmpdata),function(j){
    seq <- substring(as.character(humtmpdata$myseq[j]),seq(1,(nchar(humtmpdata$myseq[j])-2),by=3),seq(3,nchar(humtmpdata$myseq[j]),by=3));
    a <- as.data.frame(table(seq));
    names(a) <- c("codon","fre");
    a$codon <- as.vector(a$codon)
    b <- virusyij
    b$frehum <- a$fre[match(b$codon,a$codon)]
    
    b$frehum[which(is.na(b$frehum))] <- 0
    b %>% group_by(amino) %>% dplyr::mutate(nhum=sum(frehum)) %>% group_by(amino,codon) %>% dplyr::mutate(yijhum=frehum/nhum) -> yijhum
    yijhum$yijhum[which(yijhum$yijhum=="NaN")] <- 0
    yijhum %>% as.data.frame() %>%
      dplyr::summarize(sijhum=(sum((yijhum-yijvirus)^2)^0.5)) %>% 
      cbind(humtmpdata[j,c(-7)])
  }) %>% rbind.fill()   
  sijofhumgene %>% cbind(data.frame(stringsAsFactors = F,mysample,samptmp))
}) %>% rbind.fill()

save(datares,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.Rdata")


########fig3 d and e
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.Rdata")
source("~/Rfunction/style.print.R")
datares %>%
  dplyr::filter(samptmp != "HEK293T" & mysample != "GSM8401429_NI812_trimmed_cds_tpm") %>%
  ggplot(aes(x=DP.exp,y=sijhum))+
  geom_point(shape = 1,size=0.5)+
  geom_smooth(method = "lm",se=F)+
  facet_grid(~mysample)+
  style.print()
datares %>%
  dplyr::filter(samptmp != "HEK293T" & mysample != "GSM8401429_NI812_trimmed_cds_tpm") %>%
  ggplot(aes(x=DP.tai,y=sijhum))+
  geom_point(shape = 1,size=0.5)+
  geom_smooth(method = "lm",se=F)+
  facet_grid(~mysample)+
  style.print()

datares %>% group_by(mysample,samptmp) %>% dplyr::summarize(rho.exp=cor.test(DP.exp,sijhum,method="s")$estimate,
                                                            p.exp=cor.test(DP.exp,sijhum,method="s")$p.value,
                                                            rho.tai=cor.test(DP.tai,sijhum,method="s")$estimate,
                                                            p.tai=cor.test(DP.tai,sijhum,method="s")$p.value) %>% arrange(samptmp)  


###plot fig2
library(ggExtra)
sample <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sample.inf.csv",header = F,stringsAsFactors = F)[3:19,] %>% arrange(desc(V2))
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.Rdata")
rnaThresQuant <- c(1,0.5,0.3,0.1,0.05,0.04,0.03,0.02,0.01)
sampletmp <- unique(datares$mysample)
i= "GSM8401429_NI812_trimmed_cds_tpm"
dfcor.exp.D <- mclapply(mc.cores = 8,sampletmp,function(i){
  datares1 <- datares %>% dplyr::filter(mysample == i)
  sampletype <- sample$V2[which(sample$V1==i)]
  myres <- mclapply(mc.cores = 9,rnaThresQuant,function(x){
    thisSubset1 <- (datares1 %>% arrange(desc(exp)))[1:floor(x*nrow(datares1)),] %>% dplyr::filter(!is.infinite(FCTPM))
    ##plot
    p1 <- thisSubset1 %>%
      ggplot(aes(x=sijhum,y=FCTPM))+
      geom_point(shape=1,size=0.2)+
      scale_y_continuous(trans = "log2")+
      geom_smooth(method = "lm",se=F)+
      labs(x="Dopt",y="Repression level")+style.print()
    ggMarginal(p1, margins = "y", type = "histogram", groupColour = F, groupFill = F)
    cor.test(thisSubset1$DP.exp,(thisSubset1$FCTPM),method = "s",alternative = "greater")
    cor.test(thisSubset1$DP.tai,(thisSubset1$FCTPM),method = "s",alternative = "greater")
    cor.test(thisSubset1$sijhum,(thisSubset1$FCTPM),method = "s",alternative = "greater")
    res <- lapply(1:1000, function(y){
      set.seed(y)
      thisSubset <- thisSubset1[sample(1:nrow(thisSubset1),nrow(thisSubset1),replace = T),]
      set.seed(y)
      randomSubset <- datares1[sample(1:nrow(datares1),nrow(thisSubset1),replace = T),] %>% dplyr::filter(!is.infinite(FCTPM))
      pcorObj1 <- cor.test(thisSubset$sijhum,thisSubset$FCTPM,method="s")
      pcorObj2 <- cor.test(thisSubset$DP.exp,thisSubset$FCTPM,method="s")
      pcorObj3 <- cor.test(thisSubset$DP.tai,thisSubset$FCTPM,method="s")
      
      pcorObj4 <- cor.test(randomSubset$sijhum,randomSubset$FCTPM,method="s")
      pcorObj5 <- cor.test(randomSubset$DP.exp,randomSubset$FCTPM,method="s")
      pcorObj6 <- cor.test(randomSubset$DP.tai,randomSubset$FCTPM,method="s")
      
      data.frame(stringsAsFactors = F,
                 thres = rep(paste(x*100,"%",sep = ""),6),
                 estimate = c(pcorObj1$estimate,pcorObj2$estimate,pcorObj3$estimate,pcorObj4$estimate,pcorObj5$estimate,pcorObj6$estimate),
                 p.val = c(pcorObj1$p.value,pcorObj2$p.value,pcorObj3$p.value,pcorObj4$p.value,pcorObj5$p.value,pcorObj6$p.value),
                 
                 len = rep(nrow(thisSubset),6),
                 corWith = c("real","real","real","random","random","random"),
                 type2 =c("sij","exp","tai","sij","exp","tai"))
    }) %>% rbind.fill()
    res %>% group_by(thres,corWith,len,type2) %>% dplyr::summarize(mrho=mean(estimate),sd=sd(estimate)) %>% 
      merge(data.frame(stringsAsFactors = F,
                       pvalue=c(wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="exp"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="exp"))$estimate,alternative = "greater")$p.value,
                                wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="sij"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="sij"))$estimate,alternative = "greater")$p.value,
                                wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="tai"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="tai"))$estimate,alternative = "greater")$p.value),
                       type2=c("exp","sij","tai")),by.x="type2",by.y="type2")
  }) %>% rbind.fill() %>% cbind(data.frame(stringsAsFactors = F,sampledd = i,sampletype))
}) %>% rbind.fill()



dfcor.exp.D$thres <- factor(dfcor.exp.D$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.exp.D$corWith <- factor(dfcor.exp.D$corWith,levels = c("real","random"))
dfcor.exp.D$sampledd <- factor(dfcor.exp.D$sampledd,levels = sample$V1)

save(dfcor.exp.D,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/plot.trans.Rdata")

load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/plot.trans.Rdata")
source("~/Rfunction/style.print.R")
dfcor.exp.D %>%
  dplyr::filter(sampletype != "HEK293T" & sampledd == "GSM8401429_NI812_trimmed_cds_tpm" & type2=="exp") %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (FCTPM ~ metric)",limits = c(-0.15,0.6),breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes (top expression percentile)") +
  facet_grid(sampledd~type2,
             labeller = labeller(sampledd = function(x) substr(x, 1, 16)))+
  #geom_hline(yintercept = c(0.1, 0.2, .3, .4, .5))+
  style.print()

dfcor.exp.D$type2 <- factor(dfcor.exp.D$type2,levels = c("exp","tai","sij"))
dfcor.exp.D %>%
  dplyr::filter(sampletype != "HEK293T" & sampledd != "GSM8401429_NI812_trimmed_cds_tpm") %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (FCTPM ~ metric)",limits = c(-0.2,0.4),breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes (top expression percentile)") +
  facet_grid(sampledd~type2,
             labeller = labeller(sampledd = function(x) substr(x, 1, 16)))+
  #geom_hline(yintercept = c(0.1, 0.2, .3, .4, .5))+
  style.print()


dfcor.exp.D %>%
  dplyr::filter(sampletype != "HEK293T" & sampledd != "GSM8401429_NI812_trimmed_cds_tpm" & pvalue>0.05) %>% arrange(type2)
dfcor.exp.D %>%
  dplyr::filter(sampletype != "HEK293T" & pvalue<0.05 & pvalue>0.01) %>% arrange(type2)
dfcor.exp.D %>%
  dplyr::filter(sampletype != "HEK293T" & pvalue<0.01 & pvalue>0.001) %>% arrange(type2)


dfcor.exp.D %>%
  dplyr::filter(sampletype == "HEK293T") %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (FCTPM ~ metric)",limits = c(-0.5,0.2),breaks = c(-0.4,-0.2,0,0.2,0.4)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes (top expression percentile)") +
  facet_grid(type2~sampledd,
             labeller = labeller(sampledd = function(x) substr(x, 1, 16)))+
  style.print()





