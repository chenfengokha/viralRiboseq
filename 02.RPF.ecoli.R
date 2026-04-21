library(Biostrings)
library(plyr)
library(dplyr)
library(parallel)
library(ggplot2)
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.Ecoli.Rdata")
xijofD <- xijofD[,c(1:3,5)]
names(xijofD)[4] <- "xijoftai"
##############D
GmRcds <- read.csv("/home/chenfeng/project/HGT/exp/riboseq/1.reference/norm/GmR.csv",stringsAsFactors = F,header = T)
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)[-60,]
GmRyij <- lapply(1:4, function(x){
  id <- GmRcds$id[x]
  heter <- GmRcds$heter[x]
  tmpseq <- GmRcds$sequence[x]
  
  mydf0 <- substring(as.character(tmpseq),seq(1,(nchar(tmpseq)-2),by=3),seq(3,nchar(tmpseq),by=3))
  mydf <- mydf0 %>% table() %>% as.data.frame()
  names(mydf)[1] <- c("codon")
  mydf$codon <- as.vector(mydf$codon)
  
  condontmp <- codon
  condontmp$yfre <- mydf$Freq[match(condontmp$codon,mydf$codon)]
  condontmp$yfre[which(is.na(condontmp$yfre))] <- 0
  condontmp %>% group_by(aa) %>% dplyr::mutate(ny=sum(yfre)) %>% group_by(aa,codon) %>% dplyr::mutate(yij=yfre/ny) -> condontmp
  condontmp$yij[which(condontmp$yij=="NaN")] <- 0
  condontmp %>% cbind(data.frame(stringsAsFactors = F,id,heter))
}) %>% rbind.fill()

xijofD$yijGmR_0.127_2.1 <- (GmRyij %>% dplyr::filter(heter == "GmR_0.127_2.1"))$yij[match(xijofD$codon,(GmRyij %>% dplyr::filter(heter == "GmR_0.127_2.1"))$codon)]
xijofD$yijGmR_0.278_1.1 <- (GmRyij %>% dplyr::filter(heter == "GmR_0.278_1.1"))$yij[match(xijofD$codon,(GmRyij %>% dplyr::filter(heter == "GmR_0.278_1.1"))$codon)]
xijofD$yijGmR_0.778_2.1 <- (GmRyij %>% dplyr::filter(heter == "GmR_0.778_2.1"))$yij[match(xijofD$codon,(GmRyij %>% dplyr::filter(heter == "GmR_0.778_2.1"))$codon)]
xijofD$yijGmR_0.878_2.1 <- (GmRyij %>% dplyr::filter(heter == "GmR_0.878_2.1"))$yij[match(xijofD$codon,(GmRyij %>% dplyr::filter(heter == "GmR_0.878_2.1"))$codon)]
allcds <- readDNAStringSet("/home/chenfeng/project/HGT/exp/riboseq/1.reference/norm/all.norm.cds.fa")

all.DP.sij <- mclapply(mc.cores = 10,1:length(allcds),function(x){
  tmpseq <- as.character(allcds[[x]])
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
    cbind(data.frame(stringsAsFactors = F, 
                     sij_0.127=sum((xij$yijGmR_0.127_2.1-xij$yij)^2)^0.5,
                     sij_0.278=sum((xij$yijGmR_0.278_1.1-xij$yij)^2)^0.5,
                     sij_0.778=sum((xij$yijGmR_0.778_2.1-xij$yij)^2)^0.5,
                     sij_0.878=sum((xij$yijGmR_0.878_2.1-xij$yij)^2)^0.5)) %>%
    cbind(data.frame(stringsAsFactors = F,genename=names(allcds)[x]))
  
}) %>% rbind.fill()

save(all.DP.sij,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/all.DP.sij.ecoli.Rdata")

cor.test(all.DP.sij$DP.exp,all.DP.sij$sij_0.127,method = "s")
cor.test(all.DP.sij$DP.exp,all.DP.sij$sij_0.278,method = "s")
cor.test(all.DP.sij$DP.exp,all.DP.sij$sij_0.778,method = "s")
cor.test(all.DP.sij$DP.exp,all.DP.sij$sij_0.878,method = "s")


#################################################
####RPF per gene

control.tpm <- read.csv("/home/chenfeng/project/HGT/exp/riboseq/4.ribo-seq/5.riboparser/10.quantification/WT_cds_tpm_quant.txt",sep = "\t",stringsAsFactors = F)
myfile <- system("ls /home/chenfeng/project/HGT/exp/riboseq/4.ribo-seq/5.riboparser/10.quantification/sce*cds_tpm_quant.txt",intern=T)
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/all.DP.sij.ecoli.Rdata")


datares <- mclapply(mc.cores = 4,1:length(myfile),function(x){
  mysample <- strsplit(strsplit((myfile)[x],"/sce.")[[1]][2],"_cds_tpm_quant.txt")[[1]][1]
  
  expriment.tpm <- read.csv(myfile[x],sep = "\t",stringsAsFactors = F)
  names(expriment.tpm)[2] <- "TPM.experiment"
  expriment.tpm$TPM.control <- control.tpm$WT.fastq.gz_cds_tpm[match(expriment.tpm$name,control.tpm$name)]
  
  expriment.tpm %>% merge(all.DP.sij,by.x="name",by.y="genename") %>%
    dplyr::mutate(FCTPM=TPM.experiment/TPM.control) %>% 
    cbind(data.frame(stringsAsFactors = F,mysample))
}) %>% rbind.fill()

save(datares,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.ecoli.TPMFC.Rdata")
##############
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.ecoli.TPMFC.Rdata")
source("~/Rfunction/style.print.R")
tmpdd <- datares %>%
  dplyr::filter(mysample == "44_2") 
tmpdd %>%
  ggplot(aes(x=DP.exp,y=sij_0.878))+
  geom_point(shape = 1,size=0.5)+
  geom_smooth(method = "lm",se=F)+
  facet_grid(~mysample)+
  style.print()
cor.test(tmpdd$DP.exp,tmpdd$sij_0.878,method = "s")$p.value


###plot
#sample <- read.csv("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/4.ribo-seq/5.riboparser/10.quantification/sample.inf.csv",header = F,stringsAsFactors = F)[3:19,] %>% arrange(desc(V2))
load("~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/alldp.sij.ecoli.TPMFC.Rdata")
rnaThresQuant <- c(1,0.5,0.3,0.1,0.05,0.03)

sampleD <- data.frame(stringsAsFactors = F,id=c("14_3","4_2","38_1","44_2"),D=c("sij_0.127","sij_0.278","sij_0.778","sij_0.878"))
sampletmp <- unique(sampleD$id)
i = sampletmp[2]
dfcor.exp.D <- mclapply(mc.cores = 4,sampletmp,function(i){
  datares1 <- datares %>% dplyr::filter(mysample == i) %>% dplyr::filter(!(substr(name,1,3) %in% c("YFP")))
  sampletype <- sampleD$D[which(sampleD$id==i)]
  ranktmp <- which(colnames(datares1) == sampletype)
  myres <- mclapply(mc.cores = 9,rnaThresQuant,function(x){
    thisSubset1 <- (datares1 %>% arrange(desc(TPM.experiment)))[1:floor(x*nrow(datares1)),] %>% dplyr::filter(!is.infinite(FCTPM))
    ##plot
    # p1 <- thisSubset1 %>%
    #   ggplot(aes(x=DP.tai,y=FCTPM))+
    #   geom_point(shape=1,size=0.2)+
    #   scale_y_continuous(trans = "log2")+
    #   geom_smooth(method = "lm",se=F)+
    #   labs(x="D",y="Repression level")+style.print()
    # ggMarginal(p1, margins = "y", type = "histogram", groupColour = F, groupFill = F)
    # cor.test(thisSubset1$DP.exp,(thisSubset1$FCTPM),method = "s",alternative = "less")
    # cor.test(thisSubset1$DP.tai,(thisSubset1$FCTPM),method = "s",alternative = "less")
    # cor.test(thisSubset1$sij_0.878,(thisSubset1$FCTPM),method = "s",alternative = "greater")
    res <- lapply(1:1000, function(y){
      set.seed(y)
      thisSubset <- thisSubset1[sample(1:nrow(thisSubset1),nrow(thisSubset1),replace = T),]
      set.seed(y)
      randomSubset <- datares1[sample(1:nrow(datares1),nrow(thisSubset1),replace = T),] %>% dplyr::filter(!is.infinite(FCTPM))
      pcorObj1 <- cor.test(thisSubset[,ranktmp],thisSubset$FCTPM,method="s")
      pcorObj2 <- cor.test(thisSubset$DP.exp,thisSubset$FCTPM,method="s")
      pcorObj3 <- cor.test(thisSubset$DP.tai,thisSubset$FCTPM,method="s")
      
      pcorObj4 <- cor.test(randomSubset[,ranktmp],randomSubset$FCTPM,method="s")
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
                       pvalue=c(wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="exp"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="exp"))$estimate,alternative = "less")$p.value,
                                wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="sij"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="sij"))$estimate,alternative = "greater")$p.value,
                                wilcox.test((res %>% dplyr::filter(corWith == "real" & type2=="tai"))$estimate,(res %>% dplyr::filter(corWith == "random" & type2=="tai"))$estimate,alternative = "less")$p.value),
                       type2=c("exp","sij","tai")),by.x="type2",by.y="type2")
  }) %>% rbind.fill() %>% cbind(data.frame(stringsAsFactors = F,sampledd = i,sampletype))
}) %>% rbind.fill()


source("~/Rfunction/style.print.R")
dfcor.exp.D$thres <- factor(dfcor.exp.D$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.exp.D$corWith <- factor(dfcor.exp.D$corWith,levels = c("real","random"))
dfcor.exp.D$sampletype <- factor(dfcor.exp.D$sampletype,levels = c("sij_0.127","sij_0.278","sij_0.778","sij_0.878"))

#save(dfcor.exp.D,file = "~/project/chenfengdata5/human.ribo/riboseq/03.riboseqhumtissue/0.R/00.Rdata/plot.trans.ecoli.Rdata")
dfcor.exp.D$type2 <- factor(dfcor.exp.D$type2,levels = c("exp","tai","sij"))

##"sij_0.127","sij_0.278","sij_0.778","sij_0.878"
dfcor.exp.D %>%
  #dplyr::filter(sampletype=="sij_0.127") %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (FCTPM ~ metric)",limits = c(-0.3,0.7),breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes (top expression percentile)") +
  facet_grid(sampletype~type2,
             labeller = labeller(sampletype = function(x) substr(x, 5, 9)))+
  #geom_hline(yintercept = c(0.1, 0.2, .3, .4, .5))+
  style.print()





