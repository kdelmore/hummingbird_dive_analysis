## example script used to generate bootstrapped values for bchu
## ran this on hpc and 10 versions of it to save time

library(lme4)
library(data.table)
library(qgraph) 
library(RColorBrewer)
library(dplyr)

## bring in correlation matrix from phenotypic_network_1.R
cor.orig <- read.csv("bchu_cov1_rand_qnorm_july_24_2023_add_dopp.csv") ## generated in phenotypic_network_1.R
#cor.orig <- read.csv("bchu_cov1_rand_apr6_2023.csv") ## generated in phenotypic_network_1.R
#cor.orig <- read.csv("hyb_cov1_rand_apr6_2023.csv") ## generated in phenotypic_network_1.R
cor.orig <- as.matrix(cor.orig)

## bring in raw data for bootstrapping
dat0 <- read.csv("dive_kinematics_july_04_2023.csv")
dat0$Male_ID <- as.factor(dat0$Male_ID)
dat0 <- dat0[,-c(3)] ## get rid of dive number
dat0 <- dat0[,-18] ## get rid of Dopp.1
dat0 <- dat0 %>% select(-contains("30frames"))
#dat0 <- dat0 %>% select(-contains("Doppler"))
dat0$X_length <- dat0$XEnd_Relative_to_z.0 - dat0$XStart_Relative_to_z.0
dat0$W1_length <- dat0$W1End_Relative_to_z.0-dat0$W1Start_Relative_to_z.0
dat0$Y_length <- dat0$YEnd_Relative_to_z.0-dat0$YStart_Relative_to_z.0
dat0$W2_length <- dat0$W2End_Relative_to_z.0-dat0$W2Start_Relative_to_z.0
dat0 <- dat0 %>% select(-contains("End"))

## choose taxon and have to remove some for rubies
dat0 <- subset(dat0,dat0$Taxon=="bchu") 

## get confidence intervals for bootstraps of correlation matrix
t=1000 # 100000 is what they used; number of iterations for bootstrapping, play with the number of iterations to make sure the network topology stabilizes (for our paper we used 100k to be certain of a single topology)
diag(cor.orig)=0
boot.corr=array(dim=c(nrow(cor.orig),ncol(cor.orig),t))
cors=vector(length=t)

  # loop to create t=iterations datasets composed of k trait observations from n individuals resampled with replacement from our original dataset
  for (k in 1:t)
  {
    s = sample(1:nrow(dat0),nrow(dat0),replace=T)
    new_dat0 = dat0[s,]
    mat2<-matrix(0,(dim(new_dat0)[2]-2),(dim(new_dat0)[2]-2))
    for(i in 3:dim(new_dat0)[2]){
      for(j in 3:dim(new_dat0)[2]){
        dat<-na.omit(new_dat0[,c(2,i,j)])
        mi<-predict(lmer(dat[,2]~(1|dat[,1]),REML=F))
        mj<-predict(lmer(dat[,3]~(1|dat[,1]),REML=F))
        
        mi2 <- cbind(dat[,c(1,2)], mi)
        mj2 <- cbind(dat[,c(1,3)], mj)
        
        sxx<-t(mi2[,3]-mi2[,2])%*%(mi2[,3]-mi2[,2])
        syy<-t(mj2[,3]-mj2[,2])%*%(mj2[,3]-mj2[,2])
        sxy<-t(mi2[,3]-mi2[,2])%*%(mj2[,3]-mj2[,2])
        #mat2[i,j]<-sxy/sqrt(sxx)/sqrt(syy)
        mat2[(i-2),(j-2)]<-sxy/sqrt(sxx)/sqrt(syy)
      }}
    colnames(mat2) <- colnames(new_dat0[,3:length(new_dat0)])
    rownames(mat2) <- colnames(new_dat0[,3:length(new_dat0)])
    
    boot.corr[,,k]=mat2
    cors[k]=cor.test(mat2[which(upper.tri(mat2)==TRUE)],cor.orig[which(upper.tri(cor.orig)==TRUE)])$estimate
  }
  
avg.cor=apply(boot.corr,c(1,2),mean,na.rm=T) # calculate average correlation matrix from 1000 resampling iterations
#calculate bootstrap confidence intervals for each edge in the original dataset
lower.cor=apply(boot.corr, c(1,2), quantile, probs=c(0.025))
upper.cor=apply(boot.corr, c(1,2), quantile, probs=c(0.975))
boot.cor.output = list(orig.cor=cor.orig,avg.cor=avg.cor,lower.cor=lower.cor,upper.cor=upper.cor)

saveRDS(boot.cor.output, file="1_boot.cor.output_bchu_1.RData")
saveRDS(boot.corr,file="1_boot.corr_bchu_1.RData")
