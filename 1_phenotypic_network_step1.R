## step 1 = adjust for individual id so we can use dive as our unit of analysis

library(dplyr)
library(lme4)
quantNorm = function(x){qnorm(rank(x, ties.method = "average")/(length(x)+1))}

## load data, remove some variables and estimate others
dat0 <- read.csv("C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis_july24/dive_kinematics_july_04_2023.csv")
dat0$Male_ID <- as.factor(dat0$Male_ID)
dat0 <- dat0[,-c(3)] ## get rid of dive number
dat0 <- dat0[,-18] ## get rid of Dopp.1
dat0$X_length <- dat0$XEnd_Relative_to_z.0 - dat0$XStart_Relative_to_z.0
dat0$W1_length <- dat0$W1End_Relative_to_z.0-dat0$W1Start_Relative_to_z.0
dat0$Y_length <- dat0$YEnd_Relative_to_z.0-dat0$YStart_Relative_to_z.0
dat0$W2_length <- dat0$W2End_Relative_to_z.0-dat0$W2Start_Relative_to_z.0
dat0$X_length <- dat0$XEnd_Relative_to_z.0 - dat0$XStart_Relative_to_z.0
dat0 <- dat0 %>% select(-contains("End"))

## choose taxon
dat0 <- subset(dat0,dat0$Taxon=="rthu") ## 45 dives and 32 variables
## or
dat0 <- subset(dat0,dat0$Taxon=="bchu") ## 42 dives and 32 variables

## adjust mean by predicted value of birdID as random effect
temp <- apply(dat0[,3:ncol(dat0)],2,quantNorm)
dat0 <- data.frame(cbind(dat0[,1:2],temp))

mat2<-matrix(0,(dim(dat0)[2]-2),(dim(dat0)[2]-2))
for(i in 3:dim(dat0)[2]){
for(j in 3:dim(dat0)[2]){
dat<-na.omit(dat0[,c(2,i,j)])
mi<-predict(lmer(dat[,2]~(1|dat[,1]),REML=F))
mj<-predict(lmer(dat[,3]~(1|dat[,1]),REML=F))

mi2 <- cbind(dat[,c(1,2)], mi)
mj2 <- cbind(dat[,c(1,3)], mj)

sxx<-t(mi2[,3]-mi2[,2])%*%(mi2[,3]-mi2[,2])
syy<-t(mj2[,3]-mj2[,2])%*%(mj2[,3]-mj2[,2])
sxy<-t(mi2[,3]-mi2[,2])%*%(mj2[,3]-mj2[,2])
mat2[(i-2),(j-2)]<-sxy/sqrt(sxx)/sqrt(syy)
}}
colnames(mat2) <- colnames(dat0[,3:length(dat0)])
rownames(mat2) <- colnames(dat0[,3:length(dat0)])

write.table(mat2,"/output/bchu_cov1_rand_qnorm_july_24_2023.csv",sep=",",row.names=F,quote=F)
write.table(mat2,"/output/rthu_cov1_rand_qnorm_july_24_2023.csv",sep=",",row.names=F,quote=F)