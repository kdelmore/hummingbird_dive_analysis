library(circular)

## load data
angles <- read.csv("./input/male_head_angles_dives_24112024.csv")
angles$linear_diff <- abs(angles$AngletoFem-angles$AngletoSun)

## has to be in front
temp <- subset(angles,angles$AngletoSun<90)
infront <- subset(temp,temp$AngletoSun>-90)

## subset by subspecies
bchu_linear_diff <- subset(ind_linear_diff,ind_linear_diff$Taxon=="bchu")
rthu_linear_diff <- subset(ind_linear_diff,ind_linear_diff$Taxon=="rthu")

## run test
ks.test(bchu_linear_diff$linear_diff,"punif",0,180)
ks.test(rthu_linear_diff$linear_diff,"punif",0,180)

## plot
hist(ruby$linear_diff,xlim=c(0,180),ylim=c(0,5),col="black",xaxt='n',xlab="Orientation towards sun")
axis(side=1, at=seq(0,180,45), labels=c("0","45","90","135","180"))
hist(rthu_linear_diff$linear_diff,xlim=c(0,180),col="red",add=T,breaks=6)
hist(black$linear_diff,xlim=c(0,180),ylim=c(0,5),col="black",breaks=17,xaxt='n',xlab="Orientation towards sun")
axis(side=1, at=seq(0,180,45), labels=c("0","45","90","135","180"))
hist(bchu_linear_diff$linear_diff,xlim=c(0,180),col="purple",add=T,breaks=6)
