library(circular)

## load data
angles <- read.csv("C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis/female_head_angle_24062023.csv")

## when working with bchu
temp <- subset(angles,angles$t_relative_z>-0.49)
angles <- subset(temp,temp$t_relative_z<0)

## when working with rthu
temp <- subset(angles,angles$t_relative_z>-0.24)
angles <- subset(temp,temp$t_relative_z<0)

## convert mean angle to circular
angles$X360angle <- ifelse(angles$theta<0,angles$theta+360,angles$theta)
angles$circ <- as.circular(angles$X360angle,type=c("angles"),units=c("degrees"),zero=pi/2,rotation=c("clock"))

## subset dive data and subspecies
angles$ind_dive <- paste(angles$ind,angles$dive,sep="_")
dive <- aggregate(circ ~ ind_dive+ind+dive+taxon, angles, FUN=median.circular)
plot.circular(dive$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1)
rthu_dive <- subset(dive,taxon=="RBHU")
bchu_dive <- subset(dive,taxon=="BCHU")
plot.circular(rthu_dive$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1,col="red")
plot.circular(bchu_dive$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1,col="purple")

## summarize by individual
ind <- aggregate(circ ~ ind +taxon, dive, FUN=mean.circular)
rthu_ind <- subset(ind,taxon=="RBHU")
bchu_ind <- subset(ind,taxon=="BCHU")
plot.circular(rthu_ind$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1,col="red")
plot.circular(bchu_ind$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1,col="purple")

## check if consistent by individual
inds <- unique(dive$ind)

output<-NA

for (i in 1:length(inds)){
  
  dat.i <- subset(dive,ind==inds[i])
  n_dives <- nrow(dat.i)
  plot.circular(dat.i$circ, pch=16, cex=1, stack=TRUE,sep=0.05,shrink=1,main=inds[i])
  mean.i <- mean.circular(dat.i$circ)
  ray_value_general <- rayleigh.test(dat.i$circ)[1]
  p_val_general <- rayleigh.test(dat.i$circ)[2]
  output.i <- cbind(inds[i],n_dives,mean.i,ray_value_general,p_val_general,p_val_zero,dat.i[1,4])
  output <- rbind(output,output.i)

}

output <- as.data.frame(output)
output <- output[-1,]

## rayleigh test on all dives or after averaging dives for each individual
rayleigh.test(rthu_dive$circ)
rayleigh.test(bchu_dive$circ)
rayleigh.test(rthu_ind$circ)
rayleigh.test(bchu_ind$circ)

## plot
plot(rthu_dive$circ, stack=TRUE,sep=0.05,shrink=1)
rose.diag(rthu_dive$circ,bins=18,add=T)
plot(bchu_dive$circ,  stack=TRUE,sep=0.05,shrink=1)
rose.diag(bchu_dive$circ,bins=18,add=T)

plot(rthu_ind$circ, stack=TRUE,sep=0.05,shrink=1)
rose.diag(rthu_ind$circ,bins=18,add=T)
plot(bchu_ind$circ,  stack=TRUE,sep=0.05,shrink=1)
rose.diag(bchu_ind$circ,bins=18,add=T)