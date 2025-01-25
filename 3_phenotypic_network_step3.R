## step 3 = use bootraps to generate final figures
## modified from wilkins 2015 https://doi.org/10.1098/rspb.2015.1574

library(lme4)
library(data.table)
library(qgraph) 
library(RColorBrewer)
library(dplyr)
library(abind)

## bring in correlation matrix from phenotypic_network_1.R
cor.orig <- read.csv("C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis_july24/bchu_cov1_rand_qnorm_july_24_2023_add_dopp.csv") ## generated in phenotypic_network_1.R
cor.orig <- read.csv("C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis_july24/rthu_cov1_rand_qnorm_july_24_2023_add_dopp_w1.csv") ## generated in phenotypic_network_1.R

cor.orig <- as.matrix(cor.orig)

## bring in boot.corr.output from phenotypic_network_2.R
files = list.files (pattern="*boot.corr_bchu.RData")
files = list.files (pattern="*boot.corr_rthu.RData")

for (i in seq_along(files)) {
  if (i == 1){
    boot.corr <- readRDS(files[i])
  } else {
    dat <- readRDS(files[i])
    boot.corr <- abind(boot.corr,dat)
  }
}

avg.cor=apply(boot.corr,c(1,2),mean,na.rm=T) # calculate average correlation matrix from 1000 resampling iterations
lower.cor=apply(boot.corr, c(1,2), quantile, probs=c(0.025)) #calculate bootstrap confidence intervals for each edge in the original dataset
upper.cor=apply(boot.corr, c(1,2), quantile, probs=c(0.975))
boot.cor.output = list(orig.cor=cor.orig,avg.cor=avg.cor,lower.cor=lower.cor,upper.cor=upper.cor)

## plot bootstrap correlation CIs vs empirical correlations, and extract robust edges
par(cex.lab=1.8,mar=c(5,5,3,3),cex.axis=1.4,font.lab=2,mgp=c(3.5,1,0))  
D<-boot.cor.output
emp.edges=D$orig.cor[which(upper.tri(D$orig.cor)==TRUE)]
boot.avgs=D$avg.cor[which(upper.tri(D$avg.cor)==TRUE)]
boot.lower=D$lower.cor[which(upper.tri(D$lower.cor)==TRUE)]
boot.upper=D$upper.cor[which(upper.tri(D$upper.cor)==TRUE)]
ci.0=c("black","red")[((boot.lower>0&boot.upper>0)|(boot.lower<0&boot.upper<0))+1]
plot(emp.edges,boot.avgs, pch=19, ylim=c(-1,1), xlim=c(-1,1), xlab="Empirical Correlations", ylab="Bootstrap Correlations",col=ci.0,las=1)
for (i in 1:length(emp.edges)){
  lines(c(emp.edges[i], emp.edges[i]), c(boot.lower[i],boot.upper[i]),col=ci.0[i])
}
min.threshold<-
abline(v=0.3,lty=2)
abline(v=-0.3,lty=2)
abline(h=0, lty=1)

# find threshold where bootstrap CIs don't overlap 0.
mini=maxi=0
for(i in 1: length(emp.edges))
{
  if(emp.edges[i]<mini&ci.0[i]=="black"){mini=emp.edges[i]
  }else if(emp.edges[i]>maxi&ci.0[i]=="black"){maxi=emp.edges[i]
  }else{}
}

ci.0_full=c("black","red")[((D$lower.cor>0&D$upper.cor>0)|(D$lower.cor<0&D$upper.cor<0))+1]  
robust.indx<-which(ci.0_full=="red")
nonrobust.indx<-which(ci.0_full=="black")
boot.analy <- list(thresh.min=mini,thresh.max=maxi,robust.indx=robust.indx,nonrobust.indx=nonrobust.indx)

## remove nonrobust edges from graph 
nodes <- cor.orig
#nodes [nodes <0.75] <- 0
nodes[boot.analy$nonrobust.indx]<-0
diag(nodes)<-0

### TO MAKE BCHU NETWORK THE SAME AS RTHU NETWORK
#nodes <- nodes[-c(27:31),-c(27:31)]

## making phenotype network
## set node shapes (modalities) 
## can only use these three shapes and "diamonds"
#shps<-rep("circle",nrow(nodes))
#clrs<-rep("gray",nrow(nodes))
#nombres=colnames(nodes)

## BCHU
## circle = timing, triangle = not
shps=c("circle","circle","triangle","circle","triangle","circle","triangle","circle","circle","triangle",
       "triangle","triangle","triangle","triangle",
       "circle",
       "circle","triangle","triangle","circle",
       "circle","triangle","triangle","circle","triangle","triangle","triangle",
       "circle","triangle","triangle","triangle","circle","circle",
       "circle","circle")
## fed = kinematics, 3b5 = female, 5ec = sond
clrs=c(rep("#fde725",8),rep("#3b528b",7),rep("#5ec962",19))
nombres <- c("DiveLength", "DiveStart", "MaxVel", "MaxVelTime", "MaxXVel", "MaxXVelTime", "MaxYVel", "MaxYVelTime","ClosestFemaleTime", "ClosestFemaleDist",
"MaxDopp","MaxDoppTime","MinDopp","MinDoppTime", ##new
"DoppTime",
"W1Start","W1LowFreq","W1HighFreq","GapW1W2", ##new
"W2Start", "W2LowFreq", "W2HighFreq", "YStart", "Y#Elements","YLowFreq","YHighFreq", 
"XStart", "X#Elements","XLowFreq","XHighFreq","XLength","W1Length", ##new
"YLength", "W2Length")

## RTHU
## circle = timing, triangle = not
shps=c("circle","circle","triangle","circle","triangle","circle","triangle","circle","circle","triangle",
       "triangle","triangle","triangle","triangle",
       "circle",
       "circle","triangle","triangle","circle",
       "circle","triangle","triangle","circle","triangle","triangle","triangle",
       "circle",
       "circle","circle")
## fed = kinematics, 3b5 = female, 5ec = sond
clrs=c(rep("#fde725",8),rep("#3b528b",7),rep("#5ec962",14))
nombres <- c("DiveLength", "DiveStart", "MaxVel", "MaxVelTime", "MaxXVel", "MaxXVelTime", "MaxYVel", "MaxYVelTime","ClosestFemaleTime", "ClosestFemaleDist",
             "MaxDopp","MaxDoppTime","MinDopp","MinDoppTime", ##new
             "DoppTime",
             "W1Start","W1LowFreq","W1HighFreq","GapW1W2", ##new
             "W2Start", "W2LowFreq", "W2HighFreq", "YStart", "Y#Elements","YLowFreq","YHighFreq", 
             "W1Length", ##new
             "YLength", "W2Length")

# SOME of the graphical parameters we can specify (but see qgraph help for more)
min.edge=0.3 # minimum correlation to display
lay<-"spring" #specify arrangement of nodes; here, spring gives a force-embedded layout; can also specify "circle", "groups", or positions from another qgraph object (i.e. from Object$layout) 
lay.params<-list(max.delta=25) # Tweak the spring algorithm slightly for aesthetics (try 26 vs NULL[the default setting] and look at placement of RHue and BHue)
# max.delta defaults to the number of nodes and describes maximum displacement in the graph. Default settings in our case made it difficult to see the relationship
# between hue and other measures of color, which is of interest in our study ( see more spring algorithm parameters at > ?qgraph.layout.fruchtermanreingold )
nscale=6 #(specify scalar size of nodes; you'll want to change this depending on the size of your monitor/figure output)
lab.length<-"0000" #Normalize scaling of node labels to a string of length 5
fade.TF<-F #Should edges fade as a function of correlation strength?
neg.color<- 1 #What color should negative edges be? (here, black)
pos.color<- 1 #What color should positive edges be? (here, black; default green)

Q<-qgraph( nodes,color=clrs,layout=lay,layout.par=lay.params,labels=nombres,
           minimum=min.edge,vsize=nscale,vsize=nscale,label.norm=lab.length,
           fade=fade.TF, shape=shps,label.scale=T,negCol=neg.color,posCol=pos.color,label.color=1,cut=0,negDashed=T)

savefile<-"C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis_july24/bchu_cov1_rand_qnorm_july_24_wdopp.pdf"
savefile<-"C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/analysis_july24/rthu_cov1_rand_qnorm_july_24_wdopp_w1.pdf"

pdf(savefile,width=4.5,height=4.5)
plot(Q)
dev.off()

#__________________________________
# 5) Calculating assortativity (modularity)
#----------------------------------

require(assortnet)
# diag(cor.orig)=0 #no self-loops (i.e. set diagonal correlations of 1 to 0, as this will bias measures of assortativity)

# switch shapes and colors for this analysis because it was originally written assuming shapes were the different modes but i switched above
shps = clrs

# Empirical assortatitivity. Remember to take the absolute value of the correlation matrix, since negative values will mess up estimates of assortment.
as.filt=assortment.discrete(abs(nodes),shps,SE=TRUE)
as.filt #assortativity of filtered data

#   * This is used to calculate the random expected assortativity

assort.nodeperm.orig=vector(length=1000)
assort.nodeperm.filt=vector(length=1000)
for (j in 1:1000){
  rand.shps=sample(shps,length(shps),replace=FALSE) #randomize the shape of each node while keeping the total number of each shape the same.
  assort.nodeperm.orig[j]=assortment.discrete(abs(cor.orig),rand.shps)$r #calculate assortativity for original data. Remember to use absolute values of correlations.
  assort.nodeperm.filt[j]=assortment.discrete(abs(nodes),rand.shps)$r #calculate assortativity for filtered correlations (>=|0.3|)
}
assort.nodeperm.filt

#P-value as the number of times node-permutation lead to observed level of assortment.
p.assort=length(which(assort.nodeperm.filt>=as.filt$r))/1000
p.assort # p-value (<<.001)

#####
# Plot assortativity measures of node permutation and empirical correlation

par(mfrow=c(1,1),font=2,mar=c(5,5,2,0))
hist(assort.nodeperm.filt,xlab="Assortativity",main="")
points(as.filt$r,8,pch=25,col=1,bg=1)

#   * Here, we can clearly see what the p.assort value tells us--that our value of assortativity is far outside the distribution of values predicted from random node permutation.
#     That is, in our dataset, trait associations by modality are nonrandom, but imperfect (i.e. there are many correlations within and some correlations across modalities)
# not for bchu though!

#__________________________________
# 6) Calculating redundancy metrics
#----------------------------------
# 
### Average correlation for full, filtered network

(avg.cor.full<-mean(abs(nodes[which(upper.tri(nodes)&nodes!=0)])))

### Edge density for full, filtered network
(edge.density.full<-length(nodes[upper.tri(nodes)&nodes!=0]) /length(which(upper.tri(nodes))) )

## add jaccards
## number of edges contained in both networks, divided by all edges contained in either or both networks
bchu_tri <- as.data.frame(abs(bchu_nodes[upper.tri(bchu_nodes)]))
rthu_tri <- as.data.frame(abs(rthu_nodes[upper.tri(rthu_nodes)]))
bchu_rthu_tri <- cbind(bchu_tri,rthu_tri)
names(bchu_rthu_tri) <- c("bchu","rthu")
## top (number in both networks)
bchu_rthu_tri$add <- apply(bchu_rthu_tri, 1, function(i) sum(i > 0))
top=sum(bchu_rthu_tri$add>=2)
## bottom (number in either network)
bottom=length(bchu_nodes[upper.tri(bchu_nodes)&bchu_nodes!=0])+length(rthu_nodes[upper.tri(rthu_nodes)&rthu_nodes!=0])-sum(bchu_rthu_tri$add>=2)
bottom=colSums(bchu_rthu_tri["bchu"] > 0) + colSums(bchu_rthu_tri["rthu"] > 0)-sum(bchu_rthu_tri$add>=2)

## jaccards
empirical=top/bottom
# empirical
# bchu 
# 0.1076923 

jaccards <- NA
for (i in 1:1000){
  
df.i <- transform( bchu_rthu_tri[,c(1,2)], bchu = sample(bchu) )
df.i$add <- apply(df.i, 1, function(i) sum(i > 0))
top=sum(df.i$add>=2)
bottom=colSums(df.i["bchu"] > 0) + colSums(df.i["rthu"] > 0) - sum(df.i$add>=2)
## jaccards
jaccards.i <- top/bottom
jaccards <- rbind(jaccards,jaccards.i)

}
p=length(which(jaccards>=empirical))/1000
p
## [1] 0.019
hist(jaccards)
abline(v=empirical)
