## load library and set the working directory
## install.packages("pavo")
library(pavo)
library(dplyr)
setwd("C:/Users/kirad/Dropbox/Texas/Projects/hummingbirds/dive_spec_data/")

## bring in all files in the working directory 
specs<-getspec(where = getwd(), ext = "txt", lim = c(300,700), subdir = TRUE, subdir.names = TRUE)

## limit to dives
dives <- specs[grep("Dive", names(specs))]
dives <- cbind(specs[,1],dives)
dives <- as.rspec(dives)
dives <- procspec(dives,fixneg = "zero")

## quick check to see what the data look like
explorespec(dives[,1:length(dives)],by=9)
explorespec(shuttle[,1:length(shuttle)],by=4)

## average by bird and angle dives
#inds <- c("CE12K01","CE19K01","CE21K01","CE21K02","CE22K01","CE23K01","CE24K01","CE29K01","CE29K02","CE30K02","CE31K01","CF06K02","CF07K01","CF07K02","CF07K03","CE23K01_redo","CE22K01_redo","CE29K02_redo","CE31K01_redo")
angles <- c("-80_Dive","-60_Dive","-40_Dive","-20_Dive","_0_Dive","_20_Dive","_40_Dive","_60_Dive","_80_Dive")
inds <- c("CE24K01","CE29K01","CE29K02_redo","CE23K01_redo","EE15K01",
          "CF06K02","CF07K01","CF07K02","CF07K03")

mspecs1 <- NA

for (i in 1:length(inds)){
  
  ind.i <- inds[i]
  ind.i.df <- dives[grep(ind.i, names(dives))]
  
  for (j in 1:length(angles)){
    
    angle.j <- angles[j]
    angle.j.df <- ind.i.df[grep(angle.j, names(ind.i.df))]
    angle.j.df <- cbind(specs[,1],angle.j.df)
    angle.j.df <- as.rspec(angle.j.df)
    
    mspecs1.i <-aggspec(angle.j.df,FUN=mean)
    name.i <- paste(ind.i,angle.j,sep="")
    names(mspecs1.i) <- c("wl",paste(ind.i,angle.j))
    mspecs1 <- cbind(mspecs1,mspecs1.i)
    
  }
}

mspecs1 <- mspecs1 %>% select(-contains("wl"))
mspecs1 <- mspecs1[,-1]
mspecs1 <- cbind(specs[,1],mspecs1)
mspecs1 <- as.rspec(mspecs1)

angles <- c("Dive_-80","Dive_-60","Dive_-40","Dive_-20","Dive_0","Dive_20","Dive_40","Dive_60","Dive_80")
inds <- c("DE15K01","BE10K07","DE15K04","DE15K05","DE15K02","BE10K06",
          "EE19K01","DF09K01","DF10K01","DF09K03","DF09K02","AE18K03","AE18K04","AE18K11","AE18K01")

mspecs2 <- NA

for (i in 1:length(inds)){
  
  ind.i <- inds[i]
  ind.i.df <- dives[grep(ind.i, names(dives))]
  
  for (j in 1:length(angles)){
    
    angle.j <- angles[j]
    angle.j.df <- ind.i.df[grep(angle.j, names(ind.i.df))]
    angle.j.df <- cbind(specs[,1],angle.j.df)
    angle.j.df <- as.rspec(angle.j.df)
    
    mspecs2.i <-aggspec(angle.j.df,FUN=mean)
    name.i <- paste(ind.i,angle.j,sep="")
    names(mspecs2.i) <- c("wl",paste(ind.i,angle.j))
    mspecs2 <- cbind(mspecs2,mspecs2.i)
    
  }
}

mspecs2 <- mspecs2 %>% select(-contains("wl"))
mspecs2 <- mspecs2[,-1]
mspecs2 <- cbind(specs[,1],mspecs2)
mspecs2 <- as.rspec(mspecs2)

mspecs <- merge(mspecs1,mspecs2,by="wl")
mspecs <- as.rspec(mspecs)
rm(mspecs1,mspecs1.i,mspecs2,mspecs2.i,ind.i,ind.i.df,angle.j,angle.j.df,i,j,name.i,inds,angles)

## apply smoothing term
spec.sm<-procspec(mspecs,opt='smooth',span=0.2)
explorespec(spec.sm,by=9)

#calculate photon catches at each photoreceptor, here using the average avian UV visual system 
# first need to quantify receptor excitation and then consider how the signal is being processed
vm.spec.sm <- vismodel(spec.sm, visual="avg.uv",relative=T,achromatic = 'st.dc')

#calculate position in tetra-hedral colour space (see Stoddard & Prum 2008 for details)
tcs.all<-colspace(vm.spec.sm)
plot(tcs.all,
     pch = 21,
     bg = spec2rgb(spec.sm),
     perspective = TRUE,
     range = c(1, 2),
     cex = 1.5
)

##
names <- rownames(tcs.all)
rownames(tcs.all) <- NULL
data <- cbind(names,tcs.all)
data$names <- gsub("_","",data$names)
data$names <- gsub("Dive","",data$names)
foo <- data.frame(do.call('rbind', strsplit(as.character(data$names),' ',fixed=TRUE)))
names(foo) <- c("ind","angle")
data <-cbind(foo,data)

black <- c("CE24K01","CE29K01","CE29K02_redo","CE23K01_redo","EE15K01","DE15K01","BE10K07","DE15K04","DE15K05","DE15K02","BE10K06")
ruby <- c("CF06K02","CF07K01","CF07K02","CF07K03","EE19K01","DF09K01","DF10K01","DF09K03","DF09K02","AE18K03","AE18K04","AE18K11","AE18K01")

black_df <- data[data$ind %in% black,]
ruby_df <- data[data$ind %in% ruby,]
black_df$angle <- gsub(".x","",black_df$angle)
black_df$angle <- gsub(".y","",black_df$angle)
black_df$angle <- as.numeric(black_df$angle)
ruby_df$angle <- gsub(".x","",ruby_df$angle)
ruby_df$angle <- gsub(".y","",ruby_df$angle)
ruby_df$angle <- as.numeric(ruby_df$angle)

write.table(data,"./output/color_by_angle.csv",sep=",",row.names = F)

## plot
vars <- ruby_df[,c(2,4:7,20)]
melt <- melt(setDT(vars),id="angle")
melt$angle <- as.numeric(melt$angle)
ggplot(melt,aes(x=angle,y=value,group=variable))+geom_smooth(aes(col=variable,fill=variable))+geom_point(aes(col=variable),size=1)+
  scale_color_manual(values = c("#E76BF3","#00B0F6","#00BF7D","#F8766D","#A3A500"))+
  scale_fill_manual(values = c("#E76BF3","#00B0F6","#00BF7D","#F8766D","#A3A500"))+
  ylim(0,0.8)+
  ylab("mean relative cone stimulation")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")

vars <- black_df[,c(2,4:7,20)]
melt$angle <- as.numeric(melt$angle)
melt <- melt(setDT(vars),id="angle")
ggplot(melt,aes(x=angle,y=value,group=variable))+geom_smooth(aes(col=variable,fill=variable))+geom_point(aes(col=variable),size=1)+
  scale_color_manual(values = c("#E76BF3","#00B0F6","#00BF7D","#F8766D","#A3A500"))+
  scale_fill_manual(values = c("#E76BF3","#00B0F6","#00BF7D","#F8766D","#A3A500"))+
  ylim(0,0.8)+
  ylab("mean relative cone stimulation")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")

#find max value
ruby_df[as.logical(ave(ruby_df$l, ruby_df$ind, FUN = function(x) x == max(x))),]
black_df[as.logical(ave(black_df$s, black_df$ind, FUN = function(x) x == max(x))),]
