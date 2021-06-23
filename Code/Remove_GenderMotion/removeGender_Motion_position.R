#remove geder and motion for postion
# input raw table for each vertex containg correlation values
#output: tables with gender and motion removed
library(readxl)
library(xlsx)
library(crunch)

readRenviron("../settings.Renviron")
projectfolder<-Sys.getenv("projectfolder")

datapath<-paste(projectfolder,'7NETS_vertex/10_PositionVar_cosine',sep="",collapse=NULL)
datafiles<-list.files(datapath,pattern='csv.gz',full.names=TRUE)
files<-list.files(datapath,pattern='csv.gz',full.names=FALSE)
sexes<-read_excel('../../Deliveries/sexes.xlsx',1)
gzfile<-gzfile('../../Deliveries/MeanMotionRun1LR.csv.gz')
motion<-read.csv(gzfile)
gzfile<-gzfile('../../Deliveries/MeanMotionRun2LR.csv.gz')
motion2<-read.csv(gzfile)
tar_path<-paste(projectfolder,'7NETS_vertex/11_PositionVar_cosine_GMR/',sep="",collapse=NULL)
dir.create(tar_path)



for (i in 1:length(files)){
 #print(files[i])
 gzf<-gzfile(datafiles[i])
 data<-read.csv(gzf,header=T)
 coords<-data[,5:7]
 data<-cbind(motion,motion2,sexes,data[,1:4],coords)
 data_pair<-data[,-c(1,2,4,5,7,9)]                       
 linearMod <-lm(as.matrix(data_pair[,7:9])~Gender+Rel_Mean_Motion+Rel_Mean_Motion2,data_pair)
 residuals=residuals(linearMod)
 newdata_pair=cbind(data_pair[,4:6],residuals)
 outfile=paste(tar_path,files[i],sep="")
 write.csv.gz(newdata_pair,outfile)
}
