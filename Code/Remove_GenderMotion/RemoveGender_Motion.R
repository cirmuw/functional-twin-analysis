#remove gender and motion for funtional connectivty before and after functionl aligment
# input raw table for each vertex containg correlation values
#output: tables with gender and motion removed
library(readxl)
library(xlsx)
library(crunch)

readRenviron("../settings.Renviron")
projectfolder<-Sys.getenv("projectfolder")

datapath<-paste(projectfolder,'7NETS_vertex/2_7nets_ROIs',"",collapse=Null)
#datapath<-paste(projectfolder,'7NETS_vertex/7_7nets_ROIs_cosine',"",collapse=Null)# change directory according to entagled or disentagled setting
datafiles<-list.files(datapath,pattern='csv.gz',full.names=TRUE)
files<-list.files(datapath,pattern='csv.gz',full.names=FALSE)
sexes<-read_excel('../../Deliveries/sexes.xlsx',1)
gzfile<-gzfile('../../Deliveries/MeanMotionRun1RL.csv.gz')
motion<-read.csv(gzfile)
gzfile<-gzfile('../../Deliveries/MeanMotionRun2RL.csv.gz')
motion2<-read.csv(gzfile)
tar_path<-paste(projectfolder,'7NETS_vertex/3_7nets_Fisher_GenderMotionEffectsRemoved/',sep="",collapse=NULL)
#tar_path<-paste(projectfolder,'7NETS_vertex/8_7nets_Fisher_GenderMotionEffectsRemoved_cosine/') # change directory according to entagled or disentagled setting
dir.create(tar_path)


for (i in 1:length(files)){
 #print(files[i])
 gzf<-gzfile(datafiles[i])
 data<-read.csv(gzf,header=T)
 ztrans<-atanh(data[,5:17,drop=FALSE])
 data<-cbind(motion,motion2,sexes,data[,1:4],ztrans)
 data_pair<-data[,-c(1,2,4,5,7,9)]                       
 linearMod <-lm(as.matrix(data_pair[,7:19])~Gender+Rel_Mean_Motion+Rel_Mean_Motion2,data_pair)
 residuals=residuals(linearMod)
 newdata_pair=cbind(data_pair[,4:6],residuals)
 outfile=paste(tar_path,files[i],sep="")
 write.csv.gz(newdata_pair,outfile)
}
