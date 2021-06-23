#twin model for position variability:estimates genetic, commmon and environmental contributions for each vertex
#input: tables
#output: contributions for each vertex


library(OpenMx)
library(umx)
mxOption(NULL, 'Number of Threads',2)
umx_set_cores()
#mxOption(NULL,"Default optimizer","SLSQP")# uncomment to use differnt optimizer
library(readxl)
library(parallel)
library(data.table)
library(GeneNet)
library(crunch)
readRenviron("../settings.Renviron")
projectfolder<-Sys.getenv("projectfolder")



#read in data 
target_dir<-paste(projectfolder,"7NETS_vertex/12_PositionVar_estimatedModels/",sep="",collapse=NULL)
dir.create(target_dir)
datapath  <-paste(projectfolder,'7NETS_vertex/11_PositionVar_cosine_GMR',sep="",collapse=NULL)
datafiles<-list.files(datapath,pattern='csv.gz',full.names=TRUE)
Files<-list.files(datapath,pattern='csv.gz',full.names=FALSE)
files<-tools::file_path_sans_ext(Files)
files<-tools::file_path_sans_ext(files)


################################################################################# used when not run for all vertices
#iterate=read.csv('tmp_nonzero_position.csv')# contains vertices where optimization got stuck and code is run again with different start values
#iterate=iterate['nonzero']# uncomment when rerunning for some vertices only, tmp_nonzero file created with find_nonzero_status.py
#tmpfiles=data.frame(Files)
#colnames(tmpfiles)=c('nonzero')
#tmp=rbind(iterate,tmpfiles)
#tmp=duplicated(tmp)
#nm=dim(iterate)
#start=nm[1]
#end=length(tmp)
#tmp=tmp[(start+1):end]
#tmptmp=which(tmp)

#################################################################################
for (i in 1:length(files)){ # comment out when rerunning for some vertices only
#for (i in tmptmp){        # uncomment when rerunning for some vertices only
gzf<-gzfile(datafiles[i])
data<-read.csv(gzf,header=T)


#split in MZ and DZ twins
MZ=data[data$Zygosity=='MZ',,drop=FALSE] #
DZ=data[data$Zygosity=='NotMZ',,drop=FALSE] #

####################################################
#split in pairs and singleton twins#
n_occur=data.frame(table(MZ$Mother_ID))
MZ_pair=MZ[MZ$Mother_ID %in% n_occur[n_occur$Freq==2,1],,drop=FALSE]
MZ_single=MZ[MZ$Mother_ID %in% n_occur[n_occur$Freq==1,1],,drop=FALSE]

n_occur=data.frame(table(DZ$Mother_ID))
DZ_pair=DZ[DZ$Mother_ID %in% n_occur[n_occur$Freq==2,1],,drop=FALSE]
DZ_single=DZ[DZ$Mother_ID %in% n_occur[n_occur$Freq==1,1],,drop=FALSE]
###################################################



#prepare data for model
MZ_even_ind=seq(2,86,2)#nr of twins
MZ_odd_ind=seq(1,86-1,2)
DZ_even_ind=seq(2,82,2)#nr of twins
DZ_odd_ind=seq(1,82-1,2)

MZ_twin1=MZ_pair[MZ_odd_ind,4:6,drop=FALSE]
MZ_twin2=MZ_pair[MZ_even_ind,4:6,drop=FALSE]
DZ_twin1=DZ_pair[DZ_odd_ind,4:6,drop=FALSE]
DZ_twin2=DZ_pair[DZ_even_ind,4:6,drop=FALSE]

T1=paste('Twin1_',seq(1,3,1),sep="")
T2=paste('Twin2_',seq(1,3,1),sep="")
colnames(MZ_twin1)=T1
colnames(MZ_twin2)=T2
colnames(DZ_twin1)=T1
colnames(DZ_twin2)=T2

##############################################

# add singleton twins to set as pairs, some are added as twin 1 others as twin 2
MZ_single1=MZ_single[1:13,4:6,drop=FALSE]
MZ_single2=MZ_single[14:26,4:6,drop=FALSE]
DZ_single1=DZ_single[1:18,4:6,drop=FALSE]
DZ_single2=DZ_single[19:37,4:6,drop=FALSE]
colnames(MZ_single1)=T1
colnames(MZ_single2)=T2
colnames(DZ_single1)=T1
colnames(DZ_single2)=T2

tmp1=data.frame(matrix(NA,ncol=3,nrow=13))
tmp2=data.frame(matrix(NA,ncol=3,nrow=13))
tmp3=data.frame(matrix(NA,ncol=3,nrow=18))
tmp4=data.frame(matrix(NA,ncol=3,nrow=19))
colnames(tmp1)=T2
colnames(tmp2)=T1
colnames(tmp3)=T2
colnames(tmp4)=T1
MZ_twin1=rbind(MZ_twin1,rbind(MZ_single1,tmp2))
MZ_twin2=rbind(MZ_twin2,rbind(tmp1,MZ_single2))
DZ_twin1=rbind(DZ_twin1,rbind(DZ_single1,tmp4))
DZ_twin2=rbind(DZ_twin2,rbind(tmp3,DZ_single2))


##############################################
MZdata=cbind(MZ_twin1,MZ_twin2)
DZdata=cbind(DZ_twin1,DZ_twin2)







#ACE Model#
start     <- diag(1,3,3)
#start <-diag(sqrt(diag(cov(DZ_twin2,use="complete"))))
#start<-start/sqrt(3) 
#startmean=colMeans(DZ_twin2,na.rm=TRUE) 
 

# Matrices declared to store a, c, and e Path Coefficients
alabel<-matrix( c('a11',NA,NA,
                  'a21','a22',NA,
                  'a31','a32','a33'),nrow=3,ncol=3, byrow=TRUE)
clabel<-matrix(c('c11',NA,NA,
                  'c21','c22',NA,
                  'c31','c32','c33'),nrow=3,ncol=3, byrow=TRUE)
elabel<-matrix(c('e11',NA,NA,
                  'e21','e22',NA,
                  'e31','e32','e33'),nrow=3,ncol=3, byrow=TRUE)

pathA <- mxMatrix(type='Lower',nrow=3,ncol=3, free=TRUE,values=start,label=alabel, name='a')
pathC<- mxMatrix(type='Lower',nrow=3,ncol=3, free=TRUE,values=start,label=clabel, name='c')
pathE <- mxMatrix(type='Lower',nrow=3,ncol=3, free=TRUE,values=start, label=elabel,name='e')

# Matrices generated to hold A, C, and E computed Variance Components
covA <- mxAlgebra( expression=a %*% t(a), name="A" )
covC <- mxAlgebra( expression=c %*% t(c), name="C" )
covE <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances
covP      <- mxAlgebra( expression=A+C+E, name="V" )
Asum      <- mxAlgebra(expression=sum(diag2vec(A))/sum(diag2vec(V)), name="Atot")
Csum      <- mxAlgebra(expression=sum(diag2vec(C))/sum(diag2vec(V)), name="Ctot")
Esum      <- mxAlgebra(expression=sum(diag2vec(E))/sum(diag2vec(V)), name="Etot")
ACdiff    <-mxAlgebra(expression=Atot-Ctot, name='diffAC')
matI      <- mxMatrix( type="Iden", nrow=3, ncol=3, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")


# Algebra for expected Mean
laMeMZ <- paste("Mean",seq(1,3,1),sep="_")
laMeMZ <-append(laMeMZ,laMeMZ)
laMeDZ<-laMeMZ
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=2*3, free=TRUE,labels=laMeMZ, name="expMeanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=2*3, free=TRUE,labels=laMeDZ,name="expMeanDZ" )


# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZ <- mxAlgebra( expression=rbind( cbind(V, A+C),
                                          cbind(A+C, V)), name="expCovMZ" )
covDZ <- mxAlgebra( expression=rbind( cbind(V, 0.5%x%A+ C),
                                          cbind(0.5%x%A+ C, V)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=MZdata,type="raw")
dataDZ    <- mxData( observed=DZdata,type="raw")

# Objective objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ",
                                  dimnames=append(T1[1:3],T2[1:3]))
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ",
                                  dimnames=append(T1[1:3],T2[1:3]))
funML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP,matI,invSD,Asum,Csum,Esum, ACdiff )
modelMZ   <- mxModel( pars, meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
fitML     <- mxFitFunctionMultigroup(c("MZ.fitfunction","DZ.fitfunction") )
AceModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )








#request confidence intervalls
AceModel <- mxModel(AceModel,mxCI(c('ACE.Atot','ACE.Ctot','ACE.Etot','ACE.diffAC')))


# Run ACE model
AceFit    <- mxRun(AceModel, intervals=T)
#confInts<-AceFit@output$confidenceIntervals
#print(AceFit@output$confidenceIntervals)
#summary(AceFit,SaturatedLikelihood=twinSatFit,SaturatedDoF=2249)







estVA     <- diag(mxEval(a%*%t(a), AceFit))             # additive genetic variance, a^2
estVC     <- diag(mxEval(c%*%t(c), AceFit))             # common environmental variance, c^2
estVE     <- diag(mxEval(e%*%t(e), AceFit))             # unique environmental variance, e^2       
estimates<-append(estVA,append(estVC,estVE))
estimates<-data.frame(estimates)
colnames(estimates)<-'estimate'



statuscode<-AceFit@output$status$code
statuscode<-rep(statuscode,3*3)    
Status<-data.frame(statuscode)
colnames(Status)<-'Status Code'


#tmptable<-confInts
tmptable<-cbind(estimates,Status)

saveRDS(AceFit,file=paste(target_dir,files[i],'.rds',sep=""))# save model
write.csv.gz(tmptable,paste(target_dir,files[i],'.csv.gz',sep=""))# save contributions in csv file
}



