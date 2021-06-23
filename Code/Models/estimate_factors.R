#twin model for functional connectivity before and after functional alignment:estimates genetic, commmon and environmental contributions for each vertex
#input: tables
#output: contributions for each vertex


library(OpenMx)
library(umx)
mxOption(NULL, 'Number of Threads',2)
umx_set_cores()
#mxOption(NULL,"Default optimizer","SLSQP") #uncomment to use other optimizer
library(readxl)
library(parallel)
library(data.table)
library(GeneNet)
library(crunch)
readRenviron("../settings.Renviron")
projectfolder<-Sys.getenv("projectfolder")



#read in data 
target_dir<-paste(projectfolder,'7NETS_vertex/4_estimatedModels/',sep="",collapse=NULL)# entangled connectivity
#target_dir<-paste(projectfolder,'7NETS_vertex/9_estimatedModels_functional_cosine/',sep="",collapse=NULL)# disentangled connectivity
dir.create(target_dir)
datapath  <-paste(projectfolder,'7NETS_vertex/3_7nets_Fisher_GenderMotionEffectsRemoved',sep="",collapse=NULL)# entangled connectivity
#datapath  <-paste(projectfolder,'7NETS_vertex/8_7nets_Fisher_GenderMotionEffectsRemoved_cosine',sep="",collapse=NULL)# disentangled connectivity
datafiles<-list.files(datapath,pattern='csv.gz',full.names=TRUE)
Files<-list.files(datapath,pattern='csv.gz',full.names=FALSE)
files<-tools::file_path_sans_ext(Files)
files<-tools::file_path_sans_ext(files)



################################################################################# needed when code is run only for some vertices
#iterate=read.csv('tmp_nonzero_anat.csv') # contains vertices where optimization got stuck and code is run again with different start values 
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
#for (i in tmptmp){  # uncomment when rerunning for some vertices only
gzf<-gzfile(datafiles[i])
data<-read.csv(gzf,header=T)


#split in MZ and DZ twins
MZ=data[data$Zygosity=='MZ',,drop=FALSE] #
DZ=data[data$Zygosity=='NotMZ',,drop=FALSE] #


#split in pairs and singleton twins#
n_occur=data.frame(table(MZ$Mother_ID))
MZ_pair=MZ[MZ$Mother_ID %in% n_occur[n_occur$Freq==2,1],,drop=FALSE]
MZ_single=MZ[MZ$Mother_ID %in% n_occur[n_occur$Freq==1,1],,drop=FALSE]

n_occur=data.frame(table(DZ$Mother_ID))
DZ_pair=DZ[DZ$Mother_ID %in% n_occur[n_occur$Freq==2,1],,drop=FALSE]
DZ_single=DZ[DZ$Mother_ID %in% n_occur[n_occur$Freq==1,1],,drop=FALSE]




#prepare data for model
MZ_even_ind=seq(2,86,2)# nr of twins
MZ_odd_ind=seq(1,86-1,2)
DZ_even_ind=seq(2,82,2)# nr of twins
DZ_odd_ind=seq(1,82-1,2)

MZ_twin1=MZ_pair[MZ_odd_ind,4:16,drop=FALSE]
MZ_twin2=MZ_pair[MZ_even_ind,4:16,drop=FALSE]
DZ_twin1=DZ_pair[DZ_odd_ind,4:16,drop=FALSE]
DZ_twin2=DZ_pair[DZ_even_ind,4:16,drop=FALSE]

T1=paste('Twin1_',seq(1,13,1),sep="")
T2=paste('Twin2_',seq(1,13,1),sep="")
colnames(MZ_twin1)=T1
colnames(MZ_twin2)=T2
colnames(DZ_twin1)=T1
colnames(DZ_twin2)=T2


# add singleton twins to set as pairs, some are added as twin 1 others as twin 2
MZ_single1=MZ_single[1:13,4:16,drop=FALSE]
MZ_single2=MZ_single[14:26,4:16,drop=FALSE]
DZ_single1=DZ_single[1:18,4:16,drop=FALSE]
DZ_single2=DZ_single[19:37,4:16,drop=FALSE]
colnames(MZ_single1)=T1
colnames(MZ_single2)=T2
colnames(DZ_single1)=T1
colnames(DZ_single2)=T2

tmp1=data.frame(matrix(NA,ncol=13,nrow=13))
tmp2=data.frame(matrix(NA,ncol=13,nrow=13))
tmp3=data.frame(matrix(NA,ncol=13,nrow=18))
tmp4=data.frame(matrix(NA,ncol=13,nrow=19))
colnames(tmp1)=T2
colnames(tmp2)=T1
colnames(tmp3)=T2
colnames(tmp4)=T1
MZ_twin1=rbind(MZ_twin1,rbind(MZ_single1,tmp2))
MZ_twin2=rbind(MZ_twin2,rbind(tmp1,MZ_single2))
DZ_twin1=rbind(DZ_twin1,rbind(DZ_single1,tmp4))
DZ_twin2=rbind(DZ_twin2,rbind(tmp3,DZ_single2))



MZdata=cbind(MZ_twin1,MZ_twin2)
DZdata=cbind(DZ_twin1,DZ_twin2)


#ACE Model#
start     <- diag(1,13,13) 
 

 

# Matrices declared to store a, c, and e Path Coefficients
alabel<-matrix( c('a11',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'a21','a22',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'a31','a32','a33',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'a41','a42','a43','a44',NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'a51','a52','a53','a54','a55',NA,NA,NA,NA,NA,NA,NA,NA,
                  'a61','a62','a63','a64','a65','a66',NA,NA,NA,NA,NA,NA,NA,
                  'a71','a72','a73','a74','a75','a76','a77',NA,NA,NA,NA,NA,NA,
                  'a81','a82','a83','a84','a85','a86','a87','a88',NA,NA,NA,NA,NA, 
                  'a91','a92','a93','a94','a95','a96','a97','a98','a99',NA,NA,NA,NA, 
                  'a101','a102','a103','a104','a105','a106','a107','a108','a109','a1010',NA,NA,NA, 
                  'a111','a112','a113','a114','a115','a116','a117','a118','a119','a1110','a1111',NA,NA, 
                  'a121','a122','a123','a124','a125','a126','a127','a128','a129','a1210','a1211','a1212',NA, 
                  'a131','a132','a133','a134','a135','a136','a137','a138','a139','a1310','a1311','a1312','a1313'),nrow=13,ncol=13, byrow=TRUE)
clabel<-matrix(c('c11',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'c21','c22',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'c31','c32','c33',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'c41','c42','c43','c44',NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'c51','c52','c53','c54','c55',NA,NA,NA,NA,NA,NA,NA,NA,
                  'c61','c62','c63','c64','c65','c66',NA,NA,NA,NA,NA,NA,NA,
                  'c71','c72','c73','c74','c75','c76','c77',NA,NA,NA,NA,NA,NA,
                  'c81','c82','c83','c84','c85','c86','c87','c88',NA,NA,NA,NA,NA, 
                  'c91','c92','c93','c94','c95','c96','c97','c98','c99',NA,NA,NA,NA, 
                  'c101','c102','c103','c104','c105','c106','c107','c108','c109','c1010',NA,NA,NA, 
                  'c111','c112','c113','c114','c115','c116','c117','c118','c119','c1110','c1111',NA,NA, 
                  'c121','c122','c123','c124','c125','c126','c127','c128','c129','c1210','c1211','c1212',NA, 
                  'c131','c132','c133','c134','c135','c136','c137','c138','c139','c1310','c1311','c1312','c1313'),nrow=13,ncol=13, byrow=TRUE)
elabel<-matrix(c('e11',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'e21','e22',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'e31','e32','e33',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'e41','e42','e43','e44',NA,NA,NA,NA,NA,NA,NA,NA,NA,
                  'e51','e52','e53','e54','e55',NA,NA,NA,NA,NA,NA,NA,NA,
                  'e61','e62','e63','e64','e65','e66',NA,NA,NA,NA,NA,NA,NA,
                  'e71','e72','e73','e74','e75','e76','e77',NA,NA,NA,NA,NA,NA,
                  'e81','e82','e83','e84','e85','e86','e87','e88',NA,NA,NA,NA,NA, 
                  'e91','e92','e93','e94','e95','e96','e97','e98','e99',NA,NA,NA,NA, 
                  'e101','e102','e103','e104','e105','e106','e107','e108','e109','e1010',NA,NA,NA, 
                  'e111','e112','e113','e114','e115','e116','e117','e118','e119','e1110','e1111',NA,NA, 
                  'e121','e122','e123','e124','e125','e126','e127','e128','e129','e1210','e1211','e1212',NA, 
                  'e131','e132','e133','e134','e135','e136','e137','e138','e139','e1310','e1311','e1312','e1313'),nrow=13,ncol=13, byrow=TRUE)


pathA <- mxMatrix(type='Lower',nrow=13,ncol=13, free=TRUE,values=start,label=alabel, name='a')
pathC<- mxMatrix(type='Lower',nrow=13,ncol=13, free=TRUE,values=start,label=clabel, name='c')
pathE <- mxMatrix(type='Lower',nrow=13,ncol=13, free=TRUE,values=start, label=elabel,name='e')

# Matrices generated to hold A, C, and E computed Variance Components
covA <- mxAlgebra( expression=a %*% t(a), name="A" )
covC <- mxAlgebra( expression=c %*% t(c), name="C" )
covE <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances
covP <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=13, ncol=13, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")


# Algebra for expected Mean
laMeMZ <- paste("Mean",seq(1,13,1),sep="_")
laMeMZ <-append(laMeMZ,laMeMZ)
laMeDZ<-laMeMZ
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=2*13, free=TRUE, labels=laMeMZ, name="expMeanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=2*13, free=TRUE, labels=laMeDZ, name="expMeanDZ" )


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
                                  dimnames=append(T1[1:13],T2[1:13]))
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ",
                                  dimnames=append(T1[1:13],T2[1:13]))
funML     <- mxFitFunctionML()

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP,matI,invSD )
modelMZ   <- mxModel( pars, meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
fitML     <- mxFitFunctionMultigroup(c("MZ.fitfunction","DZ.fitfunction") )
AceModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, fitML )





# Run ACE model
AceFit    <- mxRun(AceModel,intervals=T)



estVA     <- diag(mxEval(a%*%t(a), AceFit))             # additive genetic variance, a^2
estVC     <- diag(mxEval(c%*%t(c), AceFit))             # common environmental variance, c^2
estVE     <- diag(mxEval(e%*%t(e), AceFit))             # unique environmental variance, e^2       
estimates<-append(estVA,append(estVC,estVE))
estimates<-data.frame(estimates)
colnames(estimates)<-'estimate'



statuscode<-AceFit@output$status$code
statuscode<-rep(statuscode,3*13)       
Status<-data.frame(statuscode)
colnames(Status)<-'Status Code'


tmptable<-cbind(estimates,Status)

saveRDS(AceFit,file=paste(target_dir,files[i],'.rds',sep=""))# save model
write.csv.gz(tmptable,paste(target_dir,files[i],'.csv.gz',sep=""))#save contibutions in CSV file
}


