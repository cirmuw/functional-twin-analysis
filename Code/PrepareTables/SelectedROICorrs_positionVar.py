#script to create tabels containig x, y and z coordinates of functionally corresponding vertices (position variability) for each twin, one table per vertex
#input:id of functionally corresponding vetices of each twin to reference
#output: tables with vertex position in each subject, one table per vetex
import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s
import pickle




#paths to subject data,id of vertices without signal, surface file, parcelation, chosen rois
infile =s.HCP_information_sheet_path    #\
subjectpath1=s.HCProot+'HCP_3T_RESTA_fmri/'# used obtain subject ids
subjectpath2=s.HCProot+'HCP_3T_RESTB_fmri/'#/
source_dir=s.projectfolder+'7NETS_vertex/5_7nets_corresponding/' #  path containing id of functionally corresponding vetices of each twin to reference
target_dir=s.projectfolder+'/7NETS_vertex/10_PositionVar_cosine/'# output tables with vertex position in each subject
if not os.path.exists(target_dir):
   os.mkdir(target_dir)

zerovertexlh=np.load('../../Deliveries/0verticeslh.npy')#ids of vertices without signal
zerovertexrh=np.load('../../Deliveries/0verticesrh.npy')

surfacedirlh='../../Deliveries/fsaverage4/lh.inflated' # surface on which vertex coordinates are based
surfacedirrh='../../Deliveries/fsaverage4/rh.inflated'
lhsurf=nib.freesurfer.io.read_geometry(surfacedirlh)
rhsurf=nib.freesurfer.io.read_geometry(surfacedirrh)
lhsurf=lhsurf[0]
lhsurf=np.delete(lhsurf,zerovertexlh,0)   
rhsurf=rhsurf[0]
rhsurf=np.delete(rhsurf,zerovertexrh,0) 
surf=np.concatenate([lhsurf,rhsurf],axis=0) 

lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'
lhannot=nib.freesurfer.io.read_annot(lhparpath)
lhlabels=lhannot[0]
rhannot=nib.freesurfer.io.read_annot(rhparpath)
rhlabels=rhannot[0]
labelslh=np.delete(lhlabels,zerovertexlh,0)
labelsrh=np.delete(rhlabels,zerovertexrh,0)

lhrois=list(np.load('../../Deliveries/chosenroislh.npy'))#save id of chosen rois
rhrois=list(np.load('../../Deliveries/chosenroisrh.npy'))
lhrois=lhrois[1:]
rhrois=rhrois[1:]
nameslhrois=['l_'+str(s) for s in lhrois]
namesrhrois=['r_'+str(s) for s in rhrois]


#get assigenment of parcels to yeo nets based on color table
lhnetwork=np.zeros((9))
rhnetwork=np.zeros((9))
lhnetwork[8]=301
rhnetwork[8]=301
c1=1
c2=1

for i in range(1,301):
    if abs(lhannot[1][i][0]-lhannot[1][i-1][0])>5:
         lhnetwork[c1]=int(i)
         c1=c1+1
      
    if abs(rhannot[1][i][0]-rhannot[1][i-1][0])>5:
          rhnetwork[c2]=int(i)
          c2=c2+1

   

#Get paths to mgh-files of available subjects
xl=pd.ExcelFile(infile)
dataframe1=xl.parse('Sheet1')
isNotTwin=dataframe1['Twin_Stat']=='NotTwin'
isNotTwin=np.where(isNotTwin)[0]
dataframe2=dataframe1.drop(isNotTwin,0)
Subjects=dataframe2['Subject'].values

path1=[]
path2=[]
for i in range(Subjects.shape[0]):
    path1.append(subjectpath1+str(Subjects[i]))
    path2.append(subjectpath2+str(Subjects[i]))


truesubjects=[]

for i in range(Subjects.shape[0]):
    if os.path.isdir(path1[i])==True:
        truesubjects.append(Subjects[i])
    if os.path.isdir(path2[i])==True:
        truesubjects.append(Subjects[i])



name=['Subject','Zygosity','Mother_ID']
nonvertexdat=np.zeros((len(truesubjects),3),dtype=object)

for j in range(len(labelslh)):
      if labelslh[j]!=0:
          positionvar=[]         
          for i in range(len(truesubjects)):
              functional=pickle.load(open(source_dir+'lh_'+str(j+1)+'correspondingvertices.p','rb'))
              index=np.where(functional[1]==-1)[0]
              index=functional[0][i][index]
              index=index[0]
              coords=surf[index]
              positionframe=pd.DataFrame(coords) 
              positionframe.columns=['x','y','z']
              positionvar.append(positionframe)
             


              if j==0:
                index=dataframe2[dataframe2['Subject']==truesubjects[i]].index.tolist()
                tmp1=np.array([str(truesubjects[i]),dataframe2['Zygosity'][index].values[0], str(dataframe2['Mother_ID'][index].values[0])])
                nonvertexdat[i,:]=tmp1
                         
          nonvertextable=pd.DataFrame(data=nonvertexdat)
          nonvertextable.columns=name
          positionframe=pd.concat(positionvar,axis=0,ignore_index=True)
          table=pd.concat([nonvertextable,positionframe],axis=1)
          table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
          table.reset_index(inplace=True)
          table=table.drop('index',axis=1) 
          writefile=target_dir+'lh_'+str(j+1)+'_mean_position.csv.gz'
          table.to_csv(writefile, compression='gzip')
for j in range(len(labelsrh)):
      if labelsrh[j]!=0:
          positionvar=[]         
          for i in range(len(truesubjects)):
              functional=pickle.load(open(source_dir+'rh_'+str(j+1)+'correspondingvertices.p','rb'))
              index=np.where(functional[1]==-1)[0]
              index=functional[0][i][index]
              index=index[0]
              coords=surf[index]
              positionframe=pd.DataFrame(coords) 
              positionframe.columns=['x','y','z']
              positionvar.append(positionframe)
             

                         
          nonvertextable=pd.DataFrame(data=nonvertexdat)
          nonvertextable.columns=name
          positionframe=pd.concat(positionvar,axis=0,ignore_index=True)
          table=pd.concat([nonvertextable,positionframe],axis=1)
          table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
          table.reset_index(inplace=True)
          table=table.drop('index',axis=1) 
          writefile=target_dir+'rh_'+str(j+1)+'_mean_position.csv.gz'
          table.to_csv(writefile, compression='gzip')
print('Finished')
