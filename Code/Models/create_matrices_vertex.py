#creates arrayes containing contribuions for each vertex
#input vertex vise tables with contribution values
#output: surface maps of contributions
import nibabel as nib
import numpy as np
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s
import scipy.io
from myfuns import myinsert

mode='functional' #set the mode, 'functional' or 'position'


lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'# parcelation
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'
zerovertexlh=np.load('../../Deliveries/0verticeslh.npy')# vertices without signal
zerovertexrh=np.load('../../Deliveries/0verticesrh.npy')

lhannot=nib.freesurfer.io.read_annot(lhparpath)
lhlabels=lhannot[0]
rhannot=nib.freesurfer.io.read_annot(rhparpath)
rhlabels=rhannot[0]
labelslh=np.delete(lhlabels,zerovertexlh,0)
labelsrh=np.delete(rhlabels,zerovertexrh,0)



if mode =='functional':
  path=s.projectfolder+'7NETS_vertex/4_estimatedModels/'# path to entagled or disentagled connectivity strength
  #path=s.projectfolder+'7NETS_vertex/9_estimatedModels_functional_cosine/'
  targetdir='../../Deliveries/ACE_matrices_anat.mat'
  #targetdir='../../Deliveries/ACE_matrices_func.mat'

  lA=np.zeros((len(labelslh)))
  lC=np.zeros((len(labelslh)))
  lE=np.zeros((len(labelslh)))
  lA[:]=-2
  lC[:]=-2
  lE[:]=-2
  for i in range(len(labelslh)):
   if labelslh[i]!=0:
    data=pd.read_csv(path+'lh_vertex'+str(i+1)+'.csv.gz',compression='gzip')
    data=data['estimate']
    dataA=data.iloc[0:13]
    dataC=data.iloc[13:26]
    dataE=data.iloc[26:39]
    a=dataA.sum(axis=0)
    c=dataC.sum(axis=0)
    e=dataE.sum(axis=0)
    lA[i]=a/(a+c+e)
    lC[i]=c/(a+c+e)
    lE[i]=e/(a+c+e)
    

  rA=np.zeros((len(labelsrh)))
  rC=np.zeros((len(labelsrh)))
  rE=np.zeros((len(labelsrh)))
  rA[:]=-2
  rC[:]=-2
  rE[:]=-2
  for i in range(len(labelsrh)):
   if labelsrh[i]!=0:
    data=pd.read_csv(path+ 'rh_vertex'+str(i+1)+'.csv.gz',compression='gzip')
    data=data['estimate']
    dataA=data.iloc[0:13]
    dataC=data.iloc[13:26]
    dataE=data.iloc[26:39]
    a=dataA.sum(axis=0)
    c=dataC.sum(axis=0)
    e=dataE.sum(axis=0)
    rA[i]=a/(a+c+e)
    rC[i]=c/(a+c+e)
    rE[i]=e/(a+c+e)
    
  lA=myinsert(np.expand_dims(lA,axis=1),zerovertexlh,-2,0)
  lC=myinsert(np.expand_dims(lC,axis=1),zerovertexlh,-2,0)
  lE=myinsert(np.expand_dims(lE,axis=1),zerovertexlh,-2,0)
  rA=myinsert(np.expand_dims(rA,axis=1),zerovertexrh,-2,0)
  rC=myinsert(np.expand_dims(rC,axis=1),zerovertexrh,-2,0)
  rE=myinsert(np.expand_dims(rE,axis=1),zerovertexrh,-2,0)
  


  ACE_dict={"lA":lA,"rA":rA,"lC":lC,"rC":rC,"lE":lE,"rE":rE}
  scipy.io.savemat(targetdir,ACE_dict)







if mode =='position':
  path=s.projectsfolder+'7NETS_vertex/12_PositionVar_estimatedModels/'
  targetdir='../../Deliveries/ACE_matrices_pos.mat'

  lA=np.zeros((len(labelslh)))
  lC=np.zeros((len(labelslh)))
  lE=np.zeros((len(labelslh)))
  lA[:]=-2
  lC[:]=-2
  lE[:]=-2
  for i in range(len(labelslh)):
   if labelslh[i]!=0:
    data=pd.read_csv(path+'lh_'+str(i+1)+'_mean_position.csv.gz',compression='gzip')
    data=data['estimate']
    a=np.sum(data.iloc[0:3])
    c=np.sum(data.iloc[3:6])
    e=np.sum(data.iloc[6:9])
    lA[i]=a/(a+c+e)
    lC[i]=c/(a+c+e)
    lE[i]=e/(a+c+e)


  rA=np.zeros((len(labelsrh)))
  rC=np.zeros((len(labelsrh)))
  rE=np.zeros((len(labelsrh)))
  rA[:]=-2
  rC[:]=-2
  rE[:]=-2
  for i in range(len(labelsrh)):
   if labelsrh[i]!=0:
    data=pd.read_csv(path+ 'rh_'+str(i+1)+'_mean_position.csv.gz',compression='gzip')
    data=data['estimate']
    a=np.sum(data.iloc[0:3])
    c=np.sum(data.iloc[3:6])
    e=np.sum(data.iloc[6:9])
    rA[i]=a/(a+c+e)
    rC[i]=c/(a+c+e)
    rE[i]=e/(a+c+e)
 
  lA=myinsert(np.expand_dims(lA,axis=1),zerovertexlh,-2,0)
  lC=myinsert(np.expand_dims(lC,axis=1),zerovertexlh,-2,0)
  lE=myinsert(np.expand_dims(lE,axis=1),zerovertexlh,-2,0)
  rA=myinsert(np.expand_dims(rA,axis=1),zerovertexrh,-2,0)
  rC=myinsert(np.expand_dims(rC,axis=1),zerovertexrh,-2,0)
  rE=myinsert(np.expand_dims(rE,axis=1),zerovertexrh,-2,0)
  

  ACE_dict={"lA":lA,"rA":rA,"lC":lC,"rC":rC,"lE":lE,"rE":rE}
  scipy.io.savemat(targetdir,ACE_dict)
