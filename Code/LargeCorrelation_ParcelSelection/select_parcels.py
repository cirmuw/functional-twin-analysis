# script to select representative parcels of Yeo networks based on a reference subject
#input:time series of reference subject, Schaefer 600 parcelation for left and right hemisphere
#output: representative parcel for each network of Yeo 7 nets; left and right hemi are seperated, which gives 14 representative parcels
         #id of vertices without signal (medial wall) is also output and needed for later
import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
from myfuns import calc_vertex_correlations_combined_runs
from myfuns import calc_vertex_correlations
import matplotlib.pyplot as plt
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import settings as s

reference_subjectpath=s.HCProot+'HCP_3T_RESTA_fmri/101915/' #folder containg fmri runs
lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'#parcelation left hemi
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'#parcelation right hemi 







# load time series and resize
lh1=nib.load(reference_subjectpath+'/lh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh').get_data()
lh2=nib.load(reference_subjectpath+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh').get_data()
rh1=nib.load(reference_subjectpath+'/rh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh').get_data()
rh2=nib.load(reference_subjectpath+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh').get_data()

lh1.resize(lh1.shape[1],lh1.shape[3])
lh2.resize(lh2.shape[1],lh2.shape[3])
rh1.resize(rh1.shape[1],rh1.shape[3])
rh2.resize(rh2.shape[1],rh2.shape[3])

lh=np.concatenate((lh1,lh2),axis=1)
rh=np.concatenate((rh1,rh2),axis=1)

lh_zerovertex=~lh.any(axis=1)
rh_zerovertex=~rh.any(axis=1)
deleteindlh=np.where(lh_zerovertex)[0] #remove vertices of medial wall whitout signal
deleteindrh=np.where(rh_zerovertex)[0]

lh=np.delete(lh,deleteindlh,axis=0)
rh=np.delete(rh,deleteindrh,axis=0)

np.save('../../Deliveries/0verticeslh',deleteindlh)#save id of removed vertices for later
np.save('../../Deliveries/0verticesrh',deleteindrh)



#load schaefer parcelation
lhannot=nib.freesurfer.io.read_annot(lhparpath)
lhlabels=lhannot[0]

rhannot=nib.freesurfer.io.read_annot(rhparpath)
rhlabels=rhannot[0]

lhlabels=np.delete(lhlabels,deleteindlh,0)
rhlabels=np.delete(rhlabels,deleteindrh,0)





#get mean time series of parcels
lhframe=pd.DataFrame(np.transpose(lh))
lhframe=lhframe.groupby(lhlabels,axis=1).mean()
rhframe=pd.DataFrame(np.transpose(rh))
rhframe=rhframe.groupby(rhlabels,axis=1).mean()

#distinguish parcels based on color table
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




#find best representative parcel of network: the one which correlations on average the highest with other network parcels is chosen
lchosenrois=[]
for i in range(0,8):
   tmplh=lhframe.ix[:,lhnetwork[i]:lhnetwork[i+1]-1]
   meancorr=tmplh.corr().mean()
   chosenroi=meancorr.idxmax()
   lchosenrois.append(chosenroi)


rchosenrois=[]
for i in range(0,8):   
   tmprh=rhframe.ix[:,rhnetwork[i]:rhnetwork[i+1]-1]
   meancorr=tmprh.corr().mean()
   chosenroi=meancorr.idxmax()
   rchosenrois.append(chosenroi)


np.save('../../Deliveries/chosenroislh',lchosenrois)#save id of chosen rois
np.save('../../Deliveries/chosenroisrh',rchosenrois)
