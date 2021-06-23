#aligenment of embeddings using MSM after initial (rotation only) alignment
# input: initially aligned embeddings, are used as features for alignment with MSM
#output: functionally aligend subjects, output by MSM


import numpy as np
import nibabel as nib
import pandas as pd
from sklearn.neighbors import NearestNeighbors 
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s
import pickle
from myfuns import myinsert
import vtk
from vtk.util import numpy_support as VN
from sklearn.neighbors import NearestNeighbors 





#paths
lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'#parcellation
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'
lhmesh='../../Deliveries/fsaverage4/lh.sphere.vtk'
rhmesh='../../Deliveries/fsaverage4/rh.sphere.vtk'
tar_dir=s.projectfolder+'7NETS_vertex/5_7nets_embedding/'# root for input and output
if not os.path.exists(tar_dir+'registration/'):
   os.mkdir(tar_dir+'registration/')
if not os.path.exists(tar_dir+'registration/cosine_subject_txt_3feat/'):
   os.mkdir(tar_dir+'registration/cosine_subject_txt_3feat/')

# labels without zero signal vertices
zerovertexlh=np.load('../../Deliveries/0verticeslh.npy')
zerovertexrh=np.load('../../Deliveries/0verticesrh.npy')
lhannot=nib.freesurfer.io.read_annot(lhparpath)
lhlabels=lhannot[0]
rhannot=nib.freesurfer.io.read_annot(rhparpath)
rhlabels=rhannot[0]
labelslh=np.delete(lhlabels,zerovertexlh,0)
labelsrh=np.delete(rhlabels,zerovertexrh,0)


# assign rois to networks
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

#  registration
realigned=pickle.load(open(tar_dir+'alignedembeddings_cosine.p','rb'))
refemb=realigned[0][0]
refcoord=np.transpose(refemb)
lhrefout=tar_dir+'registration/cosine_subject_txt_3feat/lhreference.txt'
rhrefout=tar_dir+'registration/cosine_subject_txt_3feat/rhreference.txt'
lhrefcoord=refcoord[:,0:len(labelslh)]
rhrefcoord=refcoord[:,len(labelslh):4795]
lhrefcoord=myinsert(lhrefcoord,zerovertexlh,-100,1)#insert former removed vertices without signal with costumized function
rhrefcoord=myinsert(rhrefcoord,zerovertexrh,-100,1)
lhrefcoord=lhrefcoord[0:3,:]
rhrefcoord=rhrefcoord[0:3,:]
np.savetxt(lhrefout,lhrefcoord)
np.savetxt(rhrefout,rhrefcoord)

for i in range(231):#231
     souremb=realigned[0][i+1]
     sourcoord=np.transpose(souremb)
     lhsourout=tar_dir+'registration/cosine_subject_txt_3feat/lh'+str(i)+'.txt'
     rhsourout=tar_dir+'registration/cosine_subject_txt_3feat/rh'+str(i)+'.txt'
     lhsourcoord=sourcoord[:,0:len(labelslh)]
     rhsourcoord=sourcoord[:,len(labelslh):4795]
     lhsourcoord=myinsert(lhsourcoord,zerovertexlh,-100,axis=1)
     rhsourcoord=myinsert(rhsourcoord,zerovertexrh,-100,axis=1)
     lhsourcoord=lhsourcoord[0:3,:]
     rhsourcoord=rhsourcoord[0:3,:]
     np.savetxt(lhsourout,lhsourcoord)
     np.savetxt(rhsourout,rhsourcoord)
     if not os.path.exists(tar_dir+'registration/cosine_result_3feat/'+str(i)):
                os.makedirs(tar_dir+'registration/cosine_result_3feat/'+str(i))
     lhoutput=tar_dir+'registration/cosine_result_3feat/'+str(i)+'/L.'
     rhoutput=tar_dir+'registration/cosine_result_3feat/'+str(i)+'/R.'
     os.system('./msm_centos --inmesh='+lhmesh+' --indata='+lhsourout+' --refdata='+lhrefout+' -o '+lhoutput+' -f VTK')
     os.system('./msm_centos --inmesh='+rhmesh+' --indata='+rhsourout+' --refdata='+rhrefout+' -o '+rhoutput+' -f VTK')
