# find vertex to vertex correspondance after functional alignment between referene subject and twin subjects
#input:MSM output
#output: for each vertex of reference subject the id of the functinally corresponding vertix in twins plus the functionally corresponding vertices of each twin
#belonging to to the selected parcels
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
lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'#parcelation
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'
lhmesh='../../Deliveries/fsaverage4/lh.sphere.vtk'#surface file for finding correspondance of vertices
rhmesh='../../Deliveries/fsaverage4/rh.sphere.vtk'
source_dir=s.projectfolder+'7NETS_vertex/5_7nets_embedding/'#root for aligned embeddings
tar_dir=s.projectfolder+'/7NETS_vertex/5_7nets_corresponding/'# output directory
if not os.path.exists(tar_dir):
  os.mkdir(tar_dir)


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



# get corresponding vertices
reader=vtk.vtkPolyDataReader()
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()

reader.SetFileName(lhmesh)
reader.Update()
refpoly=reader.GetOutput().GetPoints().GetData()
refpointslh=VN.vtk_to_numpy(refpoly)
refpointslh=np.delete(refpointslh,zerovertexlh,0)

reader.SetFileName(rhmesh)
reader.Update()
refpoly=reader.GetOutput().GetPoints().GetData()
refpointsrh=VN.vtk_to_numpy(refpoly)
refpointsrh=np.delete(refpointsrh,zerovertexrh,0)



realigned=pickle.load(open(source_dir+'alignedembeddings_cosine.p','rb'))
corresponding_lh=[]
corresponding_rh=[]
for i in range(231):
  filename=source_dir+'registration/cosine_result_3feat/'+str(i)+'/L.sphere.reg.vtk'
  reader.SetFileName(filename)
  reader.Update()
  sourpoly=reader.GetOutput().GetPoints().GetData()
  sourpointslh=VN.vtk_to_numpy(sourpoly)
  sourpointslh=np.delete(sourpointslh,zerovertexlh,0)
  lh_nn=NearestNeighbors(n_neighbors=1,algorithm ='ball_tree',metric='euclidean').fit(sourpointslh)
  lh_distances,lh_indices=lh_nn.kneighbors(refpointslh)#ID of functionally corresponding vertices
  corresponding_lh.append(lh_indices)
   
  filename=source_dir+'registration/cosine_result_3feat/'+str(i)+'/R.sphere.reg.vtk'
  reader.SetFileName(filename)
  reader.Update()
  sourpoly=reader.GetOutput().GetPoints().GetData()
  sourpointsrh=VN.vtk_to_numpy(sourpoly)
  sourpointsrh=np.delete(sourpointsrh,zerovertexrh,0)
  rh_nn=NearestNeighbors(n_neighbors=1,algorithm ='ball_tree',metric='euclidean').fit(sourpointsrh)
  rh_distances,rh_indices=rh_nn.kneighbors(refpointsrh)
  corresponding_rh.append(rh_indices)










#chose needed vertices
lhrois=list(np.load('../../Deliveries/chosenroislh.npy'))
rhrois=list(np.load('../../Deliveries/chosenroisrh.npy'))
lhrois=lhrois[1:]
rhrois=rhrois[1:]
for j in range(len(labelslh)):
 if labelslh[j]!=0:
   
   assign=lhnetwork>labelslh[j]
   assign=np.expand_dims(assign,0)
   assign=np.where(np.any(assign,axis=0))[0][0]
   assign=assign-2
   tmplhrois=np.delete(lhrois,assign)   
   neededlabelslh_index=np.any((labelslh==tmplhrois[0], labelslh==tmplhrois[1],labelslh==tmplhrois[2],labelslh==tmplhrois[3],labelslh==tmplhrois[4],labelslh==tmplhrois[5]),axis=0)
   neededlabelslh_index[j]=True
   neededlabelsrh_index=np.any((labelsrh==rhrois[0], labelsrh==rhrois[1], labelsrh==rhrois[2], labelsrh==rhrois[3], labelsrh==rhrois[4],labelsrh==rhrois[5], labelsrh==rhrois[6]),axis=0)
   tmp=labelslh[j]
   labelslh[j]=-1
   labelslh_rd=labelslh[neededlabelslh_index]
   labelsrh_rd=labelsrh[neededlabelsrh_index]+1000
   labels=np.concatenate([labelslh_rd,labelsrh_rd],axis=0)
   labelslh[j]=tmp
  

   chosenvertices=[]
   for i in range(231):
       lh_indices=corresponding_lh[i][neededlabelslh_index]
       rh_indices=corresponding_rh[i][neededlabelsrh_index]+len(labelslh)
       ind=np.concatenate([lh_indices,rh_indices],axis=0)
       chosenvertices.append(ind)  #ID of functionally corresponding vertices
      
       souremb=realigned[0][i+1]
       sourcoord=np.transpose(souremb)
       functionalcoord=sourcoord[:,ind]
       
       

   mylist=[chosenvertices,labels]
   pickle.dump(mylist,open(tar_dir+'lh_'+str(j+1)+'correspondingvertices.p','wb'))
 



for j in range(len(labelsrh)):
 if labelsrh[j]!=0:
   
   assign=rhnetwork>labelsrh[j]
   assign=np.expand_dims(assign,0)
   assign=np.where(np.any(assign,axis=0))[0][0]
   assign=assign-2
   tmprhrois=np.delete(rhrois,assign)
   neededlabelsrh_index=np.any((labelsrh==tmprhrois[0], labelsrh==tmprhrois[1],labelsrh==tmprhrois[2],labelsrh==tmprhrois[3],labelsrh==tmprhrois[4],labelsrh==tmprhrois[5]),axis=0)
   neededlabelsrh_index[j]=True
   neededlabelslh_index=np.any((labelslh==lhrois[0], labelslh==lhrois[1], labelslh==lhrois[2], labelslh==lhrois[3], labelslh==lhrois[4],labelslh==lhrois[5], labelslh==lhrois[6]),axis=0)
   tmp=labelsrh[j]
   labelsrh[j]=-1001
   labelslh_rd=labelslh[neededlabelslh_index]
   labelsrh_rd=labelsrh[neededlabelsrh_index]+1000
   labels=np.concatenate([labelslh_rd,labelsrh_rd],axis=0)
   labelsrh[j]=tmp  



   chosenvertices=[]
   for i in range(231):
       lh_indices=corresponding_lh[i][neededlabelslh_index]
       rh_indices=corresponding_rh[i][neededlabelsrh_index]+len(labelslh)
       ind=np.concatenate([lh_indices,rh_indices],axis=0)
       chosenvertices.append(ind)  #ID of functionally corresponding vertices
       
       souremb=realigned[0][i+1]
       sourcoord=np.transpose(souremb)
       functionalcoord=sourcoord[:,ind]
       

   mylist=[chosenvertices,labels]
   pickle.dump(mylist,open(tar_dir+'rh_'+str(j+1)+'correspondingvertices.p','wb'))
   


print('finshed')

