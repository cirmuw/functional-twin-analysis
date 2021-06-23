#script to calculate for each vertex the correlation of the vertex time series with time series of selected parcels
#input:fmri series of twin subjects
#output:tabels containig correlation values for each twin, one table per vertex
        #correlation matrices for each subjects as intermediate output
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
from myfuns import subject_into_csv_vertex


WhichModelType='Original'   #'Original' # setting for antomical aligned data ('Original') or functional aligned data ('AfterEmbedding')
corr_mats_calculated=False     #set if correlatin matrix of time series have aleady been calculated


#paths to fmri data, parcellation, id of vetices without signal,chosen rois
infile =s.HCP_information_sheet_path         # \
subjectpath1=s.HCProot+'HCP_3T_RESTA_fmri/'# used to get paths to fmri runs of twin subjects
subjectpath2=s.HCProot+'HCP_3T_RESTB_fmri/'# /

zerovertexlh=np.load('../../Deliveries/0verticeslh.npy')# id of vertices without signal
zerovertexrh=np.load('../../Deliveries/0verticesrh.npy')

lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'# path to parcelation
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'

lhrois=list(np.load('../../Deliveries/chosenroislh.npy'))
rhrois=list(np.load('../../Deliveries/chosenroisrh.npy'))
lhrois=lhrois[1:]
rhrois=rhrois[1:]
nameslhrois=['l_'+str(s) for s in lhrois]
namesrhrois=['r_'+str(s) for s in rhrois]



# set additional input and output paths
if WhichModelType=='AfterEmbedding': 
   source_dir=s.projectfolder+'7NETS_vertex/5_7nets_corresponding/' # path containig id of functionally corresponding vertices of each twin to reference
   corr_mat_dir=s.projectfolder+'7NETS_vertex/6_7nets_corr_matrices_cosine/'# correlation values for each subject, intermediate output
   targetdir=s.projectfolder+'7NETS_vertex/7_7nets_ROIs_cosine/'# tables for each vertex
   if not os.path.exists(corr_mat_dir):
     os.mkdir(corr_mat_dir)
   if not os.path.exists(targetdir):
     os.mkdir(targetdir)

 
if WhichModelType=='Original':
   corr_mat_dir=s.projectfolder+'7NETS_vertex/1_7net_corr_matrices/'# correlation values for each subject, intermediate output
   targetdir=s.projectfolder+'7NETS_vertex/2_7nets_ROIs/'# tables for each vertex
   if not os.path.exists(corr_mat_dir):
     os.mkdir(corr_mat_dir)
   if not os.path.exists(targetdir):
     os.mkdir(targetdir)

   
lhannot=nib.freesurfer.io.read_annot(lhparpath)
lhlabels=lhannot[0]
rhannot=nib.freesurfer.io.read_annot(rhparpath)
rhlabels=rhannot[0]
labelslh=np.delete(lhlabels,zerovertexlh,0)
labelsrh=np.delete(rhlabels,zerovertexrh,0)


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

      

   

#Get paths to mgh-files(fmri runs) of available twin subjects
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

fmri_LH_LR_R1=[]
fmri_RH_LR_R1=[]
fmri_LH_RL_R1=[]
fmri_RH_RL_R1=[]
fmri_LH_LR_R2=[]
fmri_RH_LR_R2=[]
fmri_LH_RL_R2=[]
fmri_RH_RL_R2=[]
truesubjects=[]

for i in range(Subjects.shape[0]):
    if os.path.isdir(path1[i])==True:
        fmri_LH_LR_R1.append(path1[i]+'/lh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R1.append(path1[i]+'/rh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_LH_RL_R1.append(path1[i]+'/lh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_RH_RL_R1.append(path1[i]+'/rh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_LR_R2.append(path1[i]+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R2.append(path1[i]+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_LH_RL_R2.append(path1[i]+'/lh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_RH_RL_R2.append(path1[i]+'/rh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        truesubjects.append(Subjects[i])
    if os.path.isdir(path2[i])==True:
        fmri_LH_LR_R1.append(path2[i]+'/lh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R1.append(path2[i]+'/rh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_LH_RL_R1.append(path2[i]+'/lh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_RH_RL_R1.append(path2[i]+'/rh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_LR_R2.append(path2[i]+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R2.append(path2[i]+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_LH_RL_R2.append(path2[i]+'/lh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        #fmri_RH_RL_R2.append(path2[i]+'/rh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        truesubjects.append(Subjects[i])








#calculate correlation matrices for subjects
if corr_mats_calculated==False:
 
  for i in range(len(truesubjects)):
      
       lh1=nib.load(fmri_LH_LR_R1[i]).get_data() #choose LR or RL runs
       lh2=nib.load(fmri_LH_LR_R2[i]).get_data()
       rh1=nib.load(fmri_RH_LR_R1[i]).get_data()
       rh2=nib.load(fmri_RH_LR_R2[i]).get_data()

       lh1.resize(lh1.shape[1],lh1.shape[3])
       lh2.resize(lh2.shape[1],lh2.shape[3])
       rh1.resize(rh1.shape[1],rh1.shape[3])
       rh2.resize(rh2.shape[1],rh2.shape[3])

       lh=np.concatenate((lh1,lh2),axis=1)
       rh=np.concatenate((rh1,rh2),axis=1)
       lh=np.delete(lh,zerovertexlh,axis=0)
       rh=np.delete(rh,zerovertexrh,axis=0)
       
       if WhichModelType=='Original': 
          lhframe=pd.DataFrame(np.transpose(lh))
          lhframe=lhframe.groupby(labelslh,axis=1).mean()
          rhframe=pd.DataFrame(np.transpose(rh))
          rhframe=rhframe.groupby(labelsrh,axis=1).mean()

      
          lhframe=lhframe.iloc[:,lhrois]
          lhframe.columns=nameslhrois
          rhframe=rhframe.iloc[:,rhrois]
          rhframe.columns=namesrhrois
          frame=pd.concat([lhframe,rhframe],axis=1)
          for j in range(len(labelslh)):
           if labelslh[j]!=0:
              if not os.path.exists(corr_mat_dir+'lh_vertex_'+str(j+1)):
                 os.makedirs(corr_mat_dir+'lh_vertex_'+str(j+1))
              vertex=pd.DataFrame(np.transpose(lh[j,:]))
              vertex.columns=[str(j+1)]
              assign=lhnetwork>labelslh[j]
              assign=np.expand_dims(assign,0)
              assign=np.where(np.any(assign,axis=0))[0][0]
              assign=assign-2
              frame=frame.drop(nameslhrois[assign],1)
              frame=pd.concat([vertex,frame],1)
              framecor=frame.corr()
              outfileLR=corr_mat_dir +'lh_vertex_'+str(j+1)+'/'+ str(i) + '.csv.gz'
              framecor.to_csv(outfileLR, compression='gzip')
              frame=pd.concat([lhframe,rhframe],axis=1)

  
          for j in range(len(labelsrh)):
           if labelsrh[j]!=0:
              if not os.path.exists(corr_mat_dir+'rh_vertex_'+str(j+1)):
                 os.makedirs(corr_mat_dir+'rh_vertex_'+str(j+1))
              vertex=pd.DataFrame(np.transpose(rh[j,:]))
              vertex.columns=[str(j+1)]
              assign=rhnetwork>labelsrh[j]
              assign=np.expand_dims(assign,0)
              assign=np.where(np.any(assign,axis=0))[0][0]
              assign=assign-2
              frame=frame.drop(namesrhrois[assign],1)
              frame=pd.concat([vertex,frame],1)
              framecor=frame.corr()
              outfileLR=corr_mat_dir +'rh_vertex_'+str(j+1)+'/'+ str(i) + '.csv.gz'
              framecor.to_csv(outfileLR, compression='gzip')
              frame=pd.concat([lhframe,rhframe],axis=1)
               

       if WhichModelType=='AfterEmbedding':  
          for j in range(len(labelslh)): 
           if labelslh[j]!=0:
             if not os.path.exists(corr_mat_dir+'lh_vertex_'+str(j+1)):
                 os.makedirs(corr_mat_dir+'lh_vertex_'+str(j+1))
             functional=pickle.load(open(source_dir+'lh_'+str(j+1)+'correspondingvertices.p','rb')) 
             brain=np.concatenate([lh,rh],axis=0)
             brain=brain[np.squeeze(functional[0][i]),:]   
             brainframe=pd.DataFrame(np.transpose(brain))
             brainframe=brainframe.groupby(functional[1],axis=1).mean()

             assign=lhnetwork>labelslh[j]
             assign=np.expand_dims(assign,0)
             assign=np.where(np.any(assign,axis=0))[0][0]
             assign=assign-2
             brainframe.columns=[str(j+1)]+nameslhrois[0:assign]+nameslhrois[assign+1:7]+namesrhrois
             framecor=brainframe.corr()

             outfileLR=corr_mat_dir +'lh_vertex_'+str(j+1)+'/'+ str(i) + '.csv.gz'
             framecor.to_csv(outfileLR, compression='gzip')

          for j in range(len(labelsrh)): 
           if labelsrh[j]!=0:
             if not os.path.exists(corr_mat_dir+'rh_vertex_'+str(j+1)):
                 os.makedirs(corr_mat_dir+'rh_vertex_'+str(j+1))
             functional=pickle.load(open(source_dir+'rh_'+str(j+1)+'correspondingvertices.p','rb')) 
             brain=np.concatenate([lh,rh],axis=0)
             brain=brain[np.squeeze(functional[0][i]),:]   
             brainframe=pd.DataFrame(np.transpose(brain))
             brainframe=brainframe.groupby(functional[1],axis=1).mean()

             assign=rhnetwork>labelsrh[j]
             assign=np.expand_dims(assign,0)
             assign=np.where(np.any(assign,axis=0))[0][0]
             assign=assign-2
             brainframe.columns=[str(j+1)]+nameslhrois+namesrhrois[0:assign]+namesrhrois[assign+1:7]
             framecor=brainframe.corr()

             outfileLR=corr_mat_dir +'rh_vertex_'+str(j+1)+'/'+ str(i) + '.csv.gz'
             framecor.to_csv(outfileLR, compression='gzip')
  print('correlation matrices: done')         









# create tables for each vertex based on correlation matrices
if WhichModelType=='Original': 
#lh
  subject_into_csv_vertex(len(labelslh),'lh',truesubjects,WhichModelType,labelslh,lhnetwork,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=None)
#rh
  subject_into_csv_vertex(len(labelsrh),'rh',truesubjects,WhichModelType,labelsrh,rhnetwork,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=None)



if WhichModelType=='AfterEmbedding':
#lh
  subject_into_csv_vertex(len(labelslh),'lh',truesubjects,WhichModelType,labelslh,lhnetwork,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=source_dir)
#rh
  subject_into_csv_vertex(len(labelsrh),'rh',truesubjects,WhichModelType,labelsrh,rhnetwork,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=source_dir)



print('Finished') 































