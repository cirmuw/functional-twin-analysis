# calculates correlation matrix of subjects
#input: fmri time series
#output: correlation matrix of subjects
import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import settings as s

from myfuns import calc_vertex_correlations_combined_runs
from myfuns import calc_vertex_correlations



infile =s.HCP_information_sheet_path #\
subjectpath1=s.HCProot+'HCP_3T_RESTA_fmri/'# needed to get path to fmri runs of twin subjects
subjectpath2=s.HCProot+'HCP_3T_RESTB_fmri/'#/
corr_mat_dir=s.projectfolder # project folder, output directory root

if not os.path.exists(corr_mat_dir +'long_corr_matricesLR/'):
   os.mkdir(corr_mat_dir +'long_corr_matricesLR/')
if not os.path.exists(corr_mat_dir +'long_corr_matricesRL/'):
   os.mkdir(corr_mat_dir +'long_corr_matricesRL/')


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
        fmri_LH_RL_R1.append(path1[i]+'/lh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_RL_R1.append(path1[i]+'/rh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_LR_R2.append(path1[i]+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R2.append(path1[i]+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_RL_R2.append(path1[i]+'/lh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_RL_R2.append(path1[i]+'/rh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        truesubjects.append(Subjects[i])
    if os.path.isdir(path2[i])==True:
        fmri_LH_LR_R1.append(path2[i]+'/lh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R1.append(path2[i]+'/rh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_RL_R1.append(path2[i]+'/lh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_RL_R1.append(path2[i]+'/rh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_LR_R2.append(path2[i]+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_LR_R2.append(path2[i]+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_LH_RL_R2.append(path2[i]+'/lh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        fmri_RH_RL_R2.append(path2[i]+'/rh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh')
        truesubjects.append(Subjects[i])




#calculate correlation matrices for subjects
for i in range(len(truesubjects)):
        
       # both runs for matrices
       cor_LR=calc_vertex_correlations_combined_runs(fmri_LH_LR_R1[i],fmri_RH_LR_R1[i],fmri_LH_LR_R2[i],fmri_RH_LR_R2[i])
       cor_RL=calc_vertex_correlations_combined_runs(fmri_LH_RL_R1[i],fmri_RH_RL_R1[i],fmri_LH_RL_R2[i],fmri_RH_RL_R2[i])

       outfileLR=corr_mat_dir +'long_corr_matricesLR/' + str(i)
       outfileRL=corr_mat_dir +'long_corr_matricesRL/'+ str(i)
       
       np.save(outfileLR,cor_LR['fullcor'])
       np.save(outfileRL,cor_RL['fullcor'])
          
print('correlation matrices: done')         





       
  
    
     
    















