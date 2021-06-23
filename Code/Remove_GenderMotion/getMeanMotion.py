import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s

infile ='../../Deliveries/Subjects.csv' # orderd subjects
subjectpath1=s.HCProot+'/HCP_3T_RESTA_Motion/motion_mean/'# paths to relative mean motion of subjects
subjectpath2=s.HCProot+'/HCP_3T_RESTB_Motion/motion_mean/'#/
targetdir='../../Deliveries/'

df=pd.read_csv(infile)
Subjects=df['Subject'].values






path1=[]
path2=[]
for i in range(Subjects.shape[0]):
    path1.append(subjectpath1+str(Subjects[i]))
    path2.append(subjectpath2+str(Subjects[i]))


truesubjects=[]

for i in range(Subjects.shape[0]):
    if os.path.isdir(path1[i])==True:
        truesubjects.append(path1[i]+'/REST1_LR_Movement_RelativeRMS_mean.txt')#chose run from which to get motion
    if os.path.isdir(path2[i])==True and os.path.isdir(path1[i])==False:
        truesubjects.append(path2[i]+'/REST1_LR_Movement_RelativeRMS_mean.txt')#chose run from which to get motion
   
store=[]
for i in range(len(truesubjects)):
    motion=pd.read_csv(truesubjects[i],sep=" ",header=None)
    store.append(motion)
motiontable=pd.concat(store,axis=0,ignore_index=True)
motiontable.columns=['Rel_Mean_Motion']# Rel_Mean_Motion for 1st run, Rel_Mean_Motion2
table=pd.concat([df['Subject'],motiontable],axis=1)
writefile=targetdir+'MeanMotionRun1LR.csv.gz' #MeanMotionRun1LR.csv.gz for run 1, MeanMotionRun2LR.csv.gz for run 2
table.to_csv(writefile, compression='gzip')
