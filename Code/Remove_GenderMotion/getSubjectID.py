import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s


infile =s.HCP_information_sheet_path #file with zygosity info
subjectpath1=s.HCProot+'HCP_3T_RESTA_fmri/'# paths to HCP fMRI runs
subjectpath2=s.HCProot+'HCP_3T_RESTB_fmri/'#


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



name=['Subject','Zygosity','Mother_ID','Handedness','Age_in_Yrs']
nonvertexdat=np.zeros((len(truesubjects),5),dtype=object)
for i in range(len(truesubjects)):
    index=dataframe2[dataframe2['Subject']==truesubjects[i]].index.tolist()
    tmp1=np.array([str(truesubjects[i]),dataframe2['Zygosity'][index].values[0], str(dataframe2['Mother_ID'][index].values[0]),str(dataframe2['Handedness'][index].values[0]),str(dataframe2['Age_in_Yrs'][index].values[0])])
    nonvertexdat[i,:]=tmp1

table=pd.DataFrame(data=nonvertexdat)
table.columns=name
table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
table.reset_index(inplace=True)
#table=table.drop('index',axis=1)
table.to_csv('../../Deliveries/Subjects.csv')


