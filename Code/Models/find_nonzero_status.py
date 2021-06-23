import numpy as np
import nibabel as nib
import pandas as pd
from glob import glob
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s


source_dir=s.projectfolder+'7NETS_vertex/4_estimatedModels' #entagled
#source_dir=s.projectfolder+'9_estimatedModels_functional_cosine'#connectivity strength
#source_dir=s.projectfolder+'12_PositionVar_estimatedModels' #position
files=glob(os.path.join(source_dir, "*.csv.gz"))

count=0
nonzero=[]
for i in range(len(files)):
   data=pd.read_csv(files[i],compression='gzip')
   status=data['Status Code'][0]
   if status!=0:
      #print(files[i])
      #print(status)
      count=count+1
      tmpfile=os.path.split(files[i])
      nonzero.append(tmpfile[1])
       
print(count)
Nonzero=np.array(nonzero,dtype=object)
Nonzero=np.expand_dims(Nonzero,axis=1)
table=pd.DataFrame(Nonzero)
table.columns=['nonzero']

writefile='/tmp_nonzero_anat.csv'
#writefile='tmp_nonzero_func.csv'
#writefile='tmp_nonzero_position.csv'
table.to_csv(writefile)
