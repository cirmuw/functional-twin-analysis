#script to calculate embeddings and align them (rotation only)
#input: correlation matrices of refrence and twins
#output: aligned embeddings
import numpy as np
import nibabel as nib
import pandas as pd
import sklearn.metrics
from glob import glob
import os
import mapalign
import pickle
from mapalign import align
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import settings as s




#paths
lhparpath='../../Deliveries/lh.Schaefer2018_600Parcels_7Networks_order.annot'
rhparpath='../../Deliveries/rh.Schaefer2018_600Parcels_7Networks_order.annot'
root_dir=s.projectfolder+'long_corr_matricesLR/'
tar_dir=s.projectfolder+'7NETS_vertex/5_7nets_embedding/'
if not os.path.exists(tar_dir):
   os.mkdir(tar_dir)


files=[root_dir+str(i)+'.npy' for i in range(231)]
files.append(root_dir+'reference.npy')

# calulate embeddings and align them

reference=np.load(files[231],mmap_mode='r')
ref=2-sklearn.metrics.pairwise_distances(reference,metric='cosine')
 
refembedding,refresults=mapalign.embed.compute_diffusion_map(ref, alpha=0.5, n_components=None, diffusion_time=0,
                          skip_checks=False, overwrite=False,
                          eigen_solver=None, return_result=True)
  
  

  
embeddings=[refembedding]
for i in range(231):
    mat=np.load(files[i],mmap_mode='r')
    mat=2-sklearn.metrics.pairwise_distances(mat,metric='cosine')
    embedding,results=mapalign.embed.compute_diffusion_map(mat, alpha=0.5, n_components=refembedding.shape[1], diffusion_time=0,
                          skip_checks=False, overwrite=False,
                          eigen_solver=None, return_result=True)
    embeddings.append(embedding)
   

realigned=align.iterative_alignment(embeddings)
pickle.dump(realigned,open(tar_dir+'alignedembeddings_cosine.p','wb'))
  



