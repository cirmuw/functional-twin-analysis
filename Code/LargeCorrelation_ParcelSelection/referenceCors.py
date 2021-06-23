#script calculates and saves correlation matrix of reference subject for later use(embedding)
#input:frmi time series of reference subject
#output: correlation matrix of fmri time series

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



outfileLR=s.projectfolder + 'long_corr_matricesLR/reference'# output folder
outfileRL=s.projectfolder +'long_corr_matricesRL/reference'# output folder
reference_subjectpath=s.HCProot+'HCP_3T_RESTA_fmri/101915/' #folder containg fmri runs
files=glob(os.path.join(reference_subjectpath, "*.mgh"))

if not os.path.exists(s.projectfolder +'long_corr_matricesLR/'):
   os.mkdir(s.projectfolder +'long_corr_matricesLR/')
if not os.path.exists(s.projectfolder +'long_corr_matricesRL/'):
   os.mkdir(s.projectfolder +'long_corr_matricesRL/')


fmri_LH_LR_R1=reference_subjectpath+'/lh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_RH_LR_R1=reference_subjectpath+'/rh.rfMRI_REST1_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_LH_RL_R1=reference_subjectpath+'/lh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_RH_RL_R1=reference_subjectpath+'/rh.rfMRI_REST1_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_LH_LR_R2=reference_subjectpath+'/lh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_RH_LR_R2=reference_subjectpath+'/rh.rfMRI_REST2_LR_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_LH_RL_R2=reference_subjectpath+'/lh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
fmri_RH_RL_R2=reference_subjectpath+'/rh.rfMRI_REST2_RL_Atlas_hp2000_clean_bpss_gsr_fs4.mgh'
cor_LR=calc_vertex_correlations_combined_runs(fmri_LH_LR_R1,fmri_RH_LR_R1,fmri_LH_LR_R2,fmri_RH_LR_R2)
cor_RL=calc_vertex_correlations_combined_runs(fmri_LH_RL_R1,fmri_RH_RL_R1,fmri_LH_RL_R2,fmri_RH_RL_R2)

np.save(outfileLR,cor_LR['fullcor'])
np.save(outfileRL,cor_RL['fullcor'])
      
