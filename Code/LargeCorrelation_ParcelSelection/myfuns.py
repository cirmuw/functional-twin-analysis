#file containg customized functions in scripts
import numpy as np
import nibabel as nib
import pandas as pd
import pickle


def calc_vertex_correlations(path_lh, path_rh):
    #calculate correlation matrix of time series for each runs sepeately 

    # load fmri data for left and right hemi
    lh_fmri=nib.load(path_lh)
    rh_fmri=nib.load(path_rh)


    #get image data and resize
    lh_imagedata=lh_fmri.get_data()
    lh_imagedata.resize(lh_imagedata.shape[1],lh_imagedata.shape[3])
    rh_imagedata=rh_fmri.get_data()
    rh_imagedata.resize(rh_imagedata.shape[1],rh_imagedata.shape[3])





    # there are vetices without signal: check if thery are the only ones with zero variance
    lh_zerovertex=~lh_imagedata.any(axis=1)
    rh_zerovertex=~rh_imagedata.any(axis=1)

    lh_cov=np.cov(lh_imagedata)  
    lh_var=np.diag(lh_cov)
    lh_var=np.expand_dims(lh_var,axis=1)
    lh_varzero=~lh_var.any(axis=1)
    rh_cov=np.cov(rh_imagedata)
    rh_var=np.diag(rh_cov)
    rh_var=np.expand_dims(rh_var,axis=1)
    rh_varzero=~rh_var.any(axis=1)
    
    if np.all(lh_varzero==lh_zerovertex)==False or np.all(rh_varzero==rh_zerovertex)==False:
       print("\n###########################################################")
       print("All Vertices with zero variance of left hemi have no signal: ")
       print(np.all(lh_varzero==lh_zerovertex))
       print("All Vertices with zero variance of right hemi have no signal: ")
       print(np.all(rh_varzero==rh_zerovertex))
       print("############################################################\n")

    #calculate correlations and remove vertices without signal as seeds
    imagedata=np.concatenate((lh_imagedata,rh_imagedata),axis=0)
    zerovertex=np.concatenate((lh_zerovertex,rh_zerovertex),axis=0)
    deleteind=np.where(zerovertex)[0]
    imagedata=np.delete(imagedata,deleteind,0)
    cor=np.corrcoef(imagedata)
    
    
    
    return {'fullcor':cor,'zerovertices': deleteind, 'splitposition': lh_zerovertex.shape[0]}






















def calc_vertex_correlations_combined_runs (path_lh1, path_rh1,path_lh2, path_rh2):
   #calculate correlation matrix of time series after combining 2 runs

    # load fmri data for left and right hemi
    lh_fmri1=nib.load(path_lh1)
    rh_fmri1=nib.load(path_rh1)
    lh_fmri2=nib.load(path_lh2)
    rh_fmri2=nib.load(path_rh2)

    #get image data and resize
    lh_imagedata1=lh_fmri1.get_data()
    lh_imagedata1.resize(lh_imagedata1.shape[1],lh_imagedata1.shape[3])
    rh_imagedata1=rh_fmri1.get_data()
    rh_imagedata1.resize(rh_imagedata1.shape[1],rh_imagedata1.shape[3])
    lh_imagedata2=lh_fmri2.get_data()
    lh_imagedata2.resize(lh_imagedata2.shape[1],lh_imagedata2.shape[3])
    rh_imagedata2=rh_fmri2.get_data()
    rh_imagedata2.resize(rh_imagedata2.shape[1],rh_imagedata2.shape[3])




    # there are vetices without signal
    lh_zerovertex=~lh_imagedata1.any(axis=1)
    rh_zerovertex=~rh_imagedata1.any(axis=1)

   

    #calculate correlations and remove vertices without signal as seeds
    imagedata1=np.concatenate((lh_imagedata1,rh_imagedata1),axis=0)
    imagedata2=np.concatenate((lh_imagedata2,rh_imagedata2),axis=0)
    zerovertex=np.concatenate((lh_zerovertex,rh_zerovertex),axis=0)
    deleteind=np.where(zerovertex)[0]
    imagedata1=np.delete(imagedata1,deleteind,0)
    imagedata2=np.delete(imagedata2,deleteind,0) 
    imagedata=np.concatenate((imagedata1,imagedata2),axis=1)
    cor=np.corrcoef(imagedata)
    
    
    
    return {'fullcor':cor,'zerovertices': deleteind, 'splitposition': lh_zerovertex.shape[0]}








def subjects_into_csv_gz_7nets(nrOFNetworks,truesubjects,WhichModelType,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,targetdir2=None,functional=None,spatial=None):

  #For each ROI as seed create csv file containing connections to all other ROIs for all subjects
  labels=np.arange(nrOFNetworks) 
  labels=labels.tolist(labels)
  name=['Subject','Zygosity','Mother_ID']
  nonvertexdat=np.zeros((len(truesubjects),3),dtype=object)

  if WhichModelType=='AfterEmbedding':
       spatialvar=[]

  for j in range(len(labels)):
         tablelh=[]
         tablerh=[]
       
      
       
         for i in range(len(truesubjects)):
             infile=corr_mat_dir +str(i)+'.csv.gz'
             corr=pd.read_csv(infile,compression='gzip')
             tmp=corr[['Unnamed: 0',nameslhrois[j]]]
             tmp=tmp.drop([j])
             tmp=tmp.transpose()
             tmptablelh=tmp.drop(['Unnamed: 0'])
          
          
             tmp=corr[['Unnamed: 0',namesrhrois[j]]]
             tmp=tmp.drop([j+7])
             tmp=tmp.transpose()
             tmptablerh=tmp.drop(['Unnamed: 0'])

             tablelh.append(tmptablelh) 
             tablerh.append(tmptablerh) 
           


             if j==0:

               index=dataframe2[dataframe2['Subject']==truesubjects[i]].index.tolist()
               tmp1=np.array([str(truesubjects[i]),dataframe2['Zygosity'][index].values[0], str(dataframe2['Mother_ID'][index].values[0])])
               nonvertexdat[i,:]=tmp1
               if WhichModelType=='AfterEmbedding':
                  spatframe=pd.DataFrame(np.expand_dims(spatial[i],0)).groupby(functional[1],axis=1).mean()
                  spatframe.columns=nameslhrois+namesrhrois
                  spatialvar.append(spatframe)
                
                


         tablelh=pd.concat(tablelh,axis=0,ignore_index=True)
         tablerh=pd.concat(tablerh,axis=0,ignore_index=True)
         roinames=nameslhrois+namesrhrois
         del roinames[j]
         tablelh.columns=roinames
         roinames=nameslhrois+namesrhrois
         del roinames[j+7]
         tablerh.columns=roinames
         nonvertextable=pd.DataFrame(data=nonvertexdat)
         nonvertextable.columns=name
       
       
         table=pd.concat([nonvertextable,tablelh],1)
         table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
         table.reset_index(inplace=True)
         table=table.drop('index',axis=1)
         writefile=targetdir+'lh_net'+str(j+1)+'.csv.gz'
         table.to_csv(writefile, compression='gzip')
      
       
         table=pd.concat([nonvertextable,tablerh],1)
         table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0,ascending=[True,True])
         table.reset_index(inplace=True)
         table=table.drop('index',axis=1)
         writefile=targetdir+'rh_net'+str(j+1)+'.csv.gz'
         table.to_csv(writefile, compression='gzip')
    
  if WhichModelType=='AfterEmbedding':
     spatialframe=pd.concat(spatialvar,axis=0,ignore_index=True)
     table=pd.concat([nonvertextable,spatialframe],axis=1)
     table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
     table.reset_index(inplace=True)
     table=table.drop('index',axis=1)
     writefile=targetdir2+'mean_distances.csv.gz'
     table.to_csv(writefile, compression='gzip')






def subjects_into_csv_gz_specifiedROIs(ROIs,NetNr,hemi,truesubjects,WhichModelType,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=None,targetdir2=None):
#For each ROI as seed create csv file containing connections to all other ROIs for all subjects
  
  name=['Subject','Zygosity','Mother_ID']
  nonvertexdat=np.zeros((len(truesubjects),3),dtype=object)
  


  for j in range(ROIs[0],ROIs[1]):
         table=[]
         
         if WhichModelType=='AfterEmbedding_notSelected':
             spatialvar=[]         
         for i in range(len(truesubjects)):
             infile=corr_mat_dir+hemi+'_'+str(j)+'/' +str(i)+'.csv.gz'
             corr=pd.read_csv(infile,compression='gzip')
             
             if hemi=='lh':
                nameslhrois[NetNr-1]='lh_'+str(j)
                tmp=corr[['Unnamed: 0',nameslhrois[NetNr-1]]]
                tmp=tmp.drop([NetNr-1])
                
             if hemi=='rh':
                namesrhrois[NetNr-1]='rh_'+str(j)
                tmp=corr[['Unnamed: 0',namesrhrois[NetNr-1]]]
                tmp=tmp.drop([NetNr-1+7])
               
             
             tmp=tmp.transpose()
             tmptable=tmp.drop(['Unnamed: 0'])
             table.append(tmptable) 
             
             if j==ROIs[0]:
               index=dataframe2[dataframe2['Subject']==truesubjects[i]].index.tolist()
               tmp1=np.array([str(truesubjects[i]),dataframe2['Zygosity'][index].values[0], str(dataframe2['Mother_ID'][index].values[0])])
               nonvertexdat[i,:]=tmp1
             
             if WhichModelType=='AfterEmbedding_notSelected':
                  if hemi=='lh':
                     functional=pickle.load(open(source_dir+'lh_'+str(j)+'_'+'correspondingvertices.p','rb'))
                     spatial= pickle.load(open(source_dir+'lh_'+str(j)+'_'+'FuncSpatDifference.p','rb'))
                  if hemi=='rh':
                     functional=pickle.load(open(source_dir+'rh_'+str(j)+'_'+'correspondingvertices.p','rb'))
                     spatial= pickle.load(open(source_dir+'rh_'+str(j)+'_'+'FuncSpatDifference.p','rb'))
                  spatframe=pd.DataFrame(np.expand_dims(spatial[i],0)).groupby(functional[1],axis=1).mean()
                  spatframe.columns=nameslhrois+namesrhrois
                  spatialvar.append(spatframe)
                
                


         table=pd.concat(table,axis=0,ignore_index=True)
         roinames=nameslhrois+namesrhrois
         if hemi=='lh':
            del roinames[NetNr-1]
         if hemi=='rh':
            del roinames[NetNr-1+7]
         table.columns=roinames
         nonvertextable=pd.DataFrame(data=nonvertexdat)
         nonvertextable.columns=name
       
       
         table=pd.concat([nonvertextable,table],1)
         table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
         table.reset_index(inplace=True)
         table=table.drop('index',axis=1)
         if hemi=='lh':
            writefile=targetdir+'lh_roi'+str(j)+'.csv.gz'
         if hemi=='rh':
            writefile=targetdir+'rh_roi'+str(j)+'.csv.gz'
         table.to_csv(writefile, compression='gzip')
      
       
         if WhichModelType=='AfterEmbedding_notSelected':
            spatialframe=pd.concat(spatialvar,axis=0,ignore_index=True)
            table=pd.concat([nonvertextable,spatialframe],axis=1)
            table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
            table.reset_index(inplace=True)
            table=table.drop('index',axis=1)
            if hemi=='lh':
               writefile=targetdir2+'lh_'+str(j)+'_mean_distances.csv.gz'
            if hemi=='rh':
               writefile=targetdir2+'rh_'+str(j)+'_mean_distances.csv.gz'
            table.to_csv(writefile, compression='gzip')

def myinsert(array,position,value,axis):
     m,n=array.shape
     positionlength=len(position)
     position=np.concatenate((position,np.array([np.nan])),axis=0)
     counta=0
     countp=0
     if axis==0:
        newarray=np.zeros((m+positionlength,n),dtype=float)
        for i in range((m+positionlength)):
          if i==position[countp]:
            newarray[i,:]=value
            countp=countp+1
          else:
            newarray[i,:]=array[counta,:]
            counta=counta+1
          
     if axis==1:
        newarray=np.zeros((m,n+positionlength),dtype=float)
        for i in range((n+positionlength)):
          if i==position[countp]:
            newarray[:,i]=value
            countp=countp+1
          else:
            newarray[:,i]=array[:,counta]
            counta=counta+1
          
     return newarray

def subject_into_csv_vertex(nr_vertices,hemi,truesubjects,WhichModelType,labels,networkborders,corr_mat_dir,targetdir,nameslhrois,namesrhrois,dataframe2,source_dir=None,targetdir2=None):
  name=['Subject','Zygosity','Mother_ID']
  nonvertexdat=np.zeros((len(truesubjects),3),dtype=object)
  lhroisnames=np.array(nameslhrois)
  rhroisnames=np.array(namesrhrois)

  for j in range(nr_vertices):
    if labels[j]!=0: 
         table=[]
         
         if WhichModelType=='AfterEmbedding':
             spatialvar=[]         
         for i in range(len(truesubjects)):
           
             infile=corr_mat_dir+hemi+'_vertex_'+str(j+1)+'/' +str(i)+'.csv.gz'
             corr=pd.read_csv(infile,compression='gzip')
                                       
             tmp=corr[['Unnamed: 0',str(j+1)]]
             tmp=tmp.drop([0])             
             tmp=tmp.transpose()
             tmptable=tmp.drop(['Unnamed: 0'])
             table.append(tmptable) 
             
             assign=networkborders>labels[j]
             assign=np.expand_dims(assign,0)
             assign=np.where(np.any(assign,axis=0))[0][0]
             assign=assign-2
             if hemi=='lh':
                roinames=np.concatenate([np.delete(lhroisnames,assign),rhroisnames],0)
             if hemi=='rh':
                roinames=np.concatenate([lhroisnames,np.delete(rhroisnames,assign)],0)
             roinames=roinames.tolist()

             if j==0:
               index=dataframe2[dataframe2['Subject']==truesubjects[i]].index.tolist()
               tmp1=np.array([str(truesubjects[i]),dataframe2['Zygosity'][index].values[0], str(dataframe2['Mother_ID'][index].values[0])])
               nonvertexdat[i,:]=tmp1
             

           
             if WhichModelType=='AfterEmbedding':
                  if hemi=='lh':
                     functional=pickle.load(open(source_dir+'lh_'+str(j+1)+'correspondingvertices.p','rb'))
                     spatial= pickle.load(open(source_dir+'lh_'+str(j+1)+'FuncSpatDifference.p','rb'))
                  if hemi=='rh':
                     functional=pickle.load(open(source_dir+'rh_'+str(j+1)+'correspondingvertices.p','rb'))
                     spatial= pickle.load(open(source_dir+'rh_'+str(j+1)+'FuncSpatDifference.p','rb'))
                  spatframe=pd.DataFrame(np.expand_dims(spatial[i],0)).groupby(functional[1],axis=1).mean()
                  spatframe.columns=[str(j+1)]+roinames
                  spatialvar.append(spatframe)
                
                


         table=pd.concat(table,axis=0,ignore_index=True)
         table.columns=roinames
         nonvertextable=pd.DataFrame(data=nonvertexdat)
         nonvertextable.columns=name
       
       
         table=pd.concat([nonvertextable,table],1)
         table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
         table.reset_index(inplace=True)
         table=table.drop('index',axis=1)
         if hemi=='lh':
            writefile=targetdir+'lh_vertex'+str(j+1)+'.csv.gz'
         if hemi=='rh':
            writefile=targetdir+'rh_vertex'+str(j+1)+'.csv.gz'
         table.to_csv(writefile, compression='gzip')
      
       
         if WhichModelType=='AfterEmbedding':
            spatialframe=pd.concat(spatialvar,axis=0,ignore_index=True)
            table=pd.concat([nonvertextable,spatialframe],axis=1)
            table=table.sort_values(['Zygosity', 'Mother_ID'], axis=0, ascending=[True,True])
            table.reset_index(inplace=True)
            table=table.drop('index',axis=1)
            if hemi=='lh':
               writefile=targetdir2+'lh_'+str(j+1)+'_mean_distances.csv.gz'
            if hemi=='rh':
               writefile=targetdir2+'rh_'+str(j+1)+'_mean_distances.csv.gz'
            table.to_csv(writefile, compression='gzip')

            
            
            

def get_polygon_area_for_each_vertex(vertex_coords, polygons):
    triangle_area=np.zeros((polygons.shape[0]))
    for i in range(polygons.shape[0]):
        vertex_ind=polygons[i,1:4]
        v1=vertex_coords[vertex_ind[0],:]
        v2=vertex_coords[vertex_ind[1],:]
        v3=vertex_coords[vertex_ind[2],:]
        v_12=np.linalg.norm(v2-v1)
        v_13=np.linalg.norm(v3-v1)
        alpha=np.arccos(np.dot((v2-v1),(v3-v1))/(v_12*v_13))
        height=v_13*np.sin(alpha)
        triangle_area[i]=height*v_12/2
    polygon_area=np.zeros((vertex_coords.shape[0]))
    for i in range(vertex_coords.shape[0]):
        indices=np.any(polygons==i,axis=1)
        polygon_area[i]=np.sum(triangle_area[indices])
        
    return polygon_area
        
