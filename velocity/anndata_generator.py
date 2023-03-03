import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import anndata as ad
import os
from scipy import sparse as sp
import sys

loom_data = scv.read(sys.argv[1], cache=False)
loom_data.obs

#鏋勫缓鍚堝苟鐨刲oom鏂囦欢

#鎻愬彇鍧愭爣
B=int(sys.argv[2])
n=len(loom_data.obs)
loom_data.obs['X']=np.zeros((n,1),dtype = int)
loom_data.obs['Y']=np.zeros((n,1),dtype = int)
for i in range(0,n):
  loom_data.obs['X'][i]=int( int(loom_data.obs_names[i].split(":")[1].split("_")[0])/B)
  loom_data.obs['Y'][i]=int( int(loom_data.obs_names[i].split(":")[1].split("_")[1][:-1])/B)

loom_data.obs['ID']=loom_data.obs['X'].map(str)+"_"+loom_data.obs['Y'].map(str)

#鎻愬彇鍚堝苟鏁版嵁
d={}
n=0
for j in ['matrix', 'ambiguous', 'spliced', 'unspliced']:    
    df=pd.DataFrame()
    for i in set(loom_data.obs.ID):
        temp=loom_data[loom_data.obs.ID==i].to_df(layer=j)
        temp2=temp.apply(lambda x:x.sum())    
        temp2.name=i
        df=df.append(temp2)
    d[n]=df
    n+=1

#鏋勫缓鏂扮殑loom
obs = pd.DataFrame(d[0].index,index=d[0].index,columns=['cell_ID'])
loc = obs['cell_ID'].str.split('_', 1, expand=True).astype(int)
loc.columns = ['cx', 'cy']
obs = pd.merge(obs, loc, how='left', left_index=True, right_index=True)

loom_new=ad.AnnData(d[0],obs=obs,var=loom_data.var,layers={'matrix':d[0],'ambiguous':d[1],'spliced':d[2],'unspliced':d[3]})

loom_new.write("1.h5ad")


