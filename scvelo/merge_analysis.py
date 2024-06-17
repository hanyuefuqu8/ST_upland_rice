import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import anndata as ad
import os
import sys
from scipy import sparse as sp

path= '../dynamo/'
all_data=[]
lines=['GS69','GS74','GS106','GS112']
name10=[]
for line in open("sections.txt"):
    line=line.strip('\n')
    temp=str.split(line,'/')
    n=len(str.split(line,'/'))
    temp0=temp[n-1]
    temp1=temp0.replace('.gem','')
    name10.append(temp1)

name20=[x.replace('_','-') for x in name10]

for line in lines:
    name1=[x for x in name10 if (line in x)]
    name2=[x for x in name20 if (line in x)]
    all_data=[]
    for index in range(len(name1)-1):
        loom = path + name1[index] + '/rna_velocity.loom'
        data_path = path + name1[index] +'/out.h5ad'
        meta_path = path+'meta_path/BIN40'+name2[index] +'/'
        adata = scv.read(loom, cache=True)
        sdata=scv.read(data_path)
        adata.obs['x']=sdata.obs['x']
        adata.obs['y']=sdata.obs['y']
        adata.obsm['X_spatial']=adata.obs[['x','y']].values.astype(float)
        xs=adata.obs['x']/40
        ys=adata.obs['y']/40
        id='BIN40'+ name2[index] +'_X'+ xs.astype("int").astype("str") + '_' + ys.astype("int").astype("str")
        adata.obs=adata.obs.rename(id)
        sample_obs = pd.read_csv(os.path.join(meta_path, "cellID_obs.csv"))
        cell_umap= pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"), header=0, names=["CellID", "UMAP_1", "UMAP_2"])
        cell_clusters = pd.read_csv(os.path.join(meta_path, "cell_clusters.csv"), header=0, names=["CellID", "cluster"])
        sample_one = adata[np.isin(adata.obs.index, sample_obs)]
        sample_one_index = pd.DataFrame(sample_one.obs.index)
        sample_one_index = sample_one_index.rename(columns = {0:'CellID'})
        clusters_ordered = sample_one_index.merge(cell_clusters, on = "CellID")
        umap_ordered = sample_one_index.merge(cell_umap, on = "CellID")
        umap_ordered = umap_ordered.iloc[:,1:]
        clusters_ordered = clusters_ordered.iloc[:,1:]
        sample_one.obsm['X_umap'] = umap_ordered.values
        sample_one.obs['clusters'] = clusters_ordered.values
        adata = sample_one
        adata.obs['group']=name2[index]
        all_data.append(adata)

    adata = ad.concat(all_data, join="outer")

    scv.settings.figdir ='fig'
    scv.set_figure_params('scvelo')

    scv.pp.filter_and_normalize(adata, n_top_genes=3000)
    CR_select = pd.read_csv("../CR_select/CR_selected.csv",header=0,index_col=0)
    adatas = adata[adata.obs.index.isin(CR_select.index)]
    adata = adatas
    ident_colours=["#CB181D","#c49374", "#EAFC1C","#FB6A4A","#FCAE91", "#FEE5D9"]
    #Dynamical Modeling
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    #scv.pl.proportions(adata,save=line+'_proportion',groupby='clusters')

    scv.pl.velocity_embedding_stream(adata, basis='X_umap',size=300,color = "clusters", palette = ident_colours,save=line+"_UMAP",density=4,legend_loc='right')

    adatas = []
    for group, idx in adata.obs.groupby("group").indices.items():
        sub_adata = adata[idx].copy()
        scv.pl.velocity_embedding_stream(sub_adata, basis='spatial',color = "clusters", palette = ident_colours,size = 300,alpha =0.7,save=group+"_spatial.png",density=4,legend_loc='right')
        adatas.append(sub_adata)
