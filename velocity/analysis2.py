# The first parameter is prefix of cell_ID, the second parameter is meta path
import scvelo as scv
import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import anndata as ad
import os
import sys
from scipy import sparse as sp
adata = scv.read('1.h5ad', cache=True)
name1 = sys.argv[1]
name2 = sys.argv[2]
#adata.obsm=adata.obs[['cx','cy']]

scv.settings.figdir ='fig'
scv.set_figure_params('scvelo')
#scv.pl.proportions(adata)
scv.pp.filter_and_normalize(adata, n_top_genes=3000)
adata.obsm['X_spatial']=adata.obs[['cx','cy']].values.astype(float)

#将seurat分群信息导入
#将细胞名修改为与Seurat一致
adata.obs=adata.obs.rename(index = lambda x: name1+x)
meta_path = name2
sample_obs = pd.read_csv(os.path.join(meta_path, "cellID_obs.csv"))
cell_umap= pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"), header=0, names=["Cell ID", "UMAP_1", "UMAP_2"])
cell_clusters = pd.read_csv(os.path.join(meta_path, "cell_clusters.csv"), header=0, names=["Cell ID", "cluster"])
#取交集
sample_one = adata[np.isin(adata.obs.index, sample_obs)]
sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'Cell ID'})
clusters_ordered = sample_one_index.merge(cell_clusters, on = "Cell ID")
umap_ordered = sample_one_index.merge(cell_umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
clusters_ordered = clusters_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
sample_one.obs['clusters'] = clusters_ordered.values

adata = sample_one
#scv.pp.normalize_per_cell(adata)
#scv.pp.moments(adata)

ident_colours=["#E64B35","#74C476","#3C5488","#8491B4","#DC0000","#4DBBD5","#00441B","#7E6148"]
#Dynamical Modeling
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='X_umap',size=120,color = "clusters", palette = ident_colours,save="_UMAP")
scv.pl.velocity_embedding_stream(adata, basis='spatial',color = "clusters", palette = ident_colours,size = 120,alpha =0.7,save="_spatial")

#计算每个基因的剪接动力
#df = adata.var
#df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
#scv.get_df(adata, 'fit*', dropna=True).head()

#时序变化基因
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,save='_dynamical_latent')
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100,save='_top_genes')
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False,save='_top_gene_fit')

scv.tl.rank_dynamical_genes(adata, groupby='clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)

#for cluster in ['Crown_root_primordia', 'Parenchyma_1', 'Xylem_1', 'Xylem_2']:
#   scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, frameon=False,save=cluster)


scv.pl.proportions(adata,save='_proportion',groupby='clusters')


#ident_colours=["#E64B35","#74C476","#3C5488","#8491B4","#DC0000","#4DBBD5","#00441B","#7E6148"]
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8)
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8, save="UMAP_stream.png", figsize=(7,5), legend_fontsize = 9, show=False, title='')
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8, save="UMAP_stream.pdf", figsize=(7,5), legend_fontsize = 9, show=False, title='')



#identify important genes
#scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])

#kwargs = dict(frameon=True,add_outline='Crown_root_primordia, Pericycle_like, Parenchyma_1,Parenchyma_2')
#scv.pl.scatter(adata, df['Parenchyma_2'][:5], ylabel='Parenchyma_2', **kwargs,save="_important_genes")


#Speed and coherence
#scv.tl.velocity_confidence(adata)
#keys = 'velocity_length', 'velocity_confidence'
#scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save="_speed")

#pseudotime
#scv.pl.velocity_graph(adata, threshold=.1,save='_connection')

#x, y = scv.utils.get_cell_transitions(adata, basis='spatial', starting_cell='BIN40GS112-1-1_X249_354')
#ax = scv.pl.velocity_graph(adata, basis='spatial',c='lightgrey', edge_width=.05, show=False)
#ax = scv.pl.scatter(adata, basis='spatial',x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax,save='_trajectory_of_one_cell')

#scv.tl.velocity_pseudotime(adata)
#scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',save='_trajectory')
#scv.pl.scatter(adata, basis='spatial',color='velocity_pseudotime', cmap='gnuplot',save='_spatial_trajectory')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='clusters')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,save="_page")




#save data
adata.write('Allcelltype_dynamicModel.h5ad', compression = 'gzip')


