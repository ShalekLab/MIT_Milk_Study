import anndata

def new_adata_from_raw(adata, cluster_values, cluster_key="leiden"):
    '''
    adata: object you want to subset
    cluster_values: a list of cluster names to include in the new raw subsetted adata
    cluster_key: the column in the metadata matrix (adata.obs) which contains the cluster_values on which you want to subset
    '''
    adata_0 = anndata.AnnData(adata[adata.obs[cluster_key].isin(cluster_values)].raw.X)
    adata_0.obs_names = adata[adata.obs[cluster_key].isin(cluster_values)].obs_names
    adata_0.var_names = adata[adata.obs[cluster_key].isin(cluster_values)].raw.var_names
    adata_0.obs = adata[adata.obs[cluster_key].isin(cluster_values)].obs
    adata_0.raw = adata_0
    return adata_0



def make_new_subcluster_adata(all_cells_df, cluster_key, clusters_to_include,old_adata, old_adata_name):
    new_adata= anndata.AnnData(all_cells_df[old_adata.obs.loc[old_adata.obs[cluster_key].isin(clusters_to_include)].index].T)
    new_adata.obs_names = old_adata.obs.loc[old_adata.obs[cluster_key].isin(clusters_to_include)].index
    new_adata.var_names = all_cells_df.index
    new_adata.obs = old_adata.obs.loc[old_adata.obs[cluster_key].isin(clusters_to_include)]
    new_adata.obs[old_adata_name+"_"+cluster_key] = new_adata.obs[cluster_key]
    return new_adata
