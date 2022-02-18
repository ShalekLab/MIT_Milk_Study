import pandas as pd

def save_pseudobulk_counts(count_adata,processed_metadata=None,celltype_col=None,sample_col="sample"):
    '''
    count_adata should be an adata object with raw counts saved in the 'X' slot
    processed_metadata [a dataframe] is probably the obs of whatever adata object you have been working with, but if it is None, use the metadata from count_adata
    celltype_col is the column in the processed_metadata dataframe on which you would like to split your sample pseudobulks. If it is None, a single pseudobulk is created for each sample.
    sample_col is the other column to split your cells on (usually just called sample...)
    '''
    if processed_metadata is None:
        processed_metadata = count_adata.obs

    exp = pd.DataFrame(count_adata.X,index=count_adata.obs_names,columns=count_adata.var_names).loc[processed_metadata.index]
    exp[sample_col] = processed_metadata[sample_col]

    if celltype_col is not None:
        exp[celltype_col] = processed_metadata[celltype_col]
        pseudobulk = exp.groupby([celltype_col,sample_col]).sum()
        pseudobulk[str(celltype_col)+"_"+str(sample_col)] = [i[0]+"_"+i[1] for i in pseudobulk.index]
        pseudobulk.index = pseudobulk[str(celltype_col)+"_"+str(sample_col)]

        return pseudobulk[[c for c in pseudobulk.columns if c not in [sample_col,celltype_col,str(celltype_col)+"_"+str(sample_col)]]]
    else:
        return exp.groupby([sample_col]).sum()

def save_metadata_for_pseudobulk_deseq(processed_metadata, celltype_col = None, sample_col="sample",cell_counts=True):
    meta_cols = processed_metadata.columns
    if celltype_col is not None:
        grp_name = "TMP sample cols"
        processed_metadata[grp_name] = [processed_metadata.loc[i,celltype_col]+"_"+processed_metadata.loc[i,sample_col] for i in processed_metadata.index]
        df_index = list(processed_metadata[grp_name].unique())
    else:
        df_index = list(processed_metadata[sample_col].unique())
        grp_name = sample_col
    meta = pd.DataFrame(index=df_index, columns=meta_cols)

    for c in meta_cols:
        processed_metadata[c] = processed_metadata[c].astype(str)

        v = processed_metadata.groupby(grp_name)[c].unique()
        meta[c] = meta.index.map({i:v.loc[i][0] for i in v.index})
    if cell_counts:
        cnts = processed_metadata.groupby([grp_name]).count().iloc[:,0]
        meta["Cell_Number"]=cnts.loc[meta.index]
    return meta