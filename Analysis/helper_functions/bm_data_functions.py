import pandas as pd
def get_bm_metadata():
    #TODO add donor metadata
    return  pd.read_csv("../../Data/Supplemental Dataset 1 Metadata_BM_Study.csv", index_col=0, sep=",")

def bm_filtering(adata,timepoint):
    calc_qc_bm(adata,sample_key="sample")
    df = pd.DataFrame(index = adata.obs["sample"].unique())
    fig,ax = plt.subplots(2,3,figsize=(12,6))
    for i,stat in enumerate(['n_genes_by_counts', 'total_counts', "percent_mito"]):
        p = sc.pl.violin(adata, stat, groupby="donor",jitter=0.4, show=False,ax=ax[0,i])
        p=sc.pl.violin(adata, stat, groupby="sample",jitter=0.4, show=False,rotation=90,ax=ax[1,i])
    plt.tight_layout()
    fig.savefig("qc_plots_per_timepoint/"+timepoint+"_prefiltering.pdf",bbox_inches='tight')
    fig.savefig("qc_plots_per_timepoint/"+timepoint+"_prefiltering.png",bbox_inches='tight')
    df["prefilter_total_cells"] = adata.obs.groupby("sample").count()["donor"].loc[df.index]
    df["n_predicted_doublets"] = adata.obs.groupby("sample")["predicted_doublet"].sum().loc[df.index]
    df["prefilter_min_UMI"] = adata.obs.groupby("sample")["total_counts"].min().loc[df.index]
    df["prefilter_max_UMI"] = adata.obs.groupby("sample")["total_counts"].max().loc[df.index]
    df["prefilter_min_genes"] = adata.obs.groupby("sample")["n_genes_by_counts"].min().loc[df.index]
    df["prefilter_max_genes"] = adata.obs.groupby("sample")["n_genes_by_counts"].max().loc[df.index]

    df["prefilter_min_percent_mito"] = adata.obs.groupby("sample")["percent_mito"].min().loc[df.index]
    df["prefilter_max_percent_mito"] = adata.obs.groupby("sample")["percent_mito"].max().loc[df.index]
    max_genes = adata.obs[["sample","total_counts"]].groupby("sample").std()*2.5 + adata.obs[["sample","total_counts"]].groupby("sample").mean()
    cells_include = []
    for sample in adata.obs["sample"].unique():
        s_obs = adata.obs[adata.obs["sample"]==sample]
        print("lost "+ str(s_obs.shape[0]-len(list(s_obs.loc[adata.obs["total_counts"]<max_genes.loc[sample,"total_counts"]].index)))+" cells from "+ sample)
        cells_include += list(s_obs.loc[adata.obs["total_counts"]<max_genes.loc[sample,"total_counts"]].index)
    adata = adata[cells_include]
    adata = adata[adata.obs["predicted_doublet"]==0.0,:]
        
    adata.obs["value"] = 0
    sc.pp.filter_cells(adata, min_genes=400) 
    sc.pp.filter_cells(adata, min_counts=750)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['percent_mito'] < 0.2, :]
    df["post_filter_total_cells"] = adata.obs.groupby("sample").count()["donor"].loc[df.index]
    fig,ax = plt.subplots(2,3,figsize=(12,6))
    for i,stat in enumerate(['n_genes_by_counts', 'total_counts', "percent_mito"]):
        p = sc.pl.violin(adata, stat, groupby="donor",jitter=0.4, show=False,ax=ax[0,i])
        p=sc.pl.violin(adata, stat, groupby="sample",jitter=0.4, show=False,rotation=90,ax=ax[1,i])
    plt.tight_layout()
    fig.savefig("qc_plots_per_timepoint/"+timepoint+"_postfiltering.pdf",bbox_inches='tight')
    fig.savefig("qc_plots_per_timepoint/"+timepoint+"_postfiltering.png",bbox_inches='tight')
    df.to_csv("qc_plots_per_timepoint/"+timepoint+"_per_array_QC_stats.csv")
    return adata


def add_bm_doublet_scores(adata):
    all_doublets = pd.read_csv("../tables/all_scrublet_doublet_scores.csv",index_col=0)
    for col in all_doublets.columns:
            adata.obs[col] = all_doublets.loc[adata.obs_names,col]

def set_celltype_colors(adata, celltype_col="celltype"):
    celltype_color_file = "celltype_colors.csv"
    colors = pd.read_csv(celltype_color_file,  index_col=1)
    celltypes = adata.obs[celltype_col].unique()
    celltype_dict = {}
    for cell in colors.index:
        
        celltype_dict[cell] = colors.loc[cell,"color"]
        
    for c in celltypes:
        if c not in celltype_dict:
            print(c + " does not have a color assigned, choose one of:")
            print(celltype_dict.keys())
    set_colors_from_dict(adata, celltype_dict, celltype_col)


def add_donor_colors(adata, group_name="donor "):
    donor_color_df = pd.read_csv("donor_colors.csv",index_col=0)
    donor_color_dict = {i:donor_color_df.loc[i,"color"] for i in donor_color_df.index}
    set_colors_from_dict(adata, donor_color_dict, group_name)

def set_colors_from_dict(adata, colordict, group_name):
    adata.obs[group_name] = adata.obs[group_name].astype("category") 
    group_values_ordered = list(adata.obs[group_name].cat.categories)
    #group_values_ordered.sort()
    adata.uns[group_name+"_colors"] = [colordict[i] for i in group_values_ordered]

