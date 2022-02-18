import sklearn.cluster
from scipy.stats import zscore
from matplotlib.patches import Patch
import gseapy as gp
import numpy as np
import pandas as pd
import sys
import scanpy as sc
sys.path.insert(0,"/data/cb/nyquist/breast_milk/breastMilk/")
import heatmap_helper_functions as hh


def get_genelist_references(reference_file_path = "../../Data/",gene_sets=["GO_Biological_Process_2021"]):
    
    genelist_references = {}
    for s in gene_sets:
        genelist_references[s] = {}
        genelist_reference_file = open(reference_file_path+s+".txt")
    
        for l in genelist_reference_file:
            m = l.split("\t")
            genelist_references[s][m[0]] = m[1:]
    return genelist_references


def make_ordered_exp(epi_celltype_exp,celltypes,metadata,adata,celltype_col="celltype",lognorm=True,filter_expression=True):
    if type(celltypes)!=list:
        celltypes = [celltypes]
    exp = epi_celltype_exp[epi_celltype_exp[celltype_col].isin(celltypes)]
    exp.index=exp["sample"]
    exp = exp.iloc[:,2:]
    #exp =exp.dropna()
    # map expression to time post partum metadata
    exp["time_post_partum_days"] = exp.index.map(metadata["time_post_partum_days"])
    exp = exp.loc[exp["time_post_partum_days"]<400]
    exp = exp.iloc[:,:-1]
    exp=exp.loc[adata.obs[adata.obs["Epithelial Subclusters"].isin(celltypes)].groupby(["sample"]).count().loc[exp.index,"phase"] > 10]
    
    # remove genes not expressed
    exp=exp.loc[:,exp.sum(axis=0)>0]
    if lognorm:
        #sample normalize
        exp_norm = exp.div(exp.sum(axis=1),axis=0)*1000
        # log
        #exp_log=np.log(exp+1)
        exp_lognorm = np.log(exp_norm+1)
        #order exp by time post partum
    else:
        exp_lognorm = exp
    ordered_exp = exp_lognorm.iloc[exp_lognorm.index.map(metadata["time_post_partum_days"]).argsort()]
    return exp,ordered_exp

def heatmap_and_clusters_by_time(epi_celltype_exp, des_res, celltype,metadata,adata,minlfc=0.005,minmean=20,vmax=3,vmin=-2, min_pts = .1):
    directory = "time_series_heatmaps/"
    exp,ordered_exp = make_ordered_exp(epi_celltype_exp, celltype,metadata,adata)
    
    if "rank_genes_groups" not in adata_all_epi.uns or adata_all_epi.uns["rank_genes_groups"]["params"]["groupby"] != "Epithelial Subclusters" or "pts" not in adata_all_epi.uns["rank_genes_groups"]:
            sc.tl.rank_genes_groups(adata_all_epi, groupby="Epithelial Subclusters", pts=True)
    des_res_reduced = des_res.loc[des_res["padj"]<.05]
    des_res_reduced = des_res_reduced.loc[des_res_reduced["log2FoldChange"].abs()>minlfc]
    des_res_reduced = des_res_reduced.loc[des_res_reduced["baseMean"].abs()>minmean]
    #g = [i.replace(".","_") for i in des_res_reduced.index]
    overlap_genes = list(set(des_res_reduced.index).intersection(set(adata_all_epi.uns["rank_genes_groups"]["pts"].index)))
    #des_res_reduced.index =  [i.replace(".","_") for i in des_res_reduced.index]
    des_res_reduced = des_res_reduced.loc[overlap_genes]
    des_res_reduced["pts"] = adata_all_epi.uns["rank_genes_groups"]["pts"].loc[des_res_reduced.index,celltype]
    des_res_reduced = des_res_reduced.loc[des_res_reduced["pts"]>min_pts]
    genes=[i for i in des_res_reduced.sort_values('log2FoldChange').index if i in ordered_exp.columns]
    #zscore each column
    z=ordered_exp.apply(zscore)
    n_clusters = 5
    labels = sklearn.cluster.KMeans(n_clusters=n_clusters).fit_predict(z.T.loc[genes])
    new_gene_order=reorder_from_labels(labels,genes)
    lut=dict(zip(list(set(labels)),("r","g","y","m","k")))
    row_colors=[lut[i] for i in labels]
    row_color_order = reorder_from_labels(labels,row_colors)
    exp.iloc[exp.index.map(metadata["time_post_partum_days"]).argsort()][new_gene_order].to_csv(directory+celltype+"_reduced_pseudobulk_expression_for_heatmap_raw.csv")
    col_colors=ordered_exp.index.map(metadata["milk_stage"]).map(hh.milk_stage_colors)
    ordered_exp[new_gene_order].to_csv(directory+celltype+"_reduced_pseudobulk_expression_for_heatmap_lognormed.csv")
    pd.DataFrame(labels, index=genes).to_csv(directory+celltype+"_time_dep_gene_cluster_labels.csv")
    g=sns.clustermap(ordered_exp.T.loc[new_gene_order],row_cluster=False,col_cluster=False,row_colors=row_color_order,col_colors=col_colors,z_score=0,vmax=vmax,vmin=vmin)
    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(handles, lut, title='Gene Cluster',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    plt.savefig(directory+celltype+"_time_dependent_genes_heatmap.pdf",bbox_inches="tight")
    plt.savefig(directory+celltype+"_time_dependent_genes_heatmap.png",bbox_inches="tight")    
    return genes,labels


def gsea_prerank_heatmaps(epi_celltype_exp, des_res, celltype,metadata,adata,gene_sets="GO_Biological_Process_2021"):
    
    outdir='time_series_heatmaps/prerank/'+celltype.replace("/","_").replace(" ","_")+'/prerank_report_hallmark'
    des_res_reduced = des_res.loc[des_res["padj"]<.05]
    genes_gsea = des_res_reduced.sort_values('log2FoldChange').index
    pre_res = gp.prerank(rnk=des_res_reduced.loc[genes_gsea,"log2FoldChange"], gene_sets=gene_sets,
                     processes=4,
                     permutation_num=100, # reduce number to speed up testing
                     outdir=outdir, format='png', seed=6)
    _,exp = make_ordered_exp(epi_celltype_exp,celltype,metadata,adata)
    z=exp.apply(zscore)
    for pathway in pre_res.res2d[pre_res.res2d["pval"] < .05].index:
        g=sns.clustermap(z.T.loc[pre_res.res2d.loc[pathway,"genes"].split(";")],row_cluster=False,col_cluster=False,vmax=3)
        plt.title(pathway)
        plt.savefig(outdir+"/"+pathway.replace("-","").replace(" ","_").replace("/","_")+"_heatmap.png",bbox_inches="tight")
    pre_res.res2d.to_csv(outdir+"/"+gene_sets+"_prerank_results.csv")
    return pre_res
    
    
def dotplot_of_pre_res(pre_res, celltype):
    ordered = pre_res.res2d.iloc[pre_res.res2d["nes"].argsort()]
    ordered = ordered[(ordered["fdr"]<.25) & (ordered["pval"]<.05)]
    plt.figure(figsize=(5,6))
    sns.scatterplot(y=ordered.index,x=ordered["nes"],size=ordered["matched_size"],hue=ordered["fdr"])
    plt.xlabel("Normalized Enrichment Score")
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.title(celltype+" Time Associated Genes")

def get_mean_scores_for_heatmap(adata,celltypes, GO_list,metadata):
    if type(celltypes)==list:
        adata.obs["tmp_celltype"] = ""
        adata.obs.loc[adata.obs["Epithelial Subclusters"].isin(celltypes),"tmp_celltype"] = "_".join(celltypes)
        celltypes.append("_".join(celltypes))
        mean_scores_per_celltype=adata.obs.groupby(["tmp_celltype","sample"])[GO_list].mean().reset_index()
        _,sec_ordered = make_ordered_exp(mean_scores_per_celltype,celltypes,metadata,adata,celltype_col="tmp_celltype",lognorm=False)

    else:
        mean_scores_per_celltype=adata.obs.groupby(["Epithelial Subclusters","sample"])[GO_list].mean().reset_index()
    
        _,sec_ordered = make_ordered_exp(mean_scores_per_celltype,celltypes,metadata,adata,celltype_col="Epithelial Subclusters",lognorm=False)
    
    combined_scores = sec_ordered.T
    combined_scores.columns = combined_scores.columns.astype(str)
    combined_scores=combined_scores[sec_ordered.index]
    return combined_scores

    
def collapse_GO_hits(GO_hits, enr,overlap_threshold = 0.6):
    if type(enr) != type(pd.DataFrame()):
        enr_df = enr.res2d
    else:
        enr_df = enr
    overlap_info = {}
    ordered_go_hits= enr_df.loc[enr_df["Term"].isin(GO_hits)].sort_values("n_genes")["Term"].values
    print(len(ordered_go_hits))
    for g in ordered_go_hits:
        found_overlap = False
        genes = enr_df.loc[enr_df["Term"]==g,"Genes"].values[0].split(";")
        #print(genes)
        genelist_len = len(genes)
        max_overlap = 0
        max_overlap_key = ""
        for s in overlap_info:
            
            overlap = len(set(genes).intersection(overlap_info[s]["Genes"]))
            if overlap>max_overlap and (1.0*overlap)/genelist_len > overlap_threshold and not found_overlap:
                found_overlap = True
                max_overlap=overlap
                max_overlap_key = s
        if found_overlap:
            overlap_info[max_overlap_key]["Genes"] = overlap_info[max_overlap_key]["Genes"].union(genes)
            overlap_info[max_overlap_key]["listnames"].append(g)
            overlap_info[max_overlap_key]["combined_scores"].append(enr_df.loc[enr_df["Term"]==g,"Combined Score"].values[0])
        if not found_overlap:
            overlap_info[g] = {}
            overlap_info[g]["Genes"] = set(genes)
            overlap_info[g]["listnames"] = [g]
            overlap_info[g]["combined_scores"] = [enr_df.loc[enr_df["Term"]==g,"Combined Score"].values[0]]
    collapsed_list = []
    for o in overlap_info:
        top = overlap_info[o]['listnames'][np.argmax(overlap_info[o]["combined_scores"])]
        overlap_info[o]["collapsed_listname"] = top
        collapsed_list.append(top)
    return overlap_info, collapsed_list

def enr_and_score_genes(adata, genelist_use,genelist_references,plots_title,gene_set="GO_Biological_Process_2021",overlap_threshold=0.6 ):
    enr=gp.enrichr(gene_list=genelist_use,gene_sets=gene_set, organism="human")
    gp.dotplot(enr.results, title=plots_title, figsize=(4,8),  top_term=20,cmap="Greens")
    enr.res2d["n_genes"] = [len(i.split(";")) for i in enr.res2d["Genes"]]
    #print(enr.res2d)
    GO_hits = list(enr.res2d.loc[(enr.res2d["n_genes"]>=4)&(enr.res2d["Adjusted P-value"]<=.05),"Term"])
    #print(GO_hits)
      
    #overlap_info, collapsed_list = collapse_GO_hits(GO_hits,enr,overlap_threshold=overlap_threshold)
    for list_name in GO_hits:
        genelist = genelist_references[gene_set][list_name]
        if list_name not in adata.obs.columns:
            sc.tl.score_genes(adata, genelist, score_name=list_name,use_raw=True)
    return enr,GO_hits

def remove_uncorrelated_scores(adata, GO_hits,celltype, corr_direction,metadata):
    mean_scores_per_celltype=adata.obs.groupby(["Epithelial Subclusters","sample"])[GO_hits].mean().reset_index()
    _,sec_ordered = make_ordered_exp(mean_scores_per_celltype,celltype,metadata,adata,celltype_col="Epithelial Subclusters",lognorm=False)

    combined_scores = sec_ordered.T
    combined_scores.columns = combined_scores.columns.astype(str)
    combined_scores=combined_scores[sec_ordered.index]
   

    if corr_direction == "up":
        new_increasing_list = sec_ordered.columns[np.corrcoef(sec_ordered.index.map(metadata["time_post_partum_days"]), sec_ordered.T)[0][1:] > .1]
    else:
        new_increasing_list = sec_ordered.columns[np.corrcoef(sec_ordered.index.map(metadata["time_post_partum_days"]), sec_ordered.T)[0][1:] < -.1]
    return new_increasing_list

def collapsed_enrichr_analysis(adata, genelist_use, celltype, plots_title,genelist_references,metadata,gene_set="GO_Biological_Process_2021",overlap_threshold=0.6,corr_direction="",go_res_dir = ""):
    
    enr,GO_hits = enr_and_score_genes(adata, genelist_use,genelist_references,plots_title,gene_set=gene_set,overlap_threshold=overlap_threshold )
    go_used = enr.res2d[enr.res2d["Term"].isin(GO_hits)]
    go_used.index = go_used["Term"]
    if corr_direction != "":
        GO_hits = remove_uncorrelated_scores(adata, GO_hits,celltype, corr_direction,metadata)
                                
    #mean_scores_per_celltype=adata.obs.groupby(["Epithelial Subclusters","sample"])[GO_hits].mean().reset_index()
    overlap_info, collapsed_list = collapse_GO_hits(GO_hits,enr)
    
    
    go_used["collapsed_into_geneset"] = ""
    for o in overlap_info:
        top = overlap_info[o]['listnames'][np.argmax(overlap_info[o]["combined_scores"])]
        for g in overlap_info[o]['listnames']:
            go_used.loc[g,"collapsed_into_geneset"] = top

    
    go_used.to_csv(go_res_dir +"/enrichr_collapsed_"+plots_title+".csv")


    return enr,collapsed_list,GO_hits

def GO_term_reduced_heatmap(adata,celltype,collapsed_list,epi_sub_colors,metadata):
    combined_scores = get_mean_scores_for_heatmap(adata,celltype, collapsed_list)
    
    col_colors=combined_scores.columns.map(metadata["milk_stage"]).map(hh.milk_stage_colors)
    row_colors = [epi_sub_colors[celltype]]*len(collapsed_list)
    g=sns.clustermap( combined_scores,row_cluster=False,col_cluster=False,col_colors=col_colors,row_colors=row_colors,z_score=0,figsize=(10,10),yticklabels=True)
    #g.savefig(go_res_dir+"/"+plots_title+"_enrichr_heatmap.pdf",bbox_inches="tight")


def make_GO_term_metadata(adata_all_epi, group_column = "Epithelial Subclusters"):
    go_kwds = [i for i in adata_all_epi.obs.columns if "GO" in i]
    celltype_means = adata_all_epi.obs[[group_column,]+list(go_kwds)].groupby([group_column]).mean()
    # take the ones that vary the most - will probably still need to weed out from these though
    stds = celltype_means.std()
    max_score_ids = celltype_means.idxmax()

      # maybe the tie break for parent/child conflicts can be keep both if they are up in different clusters but drop wht one with the lower std if they are up in the same cluster
    # what we want to know is, if one is more specific and the other is just a more general version of that or if the other is capturing somethign different

    min_score_ids = celltype_means.idxmin()
    fold_changes = [(celltype_means.loc[max_score_ids[s],s]-celltype_means.loc[min_score_ids[s],s])/np.abs(celltype_means.loc[min_score_ids[s],s]) for s in celltype_means.columns]
    geneset_metadata = pd.DataFrame(index=celltype_means.columns)
    geneset_metadata["max score celltype"] = max_score_ids
    geneset_metadata["std"] = stds
    geneset_metadata["min score celltype"] = min_score_ids
    geneset_metadata["fold_change"] = fold_changes
    return geneset_metadata
    
    
def build_parents_and_children_dict(geneset_metadata, genelist_references,overlap_threshold=0.4,gene_set="GO_Biological_Process_2021"):
    '''
    make a dictionary of geneset name : {"parents":set of parents with overlaps in the reference > 0.4, "children":set of children with overlaps in the reference > 0.4} based off the genes in the original GO references
    
    geneset metadata is a dataframe whose index includes the list of all GO terms which have been gene set scored in the adata you are using
    
    genelist_references is a dict of gene set to geneset names to list of genes
    '''
    geneset_build = list(geneset_metadata.index)
    parents_and_children = {i:{"parents":set(),"children":set()} for i in geneset_build}

    
    for l,g1 in enumerate(geneset_build):

        for g2 in geneset_build[l+1:]:
            genelist1 = set(genelist_references[gene_set][g1])
            genelist2 = set(genelist_references[gene_set][g2])
            ref_overlap = len(genelist1.intersection(genelist2))-2
            if ref_overlap>1 and (1.0*ref_overlap)/len(genelist1) > overlap_threshold:
                parents_and_children[g1]["parents"].add(g2)
                parents_and_children[g2]["children"].add(g1)
            elif ref_overlap>1 and (1.0*ref_overlap)/len(genelist2) > overlap_threshold:
                parents_and_children[g2]["parents"].add(g1)
                parents_and_children[g1]["children"].add(g2)
    return parents_and_children


def find_children(parents_and_children, c):
    all_children = set()
    children_to_check = set([c,])
    while len(children_to_check) >0:
        c = children_to_check.pop()
        all_children.add(c)
        new_children = parents_and_children[c]["children"] - all_children
        all_children = all_children.union(new_children)
        children_to_check = children_to_check.union(new_children)
    #for cl in parents_and_children[c]["children"]:
    #    all_children.add(cl)
    #    all_children = all_children.union(find_children(parents_and_children, cl))
    return all_children



def collapse_by_parents(parents_and_children,geneset_metadata,pathways_use):
    reduced_by_parents = set()
    #print("hi")
    for gl,a in parents_and_children.items():
        #print(gl)
        if len(a["parents"]) == 0 and len(a["children"])>= 1:
            #print(gl)
            if len(a["children"])<5:
                check_children = a["children"]
            
            else:
                check_children = [gl]
            for c in check_children:
                all_children = find_children(parents_and_children, c)

                all_poss_paths = set([c,]).union(all_children)

                all_poss_paths = list(all_poss_paths.intersection(pathways_use))
                all_added = list(geneset_metadata.loc[all_poss_paths].groupby("max score celltype")["std"].idxmax().values)
                reduced_by_parents = reduced_by_parents.union(set(all_added)) 
        elif len(a["parents"]) == 0 and gl in pathways_use and len(a["children"])<5:
            reduced_by_parents.add(gl)
    return reduced_by_parents