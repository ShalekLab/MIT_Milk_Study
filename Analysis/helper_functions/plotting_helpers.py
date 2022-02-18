import seaborn as sns
import pandas as pd
from collections import defaultdict
from matplotlib import colors
import matplotlib.pylab as plt
from scipy.stats import zscore
import scanpy as sc
import matplotlib.pyplot as plt

#______ UTILS________
def reorder_from_labels(labels, index):
    # order based on labels:
    clust_to_sample = defaultdict(list)
    for i,s in enumerate(index):
        clust_to_sample[labels[i]] += [s]
    
    new_order = []
    for clust,samps in clust_to_sample.items():
        new_order += samps
    return new_order
def reorder_from_multiple_labels(labels_df,index,labels_order):
    clust_to_sample = defaultdict(list)
    cur_label = labels_order[0]
    labels_order = labels_order[1:]
    
    for i,s in enumerate(index):
        
        clust_to_sample[labels_df.loc[s,cur_label]] += [s]
    
    new_order = []
    # impose an order on the samples
    clusts = sorted(clust_to_sample.keys())
    for clust in clusts:
        samps = clust_to_sample[clust]
        if len(labels_order) == 0: # base case, just reordering on one label
            new_order += samps
        else:
            new_order += reorder_from_multiple_labels(labels_df, samps,labels_order)
    return new_order

def order_labels(df, colname, correct_order):
    order = []
    for v in correct_order:
        if v in df[colname].values:
            order+=[v]
    print(order)
    return df.set_index(colname).loc[order]

def make_proportions_df(adata, x_value, color_value, hue):
    
    
    tmp = adata.obs.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)

    m=tmp.divide(tmp.sum(axis=1), axis=0)
    props = []

    i=0

    for sample in m.index:

        for celltype in m.columns:
            vals = [sample,m.loc[sample,celltype],celltype,adata.obs.loc[adata.obs[x_value]==sample,hue].unique()[0]]
            props.append(vals)
            i+=1
    props_df = pd.DataFrame(props,columns=[x_value,x_value+"_proportion",color_value,hue])
    props_df[hue]=props_df[hue].astype("category")
    return props_df

def qcplots(gran_adata, groupby="leiden", gs4=None,fig=None, donor_colname = "M.Number",sample_colname="sample",include_stackedbars=True):
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker
    if gs4 is None:
        if include_stackedbars:
            fig=plt.figure(figsize=(7,15))
            gs4 = gridspec.GridSpec(6,1)
        else:
            fig=plt.figure(figsize=(7,11))
            gs4 = gridspec.GridSpec(4,1)
    ax_tc = fig.add_subplot(gs4[0, 0])
    #else:
        #gs4 = ax.get_subplotspec()
        #ax_tc=ax
    sc.pl.violin(gran_adata, "total_counts",groupby=groupby,rotation=90,ax=ax_tc,show=False,stripplot=False)
    ax_tc.set_xlabel("")
    ax_tc.set_xticklabels([])
    ax_tc.set_xticks([])
    ax_tc.set_ylabel("n_UMI")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1))
    ax_tc.yaxis.set_major_formatter(formatter)
    ax_mito = fig.add_subplot(gs4[1, 0])
    sc.pl.violin(gran_adata, "percent_mito",groupby=groupby,rotation=90,ax=ax_mito,show=False, stripplot=False)
    ax_mito.set_xlabel("")
    ax_mito.set_xticklabels([])
    ax_mito.set_xticks([])
    ax_mito.set_ylabel("%mito")
    ax_genes = fig.add_subplot(gs4[2, 0])
    sc.pl.violin(gran_adata, "n_genes_by_counts",groupby=groupby,rotation=90,ax=ax_genes,show=False, stripplot=False)
    ax_genes.set_xlabel("")
    ax_genes.set_xticklabels([])
    ax_genes.set_xticks([])
    ax_genes.set_ylabel("n_genes")
    formatter_g = ticker.ScalarFormatter(useMathText=True)
    formatter_g.set_scientific(True) 
    formatter_g.set_powerlimits((-1,1))
    ax_genes.yaxis.set_major_formatter(formatter_g)
    ax_doublet = fig.add_subplot(gs4[3, 0])
    sc.pl.violin(gran_adata, "doublet_scores",groupby=groupby,rotation=90,ax=ax_doublet,show=False, stripplot=False)
    ax_doublet.set_ylabel("doublet\nscores")
    if include_stackedbars:
        ax_doublet.set_xlabel("")
        ax_doublet.set_xticklabels([])
        ax_doublet.set_xticks([])
        ax_sample = fig.add_subplot(gs4[4, 0])
        normalized_stacked_bar_plot(gran_adata, groupby,sample_colname,ax=ax_sample,legend=False)
        ax_sample.set_xlabel("")
        ax_sample.set_xticklabels([])
        ax_sample.set_xticks([])
        ax_sample.set_ylabel("pct cells")
        ax_monkey = fig.add_subplot(gs4[5, 0])
        hh.normalized_stacked_bar_plot(gran_adata, groupby,donor_colname,ax=ax_monkey,legend=False)
        ax_monkey.set_ylabel("pct cells")
def draw_clustergraph(adata_all_epi, reslist,y_vals=[],orders=[]):
    '''
    plots cluster membership graph from anndata object
    Prior to running this, clustering must be run at each resolution you are interested in. 
    You must also have colors saved in adata.uns[<cluster key>_colors] for each resolution (you can force scanpy to do this by plotting a umap with color=<cluster key>)
    Also the each cluster resolution column in obs must be of type categorical (also happens when you plot it)
    
    Requires networkx to run
    
    Inputs:
    adata - an adata object that meets the requirements listed above
    reslist - a list, in order, of the cluster keys for the clustering resolutions of interest (example:["leiden_res0.5","leiden_res0.6","leiden_res0.7"])
        this list can also include other cell clusterings, for example labels that came from mergind leiden clusters, etc
    y_vals - a dictionary mapping each value in reslist to its height in the figure. Defaults to plotting them 2 points apart
    orders - a dictionary mapping each value in reslist to a list specifying the order along the x axis of the categories in that resolution, defaults to a random or alphabetical order
    '''
    import networkx as nx
    # first set up spacings and orderings
    if len(y_vals) == 0: # set the y values of each resolution 
        y_vals = dict(zip(reslist,[i*2 for i,_ in enumerate(reslist)]))
    if len(orders) ==0: # chooses an order of nodes for each resolution if not provided
        orders = {}
        for r in reslist:
            orders[r] = [str(i) for i in adata_all_epi.obs[r].unique()]
            
    # space nodes at each resolution along the full x axis
    
    # get max number of clusters
    lens = [len(l) for o,l in orders.items()] 
    maxwidth = max(lens)
    
    x_vals = {}
    # calculate the x value for each node so they are spaced along full x axis
    for o,l in orders.items(): 
        w = len(l)
        spacing = maxwidth/(w*1.0)
        x_vals[o] = {v:(i*spacing) for i,v in enumerate(l)}
        
    # calculate edges for each consecutive pair of resolutions in reslist
    respairs = [(reslist[i],reslist[i+1]) for i in range(len(reslist)-1)]
    edges = []
    for res1,res2 in respairs:
        
        # edge weights come from the proportion of cells from each cluster in the top resolution
        # that are assigned to each cluster in the next resolution
        # if no cells map between a pair of clusters, there is no edge
        layer_1_counts = adata_all_epi.obs.groupby([res1,res2]).count()["sample"].unstack()
        layer_1_props = layer_1_counts.divide(layer_1_counts.sum(axis=1),axis=0)
        edgelist = layer_1_props.stack().reset_index()
        
        # so each node gets a unique name, name nodes <resolution>.<name>
        edgelist["top"]=[res1+"."+i for i in edgelist[res1]]
        edgelist["bottom"]=[res2+"."+i for i in edgelist[res2]]
        
        # edges are saved as (sender node, receiver node, weight)
        edges+=[(edgelist.loc[i,"top"],edgelist.loc[i,"bottom"],edgelist.loc[i,0]) for i in edgelist.index]
   
    # initialize graph and add edges
    G = nx.DiGraph()
    G.add_weighted_edges_from(edges) 
    
    # set graph asthetics for all nodes
    all_sizes = {}
    all_colors = {}
    all_labels = {}
    size_mult = 1000 # sizes need a multiplier to make the nodes visible
    for r in reslist:
        
        # node sizes are proportional to the proportion of cells in each cluster at that resolution
        pcts = adata_all_epi.obs.groupby(r).count()["sample"]/adata_all_epi.obs.groupby(r).count()["sample"].sum()
        d = pcts.to_dict()
        d_update = {r+"."+str(i):size_mult*d[i] for i in d.keys()}
        all_sizes.update(d_update) 
        
        # colors come from the adata colors
        colors = dict(zip(adata_all_epi.obs[r].cat.categories,adata_all_epi.uns[r+"_colors"]))
        col_update = {r+"."+str(i):colors[i] for i in colors.keys()}
        all_colors.update(col_update)
        
        # reset the labels so they match the values in obs instead of the <resolution>.<name> node names we used to make node names unique
        all_labels.update({r+"."+str(i):i for i in d.keys()})


    # set up position of each node depending on the values calculated above for orders and y_vals
    pos = {}
    for i in G.nodes():
        for r in reslist:
            if r in i:
                res=r

        pos[i] = [x_vals[res][i.split(res)[1][1:]],y_vals[res]]
    
    # plot the graph
    #
    #print(min([G[u][v]['weight'] for u,v in G.edges]))
    #print(max([G[u][v]['weight'] for u,v in G.edges]))
    nx.draw(G,pos=pos,width= [2*G[u][v]['weight'] for u,v in G.edges],edge_vmin=-1,edge_color= [G[u][v]['weight'] for u,v in G.edges],node_color=[all_colors[i] for i in G.nodes],node_size=[all_sizes[i] for i in G.nodes],labels=all_labels,edge_cmap=plt.cm.Greys)
    for r in reslist:
        plt.text(-5,y_vals[r],r)
    return G, all_colors,all_sizes,all_labels

def draw_clustergraph_from_df(adata_df, reslist,y_vals=[],orders=[]):
    '''
    
    LOTS OF REPEATED CODE FROM draw_clustergraph
    plots cluster membership graph from dataframe object 
    Prior to running this, clustering must be run at each resolution you are interested in. 
    You must also have colors saved in adata.uns[<cluster key>_colors] for each resolution (you can force scanpy to do this by plotting a umap with color=<cluster key>)
    Also the each cluster resolution column in obs must be of type categorical (also happens when you plot it)
    
    Requires networkx to run
    
    Inputs:
    adata_df - equivalant to obs value of df that meets the requirements listed above
    reslist - a list, in order, of the cluster keys for the clustering resolutions of interest (example:["leiden_res0.5","leiden_res0.6","leiden_res0.7"])
        this list can also include other cell clusterings, for example labels that came from mergind leiden clusters, etc
    y_vals - a dictionary mapping each value in reslist to its height in the figure. Defaults to plotting them 2 points apart
    orders - a dictionary mapping each value in reslist to a list specifying the order along the x axis of the categories in that resolution, defaults to a random or alphabetical order
    '''
    import networkx as nx
    # first set up spacings and orderings
    if len(y_vals) == 0: # set the y values of each resolution 
        y_vals = dict(zip(reslist,[i*2 for i,_ in enumerate(reslist)]))
    if len(orders) ==0: # chooses an order of nodes for each resolution if not provided
        orders = {}
        for r in reslist:
            orders[r] = [str(i) for i in adata_df[r].unique()]
            
    # space nodes at each resolution along the full x axis
    
    # get max number of clusters
    lens = [len(l) for o,l in orders.items()] 
    maxwidth = max(lens)
    
    x_vals = {}
    # calculate the x value for each node so they are spaced along full x axis
    for o,l in orders.items(): 
        w = len(l)
        spacing = maxwidth/(w*1.0)
        x_vals[o] = {v:(i*spacing) for i,v in enumerate(l)}
        
    # calculate edges for each consecutive pair of resolutions in reslist
    respairs = [(reslist[i],reslist[i+1]) for i in range(len(reslist)-1)]
    edges = []
    for res1,res2 in respairs:
        
        # edge weights come from the proportion of cells from each cluster in the top resolution
        # that are assigned to each cluster in the next resolution
        # if no cells map between a pair of clusters, there is no edge
        layer_1_counts = adata_df.groupby([res1,res2]).count()['time_post_partum_weeks'].unstack()
        layer_1_props = layer_1_counts.divide(layer_1_counts.sum(axis=1),axis=0)
        edgelist = layer_1_props.stack().reset_index()
        
        # so each node gets a unique name, name nodes <resolution>.<name>
        edgelist["top"]=[res1+"."+i for i in edgelist[res1]]
        edgelist["bottom"]=[res2+"."+i for i in edgelist[res2]]
        
        # edges are saved as (sender node, receiver node, weight)
        edges+=[(edgelist.loc[i,"top"],edgelist.loc[i,"bottom"],edgelist.loc[i,0]) for i in edgelist.index]
   
    # initialize graph and add edges
    G = nx.DiGraph()
    G.add_weighted_edges_from(edges) 
    
    # set graph asthetics for all nodes
    all_sizes = {}
    all_colors = {}
    all_labels = {}
    size_mult = 1000 # sizes need a multiplier to make the nodes visible
    for r in reslist:
        
        # node sizes are proportional to the proportion of cells in each cluster at that resolution
        pcts = adata_df.groupby(r).count()['time_post_partum_weeks']/adata_df.groupby(r).count()['time_post_partum_weeks'].sum()
        d = pcts.to_dict()
        d_update = {r+"."+str(i):size_mult*d[i] for i in d.keys()}
        all_sizes.update(d_update) 
        
        # colors come from the adata colors
        colors = dict(zip(adata_df[r].cat.categories,["#EDF8E9","#BAE4B3","#74C476","#31A354","#006D2C","#02421C","#FFFFFF"]))
        col_update = {r+"."+str(i):colors[i] for i in colors.keys()}
        all_colors.update(col_update)
        
        # reset the labels so they match the values in obs instead of the <resolution>.<name> node names we used to make node names unique
        all_labels.update({r+"."+str(i):i for i in d.keys()})


    # set up position of each node depending on the values calculated above for orders and y_vals
    pos = {}
    for i in G.nodes():
        for r in reslist:
            if r in i:
                res=r

        pos[i] = [x_vals[res][i.split(res)[1][1:]],y_vals[res]]
    
    # plot the graph
    #
    #print(min([G[u][v]['weight'] for u,v in G.edges]))
    #print(max([G[u][v]['weight'] for u,v in G.edges]))
    nx.draw(G,pos=pos,width= [2*G[u][v]['weight'] for u,v in G.edges],edge_vmin=-1,edge_color= [G[u][v]['weight'] for u,v in G.edges],node_color=[all_colors[i] for i in G.nodes],node_size=[all_sizes[i] for i in G.nodes],labels=all_labels,edge_cmap=plt.cm.Greys)
    for r in reslist:
        plt.text(-1,y_vals[r],r)
    return G, all_colors,all_sizes,all_labels


#________ STACKED BAR PLOTS____________
def normalized_stacked_bar_plot(adata, x_value, color_value, palette=None, legend=True,ax=None,x_order=None,color_order=None, log=False):
    if color_value+"_colors" in adata.uns:
        color_dict = dict(zip(adata.obs[color_value].cat.categories,adata.uns[color_value+"_colors"]))
        if color_order is not None:
            palette = colors.ListedColormap([color_dict[c] for c in color_order])
        else:
            palette = colors.ListedColormap(adata.uns[color_value+"_colors"])
    if x_value=="milk stage":
        df = order_labels(adata.obs, x_value, ["early", "transitional","transitional ","mature","late"])
    else:
        df = adata.obs
    tmp = df.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)
    if x_order is not None:
        tmp = tmp.loc[x_order]
    normed_tmp = tmp.divide(tmp.sum(axis=1), axis=0)
    if log == True:
        normed_tmp = -1.0*np.log10(normed_tmp)
        normed_tmp = normed_tmp.replace(np.inf, 0)
        print(normed_tmp)
    if color_order is not None:
        #normed_tmp.columns = pd.CategoricalIndex(color_order, ordered=True, categories=color_order) 
        #normed_tmp = normed_tmp.sort_index(axis=1)
        normed_tmp = normed_tmp[color_order]
    if ax is not None:
        normed_tmp.plot(kind='bar',stacked=True,  colormap=palette, ax=ax)
    else:
        ax =normed_tmp.plot(kind='bar',stacked=True, figsize=(5,7), colormap=palette)
    plt.ylabel("proportion of cells")
    if legend:
        plt.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8))
    else:
        ax.legend().set_visible(False)
    if log == False:
        ax.set_ylim(0,1.1)


def stacked_bar_plot(adata,x_value,color_value, palette=None, legend=True, ax=None):
    if color_value+"_colors" in adata.uns:
        palette = colors.ListedColormap(adata.uns[color_value+"_colors"])
    #TODO: allow for coloring based on colors stored in adata
    #if x_value=="milk stage":
    #    df = order_labels(adata.obs, x_value, ["early", "transitional ","transitional","mature","late"])
    #else:
    df = adata.obs
    if ax is not None:

        df.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0).plot(kind='bar',stacked=True,colormap=palette, ax=ax)
    else:
        ax=df.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0).plot(kind='bar',stacked=True,colormap=palette)
    plt.ylabel("n_cells")
    if legend:
        plt.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8))
    else:
        ax.legend().set_visible(False)


def multiple_stacked_bar_plots(adata, plot_sep,x_value,color_value,normed=True,palette=None):
    nplots = len(adata.obs[plot_sep].unique())
    plot_no = 1
    legend = False
    for p in adata.obs[plot_sep].unique():
        plt.subplot(1,nplots,plot_no)
        plot_no += 1
        if plot_no==nplots:
            legend = True

        a = adata[adata.obs[plot_sep]==p]
        if normed:
            normalized_stacked_bar_plot(a, x_value, color_value, palette=palette, legend=legend)
        else:
            stacked_bar_plot(a, x_value,color_value)



def normalized_stacked_bar_plot_from_df(df, x_value, color_value, palette=None,color_uns=None,ax=None):
    if color_uns is not None:
        palette = colors.ListedColormap(color_uns)
    #if x_value=="milk stage":
    #    df = order_labels(adata.obs, x_value, ["early", "transitional","transitional ","mature","late"])
    
    tmp = df.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)
    if ax is not None:
        tmp.divide(tmp.sum(axis=1), axis=0).plot(kind='bar',stacked=True, ax=ax, colormap=palette)
    else:
        ax =tmp.divide(tmp.sum(axis=1), axis=0).plot(kind='bar',stacked=True, figsize=(5,7), colormap=palette)
    ax.set_ylabel("proportion of cells")
    plt.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8))
    ax.set_ylim(0,1.1)
def stacked_bar_with_errorbars(adata, x_value, color_value, hue,palette=None):
    props_df = make_proportions_df(adata, x_value, color_value, hue)
    mean_props = props_df.groupby([hue,color_value]).agg([np.mean,np.std])#.unstack(color_value)["mean"]#.plot(kind='bar',y="mean",yerr="std",stacked=True)

    means = mean_props[('sample_proportion', 'mean')].unstack(color_value)

    std = mean_props[('sample_proportion', 'std')].unstack(color_value)
    ind = np.arange(means.shape[0])
    bottom=np.zeros(means.shape[0])
    plt.figure(figsize=(5,10))
    if type(palette) == NoneType:
        if hue+"_colors" in adata.uns:
            palette = dict(zip(adata.obs[hue].cat.categories,adata.uns[hue+"_colors"]))
        else:
            palette=None
    for c in means.columns:
        p1 = plt.bar(ind, means[c], .5, yerr=std[c],bottom=bottom, color=palette[c])
        bottom = bottom + means[c]



#____grouped by condition plots____


def grouped_dotplot(x_condition,y_condition,genes,adata,ordered_y_condition=[],within_gene_sep = .6,
                   between_gene_sep = .8):
    tmp_obs = adata.obs

    for G in genes:

        tmp_obs[G] = adata.raw[:,G].X
        tmp_obs[G+" on"] = adata.raw[:,G].X > 0

    means = tmp_obs.groupby([x_condition,y_condition])[genes].mean().stack().reset_index()
    pcts = (tmp_obs.groupby([x_condition,y_condition])[[g+" on" for g in genes]].sum()/tmp_obs.groupby([x_condition,y_condition])[[g+" on" for g in genes]].count())#.stack()
    pcts.columns = genes
    means["pcts"]=pcts.stack().reset_index()[0]
    means.columns = [x_condition,y_condition, "gene","mean","pcts"]
    #zscore the means
    means["zmeans"] =means.groupby("gene").transform(lambda x: zscore(x,ddof=1))["mean"]
    means["x_label_name"]= [means.loc[i,"gene"]+means.loc[i,x_condition] for i in means.index]
    x_coords = []#list(range(len(genes)*len(means[x_condition].unique())))
    
    
    ordered_x_condition = means[x_condition].unique()
    x_labelnames = []
    x_coord_value = 0
    linepositions = []
    gene_label_locs = []
    for g in genes:
        x_labelnames += [g+l for l in ordered_x_condition]
        x_coords += [x_coord_value + between_gene_sep,]+ [x_coord_value+between_gene_sep + (l+1)*within_gene_sep for l in range(len(ordered_x_condition)-1)]
        added_space = between_gene_sep+(within_gene_sep*(len(ordered_x_condition)-1))
        gene_label_locs+=[x_coord_value + between_gene_sep+(within_gene_sep*((len(ordered_x_condition)-1)/2.0))]
        x_coord_value+= added_space
        linepositions += [x_coord_value + (between_gene_sep/2.0)]

    x_coord_map = dict(zip(x_labelnames,x_coords))
    means["xcoord"]= means["x_label_name"].map(x_coord_map)
    if len(ordered_y_condition) == 0:
        ordered_y_condition = means[y_condition].unique()
    y_coords = range(len(ordered_y_condition))
    y_coord_map =dict(zip(ordered_y_condition, y_coords))
    means["ycoord"] = means[y_condition].map(y_coord_map)
    figheight=len(y_coords)*.38
    figwidth=len(x_coords)*.4
    plt.figure(figsize=(figwidth,figheight))

    ax=sns.scatterplot(data=means, x= "xcoord",y="ycoord",hue="zmeans",size="pcts",palette="Reds",sizes=(0, 250))
    ax.set_xticks(x_coords)
    ax.set_xticklabels(list(ordered_x_condition)*len(genes))
    ax.set_yticks(y_coords)
    ax.set_yticklabels(ordered_y_condition)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation = 90)
    ax.set_ylim((ax.get_ylim()[0]-.3,ax.get_ylim()[1]+.3))
    #ax.set_xlim((ax.get_xlim()[0]+(between_gene_sep-within_gene_sep),ax.get_xlim()[1]-(between_gene_sep-within_gene_sep)))
    for i,g in enumerate(genes):
        plt.text(gene_label_locs[i],ax.get_ylim()[1]+.3,g,horizontalalignment='center',multialignment="center")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    for xc in linepositions[:-1]:
        plt.axvline(x=xc, color='grey')
        
def boxplot_sample_proportions(adata, x_value, color_value,hue="treatment",figsize=(10,5), plottype="box",order=None,hue_order=None,edgecolor=False,swap=False):
    tmp = adata.obs.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)

    m=tmp.divide(tmp.sum(axis=1), axis=0)
    props = []

    i=0
    if hue+"_colors" in adata.uns and not swap:
        color_dict = dict(zip(adata.obs[hue].cat.categories,adata.uns[hue+"_colors"]))
    elif color_value+"_colors" in adata.uns and swap:
        color_dict = dict(zip(adata.obs[color_value].cat.categories,adata.uns[color_value+"_colors"]))
    else:
        color_dict=None
    for sample in m.index:

        for celltype in m.columns:
            vals = [sample,m.loc[sample,celltype],celltype,adata.obs.loc[adata.obs[x_value]==sample,hue].unique()[0]]
            props.append(vals)
            i+=1
    props_df = pd.DataFrame(props,columns=[x_value,x_value+"_proportion",color_value,hue])
    #sns.boxplot(x="celltype", y="sample_proportion", hue="treatment", data=tips)
    props_df[hue]=props_df[hue].astype("category")
    plt.figure(figsize=figsize)
    if swap:
        old_hue = hue
        hue=color_value
        color_value=old_hue
        old_hue_order=hue_order
        hue_order=order
        order=old_hue_order
    if plottype=="box":

        p=sns.boxplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df, palette=color_dict,hue_order=hue_order,linewidth=3)
        
        if edgecolor==True:
            for i,box in enumerate(p.artists):
                box.set_edgecolor(box.get_facecolor())
                r,g,b,a = box.get_facecolor()
                box.set_facecolor((r,g,b,.3))
            swarm_palette=color_dict
        else:
            swarm_palette={i:"white" for i in color_dict}
        if hue_order is None:
            hue_order = adata.obs[hue].cat.categories
        sns.swarmplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,  dodge=True,palette=swarm_palette,edgecolor="white",linewidth=.7,size=3.2, hue_order=hue_order)
        plt.legend(p.artists,hue_order)
    if plottype=="bar":
        p=sns.barplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,palette=color_dict,hue_order=hue_order)
    p.set_xticklabels(p.get_xticklabels(),rotation=90)



