
import os,sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.decomposition import PCA,FastICA,SparsePCA
from sklearn.cluster import KMeans,AgglomerativeClustering

def read_loom(fdata):
    adata = sc.read_loom(fdata,validate=False)
    print('Finished loading data %s'%fdata)
    return adata

def get_tflist():
    DBDIR="/ru-auth/local/home/jpeng/scratch/softwares/pySCENIC"
    TFS = os.path.join(DBDIR,'allTFs_dmel.txt')
    tf_names = np.loadtxt(TFS,dtype=str).tolist()
    return tf_names

def get_denovo_genes(flist='denovo_candidates.names.txt'):
    genelist = []
    lines = open(flist,'r')
    for line in lines:
        genelist.append(line.strip())
    return genelist

def get_testis_biased_genes(flist="fca_testis_biased_genes.csv"):
    df = pd.read_csv(flist)
    testis_biased_genes = df["fca_sym"].tolist()
    return testis_biased_genes

def plot_cells(ratio,category="tissue"):
    fig,axes = plt.subplots(1,2,figsize=(5,3.6))
    ax = axes[0]
    sns.barplot(data=ratio,ax=ax,x="ratio",y="tissue",
                linewidth=0.66,color="grey")
    ax.axvline(x=0.21,color="blue",linestyle="--")
    ax.set_ylabel("")
    ax.set_xlabel("ratio of remaining cells")

    ax = axes[1]
    sns.barplot(data=ratio,ax=ax,x="number",
                y="tissue",linewidth=0.66,color="grey")
    #ax.axvline(x=0.21,color="blue",linestyle="--")
    #ax.set_xticks([1e2,1e3,1e4])
    ax.set_ylabel("")
    ax.set_yticklabels([])
    ax.set_xscale('log')
    ax.set_xlabel("number of remaining cells")

    fig.tight_layout()
    plt.savefig("submatrix_cell_ratio_in_%s.pdf"%category)
    return True

def stat_and_plot_cells(loom_data,expsub,category="tissue"):
    if category=="tissue":
        groups = np.unique(loom_data.obs.tissue)
        category = loom_data.obs.tissue
        ntot = len(loom_data.obs.tissue)
    elif category=="annotation_broad":
        groups = np.unique(loom_data.obs.annotation_broad)
        category = loom_data.obs.annotation_broad
        ntot = len(loom_data.obs.annotation_broad)

    ncells = {t:0 for t in groups}
    mcells = {t:0 for t in groups}

    mdict = {idx:1 for idx in expsub.index}
    for i in range(ntot):
        t = category[i]
        if i in mdict:
            mcells[t] += 1
        ncells[t] += 1

    ratio = {category:[],"ratio":[],"number":[]}
    for t in tissues:
        ratio[t] = mcells[t]/ncells[t]
        print(t,ncells[t],mcells[t],ratio[t])
        ratio[category].append(t)
        ratio["ratio"].append(ratio[t])
        ratio["number"].append(mcells[t])
 
        ratio = pd.DataFrame(ratio)

    plot_cells(ratio,category=category)


if __name__ == "__main__":

    fdata = '../../RAWDATA/r_fca_biohub_all_wo_blood_10x.loom'
    loom_data = read_loom(fdata)

    denovo_genes = get_denovo_genes()
    testis_biased = get_testis_biased_genes()

    ## subgenes ##
    tf_names = get_tflist()
    subgenes = []
    for g in tf_names + denovo_genes:
        if g in loom_data.var.index:
            subgenes.append(g)

    ## subgenes2 ##
    subgenes2 = {}
    for g in tf_names + denovo_genes + testis_biased:
        if g in loom_data.var.index:
            subgenes2[g] = 1
    subgenes2 = [g for g in subgenes2]

    ## get matrix ##
    m,n = loom_data.X.shape
    matrix = np.zeros((m,n))
    matrix += loom_data.X

    matrix = pd.DataFrame(matrix,
                            columns=loom_data.var.index,)
                            #index=loom_data.obs['annotation_broad'])
                            #index=loom_data.obs.index)

    expmatrix = matrix.loc[(matrix[valid_denovo_genes].max(axis=1)>=5) & (matrix[valid_tf_names].max(axis=1)>=20)]
    expsub = expmatrix[subgenes2].copy()
    expsub.to_csv("subgenes2.csv")

    stat_and_plot_cells(loom_data,expsub,category="tissue")
    stat_and_plot_cells(loom_data,expsub,category="annotation_broad")
