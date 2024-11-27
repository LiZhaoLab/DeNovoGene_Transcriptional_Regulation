
import os,sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

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

def filter_matrix(loom_data,tf_names,denovo_genes,testis_biased):
    ## subgenes ##
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

    df_mat = pd.DataFrame(matrix,
                            columns=loom_data.var.index,)
                            #index=loom_data.obs['annotation_broad'])
                            #index=loom_data.obs.index)

    expmat = expmat.loc[(df_mat[valid_denovo_genes].max(axis=1)>=5) 
                        & (df_mat[valid_tf_names].max(axis=1)>=20)]
    expsub = expmat[subgenes2].copy()
    expsub.to_csv("subgenes2.csv")

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

def matrix_by_feature(loom_data,matrix,denovo_genes,
                      feature='tissue',dtype="loom"):
    """
        Z-score by feature
    """
    index = adata.obs[feature]
    if dtype=="h5ad":
        genes = adata.raw.var.index
    elif dtype=="loom":
        genes = adata.var.index

    m,n = X.shape

    df = pd.DataFrame(matrix,columns=genes,index=index)
    groups = [s for s in adata.obs[feature].unique() if s!='unannotated' and s!='artefact']

    ## expression matrix ##
    ndenovo = len(denovo_genes)
    dfm = []
    for grp in groups:
        sub_df  = df.loc[grp,:]
        sub_mat = np.zeros(sub_df.shape)
        sub_mat += sub_df
        #print(sub_mat)
        print(grp,sub_mat.shape)

        #dfm[grp] = sub_mat
        sub_sum = sub_mat.sum().sum()
        expscale = np.zeros(ndenovo)
        for i in range(ndenovo):
            gene = denovo_genes[i]
            if gene in genes:
                datai = np.array(sub_mat[gene])
                #print(np.sum(datai),sub_sum)
                expscale[i] = (np.sum(datai)*100)/sub_sum
        dfm.append(expscale)
    dfm = np.array(dfm)
    dfm = pd.DataFrame(dfm,index=groups,columns=denovo_genes)
    
        ## construct zscore matrix ##
    expressed_genes = []
    for gene in denovo_genes:
        if np.sum(dfm[gene])>0:
            expressed_genes.append(gene)
    mat = dfm[expressed_genes]

    means = mat.mean(axis=0)
    var_mat  = mat-means
    var_mat2 = np.square(var_mat)
    print('Finished variance matrix calculation')
    M,N = mat.shape

    ## construct zscore matrix
    var_sum = var_mat2.sum(axis=0)
    std_column = np.sqrt(var_sum/M)
    print('Finished standard deviation calculation')

    zscore_matrix = (mat-means)/std_column
    print('Finished Zscore calculation')
    zscore_matrix.to_csv("denovo_zmat_%s.csv"%feature)

    return zscore_matrix



if __name__ == "__main__":

    fdata = './RAWDATA/r_fca_biohub_all_wo_blood_10x.loom'
    loom_data = read_loom(fdata)

    ## filter matrix
    tf_names = get_tflist()
    denovo_genes = get_denovo_genes()
    testis_biased = get_testis_biased_genes()

    matrix,expsub  = filter_matrix(loom_data,tf_names,denovo_genes,
                                    testis_biased)

    ## stat cell numbers and ratios
    stat_and_plot_cells(loom_data,expsub,category="tissue")
    stat_and_plot_cells(loom_data,expsub,category="annotation_broad")

    ## Z-scores across different tissues/celltypes
    zscore_annotation = matrix_by_feature(loom_data,matrix,denovo_genes,feature="annotation",dtype="loom")
