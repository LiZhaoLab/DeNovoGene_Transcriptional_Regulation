"""
    tissue specificity based on randomly selected genes with similar expression levels

    Here, I used a genetic algorithm with simulated annealing (SA) 
    to select genes with similar expression levels.

    The distribution difference is computed as l1norm of density, 
    which is the direct differences between the distribution of de novo genes 
    and the distribution of selected genes.

    L1_norm = sum(|Distribution(de novo)-Distribution(selected)|)
"""

## now compute the ratio of strongly biased genes in 
## randomly selected de novo genes and other genes

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os,sys

fontfile = "/ru-auth/local/home/jpeng/.cache/matplotlib/fontlist-v330.json"
if os.path.isfile(fontfile):
    print(fontfile)
    os.remove("/ru-auth/local/home/jpeng/.cache/matplotlib/fontlist-v330.json")
import matplotlib.font_manager
matplotlib.get_cachedir()
matplotlib.rcParams['font.family'] = ['arial']

def get_num_biased(zmat,zcut,fbgn2fca,fbidlist):
    num = 0
    tot = 0
    for fbid in fbidlist:
        try:
            g = fbgn2fca[fbid]
            z = zmat[g]
            zmax = np.max(z)
            if zmax>=zcut:
                num += 1
            tot += 1
        except KeyError:
            #print(fbid,g,zmax)
            pass
    return num*1.0/tot

def get_genes(flist='/ru-auth/local/home/jpeng/cactus/fdr000001_candidates.names.txt'):
    genelist = []
    lines = open(flist,'r')
    for line in lines:
        genelist.append(line.strip())
    return genelist

def convert_fca2fbgn(f="all_fca_genes.validation.csv"):
    df = pd.read_csv(f,dtype=str,na_filter=False)
    fca2fbgn = dict(zip(df["fca_sym"],df["FBgn"]))
    fbgn2fca = dict(zip(df["FBgn"],df["fca_sym"]))
    return fca2fbgn,fbgn2fca

denovo_genes = get_genes("../compare_tau/fdr000001_candidates_20230502.fbid.txt")

fca2fbgn,fbgn2fca = convert_fca2fbgn()
fca2fbgn,fbgn2fca = convert_fca2fbgn()
fca_genes_in_flyatlas = [fbid for fbid in fbgn2fca if fbid in fbgn2exp and fbid not in denovo_genes]
denovo_genes_in_flyatlas = [fbid for fbid in denovo_genes if fbid in fbgn2exp and fbid in fbgn2fca]

## Simulated Annealing Parameters ##

nsteps = 100
max_iter = 2000
num = len(denovo_genes_in_flyatlas)
print(num)
K = 0.95 # after each iteration, temperature cools down
T0 = 100. # initial temprature
bins = np.arange(0,16.5,.5)
l1norm_min_list = []
l1norm_max_list = []

for n in range(nsteps):
    list_dng = random.choices(denovo_genes_in_flyatlas,k=num)
    list_rng = random.choices(fca_genes_in_flyatlas,k=num)
    list_rng0 = [g for g in list_rng]
    dng_arr = [fbgn2exp[fbid] for fbid in list_dng]
    rng_arr = [fbgn2exp[fbid] for fbid in list_rng]
    dng_hist,edges = np.histogram(dng_arr, bins=bins, range=None, density=True, weights=None)
    rng_hist,edges = np.histogram(rng_arr, bins=bins, range=None, density=True, weights=None)
    l1norm0 = sum(abs(rng_hist-dng_hist))
    l1norm_min = 1e9
    min_iter = 0
    list_min = []
    rng_hist_min = []

    if n==0:
        hist2plot_dng  = [d for d in dng_hist]
        hist2plot_rng0 = [d for d in rng_hist]
    
    print(f"step {n}","iter 0",l1norm0)
    existed = {g:1 for g in list_rng}
    l1norm_max_list.append(l1norm0)

    # Actual Simulated Annealing #
    T = T0
    for it in range(max_iter):
        # decide which one to mutate #
        to_mut = random.randint(0,num-1)
        old_gene = list_rng[to_mut]
        gene_pool = [g for g in fca_genes_in_flyatlas if g not in existed]
        mut_gene = random.choice(gene_pool)
        
        ## mutate ##
        list_rng[to_mut] = mut_gene
        existed[mut_gene] = 1
        rng_arr = [fbgn2exp[fbid] for fbid in list_rng]
        rng_hist,edges = np.histogram(rng_arr, bins=bins, range=None, density=True, weights=None)
        l1norm1 = sum(abs(rng_hist-dng_hist))
        
        ## accept if l1norm decrease ##
        paccpt = 1.0 if l1norm1<l1norm0 else np.exp(-(l1norm1-l1norm0)/T)
        p = random.uniform(0,1)
        if p < paccpt:
            l1norm0 = l1norm1
            if l1norm0 < l1norm_min:
                l1norm_min = l1norm0
                list_min = [g for g in list_rng]
                rng_hist_min = [d for d in rng_hist]
                min_iter = it
        else:
            ## reject change
            list_rng[to_mut] = old_gene

        T = T*K

    if n==0:
        hist2plot_rng1 = [d for d in rng_hist_min]

    print(f"step {n}",f"iter {min_iter}",l1norm_min)
    l1norm_min_list.append(l1norm_min)

    np.savetxt("genes_with_similar_expression/genes_randomly_selected_%d.txt"%(n+1),list_rng,fmt="%s")
    np.savetxt("genes_with_similar_expression/genes_randomly_selected_init_%d.txt"%(n+1),list_rng0,fmt="%s")
    np.savetxt("genes_with_similar_expression/denovo_genes_randomly_selected_%d.txt"%(n+1),list_dng,fmt="%s")


## plot distributions, L1 norm (manhattan distance), etc... ###
fig,axes = plt.subplots(1,2,figsize=(5,2))
ax = axes[0]
ax.plot(bins[:-1],hist2plot_dng,color="black",linestyle=":",label="De novo genes")
ax.plot(bins[:-1],hist2plot_rng0,color="red",alpha=0.5,label="Other genes before SA")
ax.plot(bins[:-1],hist2plot_rng1,color="blue",alpha=0.5,label="Other genes after SA")
ax.tick_params(axis='both', labelsize=6)
ax.set_xlabel("log2(FPKM+1)",fontsize=8)
ax.set_ylabel("Density",fontsize=8)
ax.legend(fontsize=6)

ax = axes[1]
ax.plot(l1norm_max_list,alpha=0.5,color="red",label="L1-Norm before SA")
ax.plot(l1norm_min_list,alpha=0.5,color="blue",label="L1-Norm after SA")
ax.tick_params(axis='both', labelsize=6)
#ax.xlabel("iterations")
ax.set_xlabel("Number of SA",fontsize=8)
ax.set_ylabel("Differences (L1-Norm)",fontsize=8)
ax.legend(fontsize=6)
fig.tight_layout()
plt.savefig("genes_with_similar_expression.pdf")
plt.show()

## calculate ratio of specifically expressed genes #
# read in z-score matrix
z_tissue = pd.read_csv("all_zmat_tissue.csv",index_col=0)
z_cell = pd.read_csv("all_zmat_annotation_broad.csv",index_col=0)

ratio_tissue = {"ratio":[],"category":[]}
ratio_cell   = {"ratio":[],"category":[]}

zcut = 3.
nsteps = 100
for n in range(nsteps):
    genes_dng = np.loadtxt("genes_with_similar_expression/denovo_genes_randomly_selected_%d.txt"%(n+1),dtype=str)
    genes_rng = np.loadtxt("genes_with_similar_expression/genes_randomly_selected_%d.txt"%(n+1),dtype=str)
    genes_ini = np.loadtxt("genes_with_similar_expression/genes_randomly_selected_init_%d.txt"%(n+1),dtype=str)

    num_biased_tissue_dng = get_num_biased(z_tissue,zcut,fbgn2fca,genes_dng)
    num_biased_tissue_ini = get_num_biased(z_tissue,zcut,fbgn2fca,genes_ini)
    num_biased_tissue_rng = get_num_biased(z_tissue,zcut,fbgn2fca,genes_rng)
    num_biased_cell_dng = get_num_biased(z_cell,zcut,fbgn2fca,genes_dng)
    num_biased_cell_rng = get_num_biased(z_cell,zcut,fbgn2fca,genes_rng)
    num_biased_cell_ini = get_num_biased(z_cell,zcut,fbgn2fca,genes_ini)

    ratio_tissue["ratio"].append(num_biased_tissue_dng)
    ratio_tissue["category"].append("de novo genes")
    ratio_tissue["ratio"].append(num_biased_tissue_rng)
    ratio_tissue["category"].append("other genes")
    ratio_tissue["ratio"].append(num_biased_tissue_ini)
    ratio_tissue["category"].append("other genes 0")

    ratio_cell["ratio"].append(num_biased_cell_dng)
    ratio_cell["category"].append("de novo genes")
    ratio_cell["ratio"].append(num_biased_cell_rng)
    ratio_cell["category"].append("other genes")
    ratio_cell["ratio"].append(num_biased_cell_ini)
    ratio_cell["category"].append("other genes 0")

ratio_tissue = pd.DataFrame(ratio_tissue)
ratio_cell   = pd.DataFrame(ratio_cell)
print(ratio_tissue.loc[ratio_tissue.category=="de novo genes"].describe())
print(ratio_tissue.loc[ratio_tissue.category=="other genes"].describe())
print(ratio_cell.loc[ratio_cell.category=="de novo genes"].describe())
print(ratio_cell.loc[ratio_cell.category=="other genes"].describe())

## plot ratio of specifically expressed genes ##
fig,axes = plt.subplots(1,2,figsize=(3.9,1.5))
blue = (113,156,179)
gray = (189,189,189)
orange = (228,134,70)
ax = axes[0]
sns.violinplot(ratio_tissue,x="category",y="ratio",cut=0,linewidth=0.66,
               ax=ax,palette=['#%02x%02x%02x'%blue,'#%02x%02x%02x'%gray,'#%02x%02x%02x'%orange])
ax.tick_params(axis='both', labelsize=6)
ax.set_xlabel("Tissue specificity",fontsize=8)
ax.set_ylim(0.4,1.01)
#ax.set_ylim(0,1)
ax.set_yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.set_ylabel("Ratio",fontsize=8)

ax = axes[1]
sns.violinplot(ratio_cell,x="category",y="ratio",cut=0,linewidth=0.66,
               ax=ax,palette=['#%02x%02x%02x'%blue,'#%02x%02x%02x'%gray,'#%02x%02x%02x'%orange])
ax.tick_params(axis='both', labelsize=6)
ax.set_xlabel("Cell type specificity",fontsize=8)
ax.set_ylim(0.4,1.01)
#ax.set_ylim(0,1)
ax.set_yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.set_ylabel("")
ax.set_yticklabels([])
fig.tight_layout()
plt.savefig("genes_with_similar_expression_stronglyBiasedRatio.pdf")
plt.show()
