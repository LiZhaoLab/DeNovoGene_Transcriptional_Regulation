""" 
    De novo gene transcriptional regulation network by specific input TFs.
    For better visualizaiton, we only include the input TFs (1), their direct 
    de novo gene targets (1), and possible secondary TFs that are regulated
    directly by gene set (1), and also directly regulate gene set (2).
"""

import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import numpy as np
import os,sys,pickle,pygraphviz,random
from networkx.drawing.nx_agraph import graphviz_layout

class Regulon():
    def __init__(self,inp):
        self.name = inp.name
        self.gene2weight = {}
        for gene in inp.gene2weight:
            self.gene2weight[gene] = inp.gene2weight[gene]
        self.gene2occurrence = {}
        for gene in self.gene2weight:
            self.gene2occurrence[gene] = 1
        self.transcription_factor = inp.transcription_factor
        self.context = {inp.context:inp.score}
        self.score = inp.score
        self.nes = inp.nes
        self.orthologous_identity = inp.orthologous_identity
        self.similarity_qvalue = inp.similarity_qvalue
        self.annotation = {inp.context:inp.annotation}
        self.occurrence = 1

"""
    plot the regulatory network as a whole may result in a messy network that is difficult to interpret
"""
def get_regulon(tf,data):
    reg = {}
    for d in data:
        if tf == data[d].transcription_factor:
            reg = data[d]
    return reg

def get_regulated_denovo(reg):
    regulated = [g for g in reg.gene2occurrence if g in denovo_genes]
    return regulated

def get_regulated_tf(reg):
    regulated = [g for g in reg.gene2occurrence if g in tf_genes]
    return regulated

def get_genes(flist='denovo_candidates.names.txt'):
    genelist = []
    lines = open(flist,'r')
    for line in lines:
        genelist.append(line.strip())
    return genelist

def dump_genes(pref,genes,sym2fbgn):
    f = open("%s_fbid.txt"%pref,"w")
    for g in genes:
        fbid = sym2fbgn[g]
        f.write("%s\n"%fbid)
    f.close()

def is_reachable(G,g,tf="vis"):
    flag = False
    if g in G.nodes:
        if nx.has_path(G,g,tf):
            flag = True
    return flag

def direct_regulation(G,tf,g):
    flag = False
    if g in G.nodes:
        try:
            path = nx.shortest_path(G, source=tf, target=g)
            if len(path) == 2:
                flag == True
        except nx.NetworkXNoPath:
            pass
    return flag

def number_of_direct_regulators(G, targets):
    """
    Computes the number of direct regulators for a list of targets.
    The hypothesis is that, de novo genes has less direct regulators than other genes.
            TF1 --> Target
            TF2 --> Target
    
    Parameters:
        G (networkx.DiGraph): Directed graph representing the network.
        targets (list): List of target nodes.
        
    Returns:
        float: Average number of direct regulators.
    """
    
    t1 = time.time()
    num = []
    for node in targets:
        try:
            regulators = list(G.predecessors(node)) 
        except nx.NetworkXError as e:
            regulators = []
        num.append(len(regulators))
        
    ave_num = np.mean(num)
    t2 = time.time()
    
    return num,ave_num

def number_of_alternative_regulations(G, source, targets, cutoff = 2):
    """
    Computes the number of alternative path from source to targets.
    The hypothesis is that, de novo genes has less alternative regulations than other genes.
            Major --> TF1 --> Target
            Major --> TF2 --> Target
    
    Parameters:
        G (networkx.DiGraph): Directed graph representing the network.
        source (list): List of major regulons
        targets (list): List of target nodes.
        
    Returns:
        float: Average alternative regulations.
        
    Note: 
        Accessing all simple paths can be time consuming. In this case, I decided to ignore 
        complicated regulations by implementing a cutoff path length of 2 by default.
    """
    
    t1 = time.time()

    num = []
    for target in targets:
        for tf in source:
            paths = list(nx.all_simple_paths(G, tf, target, cutoff=2))
            num_paths = len(paths)
            if num_paths:
                num.append(num_paths)
    ave_num = np.mean(num)

    return num,ave_num

def num_regulators_groups(G,reg_tf,major_denovo,other_denovo,other_biased):
    """
        bootstrap analysis for number of direct TFs
    """
    data = {
            "major_denovo":[],
            "other_denovo":[],
            "other_biased":[],
            "random":[],
         }

    nsub = 100
    nstep = 50
    for i in range(nstep):
        sub_random = random.sample(G.nodes,nsub)
        sub_major_denovo = random.sample([g for g in major_denovo],nsub)
        sub_other_denovo = random.sample([g for g in other_denovo],nsub)
        sub_other_biased = random.sample([g for g in other_biased],nsub)
        sub_other_tf     = random.sample([g for g in fca_tf if g not in reg_tf],nsub)
        _,n0 = number_of_direct_regulators(G,sub_random)
        _,n1 = number_of_direct_regulators(G,sub_major_denovo)
        _,n2 = number_of_direct_regulators(G,sub_other_denovo)
        _,n3 = number_of_direct_regulators(G,sub_other_biased)
        _,n4 = number_of_direct_regulators(G,sub_other_tf)

        
        data["major_denovo"].append(n1)
        data["other_denovo"].append(n2)
        data["other_biased"].append(n3)
        data["random"].append(n0)

    data = pd.DataFrame(data)
    
    fig,ax = plt.subplots(figsize=(4,3))
    sns.violinplot(data,ax=ax,cut=0,scale="width",line_width=0.66,color="gray")
    ax.set_ylabel("No. direct TFs",fontsize=8)
    #ax.set_xticklabels([])
    #ax.set_yticks()
    fig.tight_layout()
    plt.savefig("network_num_direct_regulators_bootstrapped.pdf")
    
    data.to_csv("network_num_direct_regulators_bootstrapped.csv",index=False)
    return data

def num_regulations_groups(G,reg_tf,major_denovo,other_denovo,other_biased):
    """
        bootstrap analysis to compute number transcriptional regulation paths
    """
    data = {
            "major_denovo":[],
            "other_denovo":[],
            "other_biased":[],
            "random":[],
         }

    nsub = 100
    nstep = 50
    for i in range(nstep):
        t1 = time.time()
        sub_random = random.sample(G.nodes,nsub)
        sub_major_denovo = random.sample([g for g in major_denovo],nsub)
        sub_other_denovo = random.sample([g for g in other_denovo],nsub)
        sub_other_biased = random.sample([g for g in other_biased],nsub)
        sub_other_tf     = random.sample([g for g in fca_tf if g not in reg_tf],nsub)
        _,n0 = number_of_alternative_regulations(G,reg_tf,sub_random)
        _,n1 = number_of_alternative_regulations(G,reg_tf,sub_major_denovo)
        _,n2 = number_of_alternative_regulations(G,reg_tf,sub_other_denovo)
        _,n3 = number_of_alternative_regulations(G,reg_tf,sub_other_biased)
        _,n4 = number_of_alternative_regulations(G,reg_tf,sub_other_tf)

        
        data["major_denovo"].append(n1)
        data["other_denovo"].append(n2)
        data["other_biased"].append(n3)
        data["random"].append(n0)
        t2 = time.time()
        print(f"step {i+1: 2d}, time {t2-t1: .2f}s")

    data = pd.DataFrame(data)
    
    fig,ax = plt.subplots(figsize=(4,3))
    sns.violinplot(data,ax=ax,cut=0,scale="width",line_width=0.66,color="gray")
    ax.set_ylabel("No. regulations",fontsize=8)
    #ax.set_xticklabels([])
    #ax.set_yticks()
    fig.tight_layout()
    plt.savefig("network_num_regulations_bootstrapped.pdf")
    
    
    data.to_csv("network_num_regulations_bootstrapped.csv",index=False)
    return data



def regulons2network(regulons,fca_denovo,fca_tf,fca_biased,major):
    """
        analyze the whole network, not only de novo genes, but also others
    """
    attr_nodes = dict()
    attr_edges = dict()

    major_denovo = {}
    other_denovo = {}
    other_biased = {}

    reg_tf = {}
    for d in regulons:
        # TF name
        g1 = regulons[d].transcription_factor
        # The network
        reg_tf[g1] = 1
        for g2 in regulons[d].gene2occurrence:
                attr_nodes[g2] = 1
                attr_edges[(g1,g2)] = 1

        # de novo genes regulated by major TFs
        if g1 in major:
            novo = [g for g in regulons[d].gene2occurrence if g in fca_denovo]
            non_novo = [g for g in regulons[d].gene2occurrence if (g not in fca_denovo and g not in fca_tf)]
            for g2 in novo:
                major_denovo[g2] = 1
            for g2 in non_novo:
                other_biased[g2] = 1
        else:
            novo = [g for g in regulons[d].gene2occurrence if g in fca_denovo]
            non_novo = [g for g in regulons[d].gene2occurrence if (g in fca_biased and g not in fca_denovo)]
            for g2 in novo:
                other_denovo[g2] = 1
            for g2 in non_novo:
                other_biased[g2] = 1

    biased_list = [g for g in other_biased]
    for g in biased_list:
        if g in reg_tf:
            other_biased.pop(g)

    #print(major_denovo)
    ## create the graph ##
    G = nx.DiGraph()
    for g in attr_nodes:
        G.add_node(g)

    for edges in attr_edges:
        n1,n2 = edges
        G.add_edge(n1, n2)
        
    return G, major_denovo, other_denovo, other_biased, reg_tf


def regulon2network_denovo(regulon,fca_tf,fca_denovo,fca_biased,
                            major=["vis","Jra","achi"]):
    """
    plot the regulatory network, which include:
        1. major regulators
        2. directl regulated de novo gene candidates
        3. possible secondary regulators that contribute to the network
    """

    attr_nodes = dict()
    attr_edges = dict()

    major_denovo = {} ## a dictionary to store major de novo for fast lookup
    reg_tf = {} ## a dictionary to store possible secondary TFs for fast lookup

    ## loop through all regulons ##
    for d in regulon:
        # TF name
        g1 = regulon[d].transcription_factor
        
        ## only include g1 if g1 regulates de novo genes
        novo = [g for g in regulon[d].gene2occurrence if g in fca_denovo]
        
        # The network should include de novo genes as targets ##
        if len(novo)>0:
            reg_tf[g1] = 1
            for g2 in regulon[d].gene2occurrence:
                    attr_nodes[g2] = 1
                    attr_edges[(g1,g2)] = 1

            # de novo genes regulated by major TFs
            if g1 in major:
                novo = [g for g in regulon[d].gene2occurrence if g in fca_denovo]
                for g2 in novo:
                    major_denovo[g2] = 1
                            
    #print(major_denovo)
    ## create the directed graph ##
    G = nx.DiGraph()
    for g in attr_nodes:
        G.add_node(g)
            
    for edges in attr_edges:
        n1,n2 = edges
        G.add_edge(n1, n2)
    

    ## For better visualizaiton, we only include the input TFs (1), their direct
    ## de novo gene targets (1), and possible secondary TFs that are regulated
    ## directly by gene set (1), and also directly regulate gene set (2).

    ## This step is to find out these secondary TFs
    ## makesure that max_steps is set to 1.
    max_steps = 1
    subgraph_nodes = set(major_denovo)  # Start with gene set (2)
    for node in major_denovo:
        for tf in major:
        # Perform BFS from the major regulator and capture nodes within `max_steps` steps
            paths = nx.single_source_shortest_path_length(G, tf, cutoff=max_steps)
        
            # Include only nodes that can reach the target node
            if node in paths:
                subgraph_nodes.update(g for g, d in paths.items() if d <= max_steps and g in reg_tf)
        
    ## now create the subgraph ##
    subG = G.subgraph(subgraph_nodes)
    
    ## to create different size and colors for different categories of nodes
    for node in subG.nodes:
        if node in major:
            subG.nodes[node]["category"] = "major"
        elif node in fca_tf:
            subG.nodes[node]["category"] = "secondary"
        else:
            subG.nodes[node]["category"] = "denovo"
    
    # custom color palette
    palette = ['#fe218b', '#fed700', '#b8b8ff',"#eaebed"]
    
    color_map["major"] = palette[0]
    color_map["secondary"] = palette[1]
    color_map["denovo"] = palette[2]
    node_colors = [color_map[subG.nodes[node]["category"]] for node in subG.nodes]
    # Node sizes based on categories
    size_map["major"] = 360
    size_map["secondary"] = 150
    size_map["denovo"] = 150
    node_sizes = [size_map[subG.nodes[node]["category"]] for node in subG.nodes]

    ## layout of the graph ##
    pos = graphviz_layout(subG)
    
    # Draw the network
    #fig,ax = plt.subplots(figsize=(10, 10))
    ## adjust figure size based on number of TFs we want to check 
    if len(major)>60:
        width = 10
    elif len(major)>25:
        width = 8
    else:
        width = 5
    fig,ax = plt.subplots(figsize=(width, width))

    # Draw nodes in background #
    nodelist=[node for node in subG.nodes if node in major_denovo]
    node_colors = [color_map[subG.nodes[node]["category"]] for node in nodelist]
    node_sizes = [size_map[subG.nodes[node]["category"]] for node in nodelist]
    nx.draw_networkx_nodes(subG, pos, nodelist=nodelist,node_color=node_colors,node_size=node_sizes)
    
    # Draw secondary in second layer
    nodelist=[node for node in subG.nodes if (node not in major_denovo and node not in major)]
    node_colors = [color_map[subG.nodes[node]["category"]] for node in nodelist]
    node_sizes = [size_map[subG.nodes[node]["category"]] for node in nodelist]
    nx.draw_networkx_nodes(subG, pos, nodelist=nodelist,node_color=node_colors,node_size=node_sizes)
    
    # Draw major regulators in foreground
    nodelist=[node for node in subG.nodes if node in major]
    node_colors = [color_map[subG.nodes[node]["category"]] for node in nodelist]
    node_sizes = [size_map[subG.nodes[node]["category"]] for node in nodelist]
    nx.draw_networkx_nodes(subG, pos, nodelist=nodelist,node_color=node_colors,node_size=node_sizes)
    
    nx.draw_networkx_labels(subG,pos,font_size=5,
                            labels={n:n for n in subG.nodes if n not in major})
    nx.draw_networkx_labels(subG,pos,font_size=8,
                            labels={n:n for n in subG.nodes if n in major})
    nx.draw_networkx_edges(subG,pos,
                           arrowstyle="-|>",arrowsize=5,
                           edge_color="gray",alpha=1,width=0.66,
                           edgelist=[(u, v) for u, v in subG.edges],
                            )
    plt.axis('off')
    fig.tight_layout()
    
    if len(major)<10:
        plt.savefig("%s_network.pdf"%"_".join(major))
    elif len(major)<51:
        plt.savefig("top31_network.pdf")
    else:
        plt.savefig("top83_network.pdf")
    plt.show()
    return G



if __name__ =="__main__":

    df = pd.read_csv("all_fca_genes.validation.csv")
    sym2fbgn = dict(zip(df["fca_sym"],df["FBgn"]))

    df = pd.read_csv("fca_testis_biased_genes.csv")
    biased_genes = df["fca_sym"].to_list()
    fca_biased = [g for g in biased_genes if g in sym2fbgn]

    denovo_genes = get_genes("fdr000001_candidates_20230502.name.txt")
    tf_genes = get_genes("allTFs_dmel.txt")

    fca_denovo = [g for g in denovo_genes if g in sym2fbgn]
    fca_tf = [g for g in tf_genes if g in sym2fbgn]

    with open("denovo_regulons_combined.pickle", "rb") as p:
        data = pickle.load(p)

    ## plot de novo gene transcriptional regulatory network by vis
    major = ["vis"]
    G = regulon2network_denovo(data,fca_tf,fca_denovo,fca_biased,major=major)

    ## plot de novo gene transcriptional regulatory network by vis/achi
    major=["vis","achi"]
    G = regulon2network_denovo(data,fca_tf,fca_denovo,fca_biased,major=major)

    ## plot de novo gene transcriptional regulatory network by vis/achi/Jra
    major=["vis","achi","Jra"]
    G = regulon2network_denovo(data,fca_tf,fca_denovo,fca_biased,major=major)



    ## transcriptional regulation network analysis ##
    major = ["vis","achi","Jra"]
    G, major_denovo, other_denovo, other_biased, reg_tf = regulons2network(data,fca_denovo,fca_tf,fca_biased,major)
    num_regulators  = num_regulators_groups(G,reg_tf,major_denovo,other_denovo,other_biased)
    num_regulations,= num_regulations_groups(G,reg_tf,major_denovo,other_denovo,other_biased)
