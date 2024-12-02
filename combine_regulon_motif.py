"""
    The script is to combine individual Regulon analysis into a final one.
    1. The master transcription factor must appear at least 80% of the 
       predictions to ensure a high level of consensus.
    2. The target genes of the master TFs must appear at least five times to 
       make sure that they do not appear randomly.
"""

import json,zlib,base64,os,numba,pickle

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def get_genelist(flist='denovo_candidates.name.txt'):
    genelist = []
    lines = open(flist,'r')
    for line in lines:
        genelist.append(line.strip())
    return genelist

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
    
    def __repr__(self):
        return f"<Regulon name:{self.name}, gene2weight:{self.gene2weight}, gene2occurrence:{self.gene2occurrence}, transcription_factor:{self.transcription_factor}, context:{self.context}, score:{self.score}, nes:{self.nes}, occurrence:{self.occurrence}>"
    
    def update(self,inp):
        assert self.name == inp.name
        for gene in inp.gene2weight:
            if gene in self.gene2weight:
                if inp.gene2weight[gene] > self.gene2weight[gene]:
                    self.gene2weight[gene] = inp.gene2weight[gene]
                self.gene2occurrence[gene] += 1
            else:
                self.gene2weight[gene] = inp.gene2weight[gene]
                self.gene2occurrence[gene] = 1
                
        if inp.context in self.context:
                if inp.score>self.context[inp.context]:
                    self.context[inp.context] = inp.score
                    self.annotation[inp.context] = inp.annotation
        else:
            self.context[inp.context] = inp.score
            self.annotation[inp.context] = inp.annotation
        
        if inp.score > self.score:
            self.score = inp.score
            
        if inp.nes > self.nes:
            self.nes = inp.nes
        
        if inp.orthologous_identity > self.orthologous_identity:
            self.orthologous_identity = inp.orthologous_identity
        
        if inp.similarity_qvalue > self.similarity_qvalue:
            self.similarity_qvalue = inp.similarity_qvalue
            
        self.occurrence += 1
        

def combineRegulons(pref,flist,nrep=100,pcut=0.8):

    genes = get_genelist(flist)
    print(genes)

    ntrials = 0
    regulon_sets = []
    for i in range(1,nrep+1):
        f = "%s_%d.pickle"%(pref,i)
        if os.path.isfile(f):
            regulons_i = pickle.load(open(f,"rb"))
            regulon_sets += regulons_i
            ntrials += 1
        else:
            print(f"!Error> File {f} not found")

    regulons = {}
    for regi in regulon_sets:
        if regi.name in regulons:
            regulons[regi.name].update(regi)
        else:
            regulons[regi.name] = Regulon(regi)

    sorted_regulons = sorted(regulons.items(),key=lambda x:x[1].occurrence,
                                reverse=True)
    sorted_regulons = dict(sorted_regulons)
        
    ncut = ntrials*pcut
    
    regulons_final = {}
    for name in sorted_regulons:
        #print(name)
        #print(regulons[name])
        if regulons[name].occurrence >= ncut:
            regulons_final[name] = regulons[name]
            regulons_final[name].gene2occurrence = {k:v for k,v in regulons[name].gene2occurrence.items() if v>5}
            regulons_final[name].gene2weight = {k:v for k,v in regulons[name].gene2weight.items() if k in regulons_final[name].gene2occurrence}


    regulons_final = dict(sorted(regulons_final.items(),key=lambda x:x[1].occurrence,reverse=True))

    results = []
    for name in regulons_final:
        tf = regulons_final[name].transcription_factor
        gene2occurrence = ["%s|%d"%(k,v) for k,v in sorted(regulons_final[name].gene2occurrence.items(),key=lambda x:x[1],reverse=True)]
        num = len(regulons_final[name].gene2occurrence)
        g2o_str = ";".join(gene2occurrence)
        goi = [g for g in regulons_final[name].gene2occurrence if g in genes]
        goi_str = ";".join(goi)
        goi_num = len(goi)
        #print(",".join([name,str(num),g2o_str,goi_str,str(goi_num)]))
        results += [[tf,num,g2o_str,goi_str,goi_num]]

    results = sorted(results,key=lambda x:x[-1],reverse=True)
    f = open("%s_combined.csv"%pref,'w')
    f.write("tf,num_targets,gene2occurrence,gene_of_interest,num\n")
    for result in results:
        tf,num,g2o_str,goi_str,goi_num = result
        f.write(",".join([tf,str(num),g2o_str,goi_str,str(goi_num)])+"\n")
    f.close()

    with open('%s_combined.pickle'%pref, 'wb') as handle:
        pickle.dump(regulons_final, handle)
#print(reg)

if __name__ == "__main__":
    pref = "denovo_regulons"
    flist = "denovo_candidates.name.txt"
    combineRegulons(pref,flist,nrep=100,pcut=0.8)
