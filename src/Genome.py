#/usr/bin/env python
#coding:utf-8

import pandas as pd
import numpy as np


class Genome(object):
    """class containing genome of target individual"""

    def __init__(self, Size, path_init):
        self.Size = Size
        self.sec_dist = 60 # security distance of prot and genes
        
        TSS = pd.read_table(path_init+'/TSS.dat', header=0)
        TTS = pd.read_table(path_init+'/TTS.dat', header=0)
        self.gene_info = TSS.merge(TTS)
        self.prot = pd.read_table(path_init+'/prot.dat', header=0)
        self.env = pd.read_table(path_init+'/environment.dat')

    def __str__(self):
        s = "Genome Info :\n--------------\n"
        s += "Size : "+str(self.Size)+"\n"
        s += "Gene Info : \n"+str(self.gene_info)+"\n\n"
        s += "Barriers : \n"+str(self.prot)+"\n\n"    
        return s


    def nearest_obj_distance(self, p):
        """Needed to ensure that genes and barriers don't become less that sec_dist apart."""
        # dist gene, dist barrier
        d_g = self.Size
        d_b = self.Size

        for i, row in self.gene_info.iterrows():
            if abs(row['TSS_pos'] - p) < d_g:
                d_g = abs(row['TSS_pos'] - p)
            if abs(row['TTS_pos'] - p) < d_g:
                d_g = abs(row['TTS_pos'] - p)

        for i, row in self.prot.iterrows():
            if abs(row['prot_pos'] - p) < d_b:
                d_b = abs(row['prot_pos'] - p)

        
        ####### for extremities of genome cases :

        # max pos of a gene
        gmax = max(max(self.gene_info['TTS_pos']),max(self.gene_info['TSS_pos'])) 
        #idem min
        gmin = min(min(self.gene_info['TTS_pos']),min(self.gene_info['TSS_pos']))


        if p > gmax :
            if self.Size-p+gmin < d_g:
                d_g = self.Size-p+gmin

        if p < gmin :
            if p+self.Size-gmax < d_g:
                d_g = p+self.Size-gmax


        # same for barriers
        bmax = max(self.prot['prot_pos']) 
        bmin = min(self.prot['prot_pos']) 

        if p > bmax :
            if self.Size-p+bmin < d_b:
                d_b = self.Size-p+bmin

        if p < bmin :
            if p+self.Size-bmax < d_b:
                d_b = p+self.Size-bmax

        return d_g, d_b


    def good_inv_pos(self, s, e):
        """checks is inversion positions are ok."""
        
        # check if pos is in gene
        if self.pos_in_gene(s) or self.pos_in_gene(e) :
            return False 

        # check the distance between barriers and gene
        s_gd, s_bd = self.nearest_obj_distance(s)
        e_gd, e_bd = self.nearest_obj_distance(e)
        
        return True if s_gd+e_bd > self.sec_dist and s_bd+e_gd > self.sec_dist else False

    def pos_in_gene(self, p):
        """Returns True if position is in a gene"""

        for i, row in self.gene_info.iterrows():
            if p>row['TSS_pos'] and p<row['TTS_pos']:
                return True
        return False


    def get_inv_pos(self):
        """Returns 2 good positions for inversion"""

        s,e = np.sort(np.random.randint(0, self.Size, size=2))
        i=0
        while not self.good_inv_pos(s,e):
            s,e = np.sort(np.random.randint(0, self.Size, size=2))
            i+=1
        
        return s,e
    
    def find_genes(self,s,e):
        """Returns list of all genes between positions s and e"""
        
        genes_list = []
        for index,row in self.gene_info.iterrows(): #iterate over all start positions
            start_gene = row["TSS_pos"]
            if start_gene >= s and start_gene <= e:
                genes_list.append(row)
        #no need to check for end of genes because find_genes is used only
        #for inversions and inversions will never cut a gene
        return genes_list
            
    
    def inversion(self):
        """inversion between two good positions given by get_inv_pos"""
        
        start_inv,end_inv = self.get_inv_pos()
        genes_to_inv = self.find_genes(start_inv,end_inv)
        print("Before:")
        print(self.gene_info)
        print("Inversion starts at "+str(start_inv)+" and ends at "+str(end_inv))
        for gene in genes_to_inv:
            ind = gene["TUindex"]
            
            #change gene start and gene end
            self.gene_info.at[ind,"TSS_pos"] = end_inv - (gene["TSS_pos"] - start_inv)
            self.gene_info.at[ind,"TTS_pos"] = start_inv + end_inv - gene["TTS_pos"]
            
            #change gene orientation
            if gene["TUorient"] == "+":
                self.gene_info.at[ind,"TUorient"] = "-"
            else:
                self.gene_info.at[ind,"TUorient"] = "+"
        
        #reorder gene_info so genes will be in ascending order of gene start
        first_gene = genes_to_inv[0]["TUindex"]
        last_gene = genes_to_inv[-1]["TUindex"]
        self.gene_info = pd.concat([self.gene_info.iloc[0:first_gene],self.gene_info.iloc[last_gene:first_gene-1:-1],self.gene_info[last_gene+1:]])
        self.gene_info["TUindex"] = list(range(len(self.gene_info["TUindex"])))
        print("after:")
        print(self.gene_info)




        
        
# class Gene:
#     """docstring for Genome"""

#     def __init__(self, s,e,ori,prom_str):
#         self.s = s
#         self.e = e
#         self.ori = ori
#         self.prom_str = prom_str