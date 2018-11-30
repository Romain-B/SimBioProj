#/usr/bin/env python
#coding:utf-8

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import simulation2 as sim
import os
import subprocess


class Genome(object):
  """class containing genome of target individual"""

  def __init__(self, Size, path_init):
    self.Size = Size
    self.sec_dist = 60 # security distance of prot and genes
    self.indel_var = 0 # variance of indel sizes
    self.p_inv = .1 # p_inser = (1-p_inv)/2 = p_del
    self.generation = 0 # generation counter
    self.evs = [None] # events occurred by generation

    TSS = pd.read_table(path_init+'/TSS.dat', header=0)
    TTS = pd.read_table(path_init+'/TTS.dat', header=0)
    self.gene_info = TSS.merge(TTS)
    self.prot = pd.read_table(path_init+'/prot.dat', header=0)
    self.env = pd.read_table(path_init+'/environment.dat', header=None)
    self.path_init = path_init
    #self.path_to_simulator = path_to_simulator

    self.fit_bygeneration = []
    self.T0 = 0.1 # linked to fitness decrease prob. 

    print(subprocess.check_output(['cp', path_init+'/params.ini', '../sim_files/current/params.ini']))
    # running for time 0
    self.write_sim_files("../sim_files/current")
    ## Run simulation
    sim.start_transcribing('../sim_files/current/params.ini', "../sim_files/future")
    self.fit_bygeneration.append(self.fitness("../sim_files/future/save_tr_nbr.csv"))


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
      if row["TUorient"] == "+":
        if p>row['TSS_pos'] and p<row['TTS_pos']:
          return True
      else:
        if p>row['TTS_pos'] and p<row['TSS_pos']:
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
  

  def frag_length(self):
    """Returns size of fragment for deletion"""
        
    var = 0 if self.indel_var==0 else np.random.randint(-self.indel_var,self.indel_var)
    return self.sec_dist + var


  def insertion(self):
    """insertion of fragment of random length"""
        
    p = np.random.randint(0, self.Size)

    while self.pos_in_gene(p) :
      p = np.random.randint(0, self.Size)

    l = self.frag_length()

    # shift genes and barriers
    for i, row in self.gene_info.iterrows():
      if row['TSS_pos'] > p:
        self.gene_info.iloc[i, self.gene_info.columns.get_loc('TSS_pos')] += l
      if row['TTS_pos'] > p:
        self.gene_info.iloc[i, self.gene_info.columns.get_loc('TTS_pos')] += l

    for i, row in self.prot.iterrows():
      if row['prot_pos'] > p:
        self.prot.iloc[i, self.prot.columns.get_loc('prot_pos')]+= l

    self.Size += l
    #return p,l


  def barr_between_pos(self, s, e):
    """Returns True if there is a barrier between argument positions"""
        
    for i, row in self.prot.iterrows():
      if row['prot_pos'] > s and row['prot_pos'] <= e:
        return True
    return False
  
  def get_barr_between_pos(self, s, e):
    """Returns indices of barriers between argument positions"""
        
    barrs = []
    for i, row in self.prot.iterrows():
      if row['prot_pos'] > s and row['prot_pos'] <= e:
        barrs.append(i)
    return barrs

  def deletion(self):
    """deletion of fragment"""
        
    l = self.frag_length()

    # find good del position
    p1 = np.random.randint(0, self.Size)
    p2 = p1 + np.random.choice([-1,1])*l
    s, e = np.sort([p1,p2]) 

    while not self.good_inv_pos(s,e) and not self.barr_between_pos(s,e):
      p1 = np.random.randint(0, self.Size)
      p2 = p1 + np.random.choice([-1,1])*l
      s, e = np.sort([p1,p2]) 

    # shift genes and barriers
    for i, row in self.gene_info.iterrows():
      if row['TSS_pos'] > s:
        self.gene_info.iloc[i, self.gene_info.columns.get_loc('TSS_pos')] -= l
      if row['TTS_pos'] > s:
        self.gene_info.iloc[i, self.gene_info.columns.get_loc('TTS_pos')] -= l

    for i, row in self.prot.iterrows():
      if row['prot_pos'] > s:
        self.prot.iloc[i, self.prot.columns.get_loc('prot_pos')]-= l

    self.Size -= l
    #return s,l

  def find_genes(self,s,e):
    """Returns list of all genes between positions s and e"""
    
    genes_list = []
    for index,row in self.gene_info.iterrows(): #iterate over all start positions
      start_gene = row["TSS_pos"]
      if start_gene >= s and start_gene <= e:
        genes_list.append(index)
    #no need to check for end of genes because find_genes is used only
    #for inversions and inversions will never cut a gene
    return genes_list
          
  
  def inversion(self):
    """inversion between two good positions given by get_inv_pos"""
    
    start_inv,end_inv = self.get_inv_pos()
    genes_to_inv = self.find_genes(start_inv,end_inv)
    genes_copy = self.gene_info.copy()
    if len(genes_to_inv) :
      for ind in genes_to_inv:
        gene = self.gene_info.iloc[ind,:].copy()
        #change gene start and gene end
        genes_copy.at[ind,"TSS_pos"] = end_inv - (gene["TSS_pos"] - start_inv)
        genes_copy.at[ind,"TTS_pos"] = start_inv + end_inv - gene["TTS_pos"]
        
        #change gene orientation
        if gene["TUorient"] == "+":
          genes_copy.at[ind,"TUorient"] = "-"
        else:
          genes_copy.at[ind,"TUorient"] = "+"
      
      #reorder gene_info so genes will be in ascending order of gene start
      self.gene_info = genes_copy.sort_values(["TSS_pos"]).reset_index(drop=True)
      
      #change position of barriers
      for barr_ind in self.get_barr_between_pos(start_inv,end_inv):
        self.prot.at[barr_ind,"prot_pos"] = end_inv - (self.prot.iloc[barr_ind,1] - start_inv)
      #reorder barriers
      self.prot = self.prot.sort_values(["prot_pos"]).reset_index(drop=True)
      
  def write_sim_files(self, path_to_sim):
    #TSS
    self.gene_info.to_csv(path_to_sim+'/TSS.dat', sep="\t", columns=['TUindex', 'TUorient', 'TSS_pos', 'TSS_strength'], index=False)
    self.gene_info.to_csv(path_to_sim+'/TTS.dat', sep="\t", columns=['TUindex', 'TUorient', 'TTS_pos', 'TTS_proba_off'], index=False)
    self.prot.to_csv(path_to_sim+'/prot.dat', sep="\t", index=False)

    with open(self.path_init+'/tousgenesidentiques.gff', 'r') as f:
      gff = f.readlines()

    gff[3] = '##sequence-region\ttousgenesidentiques\t1\t'+str(self.Size)+'\n'
    gff[4] = 'tousgenesidentiques\tRefSeq\tregion\t1\t'+str(self.Size)+'\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n'


    with open(path_to_sim+'/tousgenesidentiques.gff', 'w') as f:
      f.writelines(gff)


  def run_generation(self, path_to_sim, plot_it=True):
    # EXACTLY ONE event per generation, determined by prob p_inv.
    
    print("\n\n----\n")

    #keep old info in case mutation is BAAAAAD.
    old_gene_info = self.gene_info.copy()
    old_prot = self.prot.copy()

    #choose event
    event = np.random.choice(['inv', 'ins', 'del'], p=[self.p_inv, (1-self.p_inv)/2, (1-self.p_inv)/2])

    #do event
    if 'inv' == event:
      self.inversion()
      print("Did an inversion.")
    elif 'ins' == event :
      self.insertion()
      print("Did an insertion.")
    elif 'del' == event:
      self.deletion()
      print("Did a deletion.")
    self.evs.append(event)

    self.write_sim_files(path_to_sim)
    
    # update gen counter
    self.generation += 1
    print("Now at generation", self.generation)

    # run simulation
    sim.start_transcribing('../sim_files/current/params.ini', "../sim_files/future")
    fit = self.fitness("../sim_files/future/save_tr_nbr.csv")

    # keep or not
    # if new fitness is inferior to old fit. 
    self.fit_bygeneration.append(fit)

    if fit < self.fit_bygeneration[self.generation-1]:
      p = np.random.exponential(self.T0)
      # with proba p, keep it. So if rd>p, reset
      if np.random.random() > p :
        self.gene_info = old_gene_info
        self.prot = old_prot
        self.fit_bygeneration[self.generation] = self.fit_bygeneration[self.generation-1]
      
 

    # plot fitness
    #print(list(range(self.generation+1)), self.fit_bygeneration)

    if plot_it :
      self.plot_sim()
      plt.pause(0.001)

  
  def run_generation_no_events(self, path_to_sim, plot_it=True):
    # To check variability of fitness for the same genome
    # EXACTLY ONE event per generation, determined by prob p_inv.
    
    print("\n\n----\n")

    #keep old info in case mutation is BAAAAAD.
    old_gene_info = self.gene_info.copy()
    old_prot = self.prot.copy()
    self.evs.append('no')

    self.write_sim_files(path_to_sim)
    
    # update gen counter
    self.generation += 1
    print("Now at generation", self.generation)

    # run simulation
    sim.start_transcribing('../sim_files/current/params.ini', "../sim_files/future")
    fit = self.fitness("../sim_files/future/save_tr_nbr.csv")

    # keep or not
    # if new fitness is inferior to old fit. 
    self.fit_bygeneration.append(fit) 

    # plot fitness
    #print(list(range(self.generation+1)), self.fit_bygeneration)

    if plot_it :
      self.plot_sim()
      plt.pause(0.001)


  def plot_sim(self):

    gens = list(range(self.generation+1))
    p1, = plt.plot(gens, self.fit_bygeneration, lw=1, c='k')
    cols = ['r', 'g','b']
    mks = ['x', 's', 'D']
    p2 = [None]*3
    for c, e in enumerate(['del', 'ins', 'inv']):
      if e in self.evs:
        idx = [i for i, j in enumerate(self.evs) if j == e]
        g = [gens[i] for i in idx]
        f = [self.fit_bygeneration[i] for i in idx]
        p2[c] = plt.scatter(g, f, marker=mks[c], c=cols[c], s=15)
    plt.legend([p1,p2[0], p2[1], p2[2]],['fitness', 'del', 'ins', 'inv'], loc=2)



  def fitness(self, path_to_sim):
    """returns fitness of individual"""
    
    #takes expression profile stored after simulation at path_to_sim
    future = pd.read_table(path_to_sim, header=None).iloc[:,0]
    future = pd.concat([future, self.gene_info['TUindex']], axis=1)
    future.sort_values(by=['TUindex'])
    future = 1.0*future.iloc[:,0]
    ideal = self.env.iloc[:,1]
    future = future/sum(future)
    fit = np.exp(-sum((future-ideal)/ideal))

    return fit
    













    
# class Gene:
#     """docstring for Genome"""

#     def __init__(self, s,e,ori,prom_str):
#         self.s = s
#         self.e = e
#         self.ori = ori
#         self.prom_str = prom_str