#/usr/bin/env python
#coding:utf-8


try :
  import matplotlib
  matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
  from draw_genome import plot_genome 
except ImportError:
  pass 

import pandas as pd
import copy as cp
import numpy as np
import simulation2 as sim
import os
from shutil import copy2
import subprocess
import csv


class Genome(object):
  """
  Class containing genome of target individual.
  Parameters used in simulations and evolution are to be fixed and given as
  input to initialize object.

  Parameters include :
    - Genome size, 
    - sec_dist (Transcription unit, security distance between genome elements)
    - indel_var (Variability of indel size)
    - p_inv (Inversion probability. Note that p_ins = p_del = (1-p_inv)/2)
    - T0 (Scaling of probability to accept negative fitness change)
    - nb_pol (number of ARN-pol used in transcription simulation, must be < to nb of genes)
  """

  def __init__(self, Size, path_init, path_output="../sim_files/output", NO_PLOT=False,
              sec_dist=60, indel_var=0, p_inv=.1, T0=0.001, nb_pol = 7):
							
    self.Size = Size
    self.sec_dist = sec_dist 
    self.indel_var = indel_var 
    self.p_inv = p_inv 
    self.change_RNAPS(nb_pol, sec_dist) # change number of PRNA and the security distance in params.ini
    self.generation = 0 # generation counter
    self.evs = [None] # events occurred by generation

    TSS = pd.read_table(path_init+'/TSS.dat', header=0)
    TTS = pd.read_table(path_init+'/TTS.dat', header=0)
    self.gene_info = TSS.merge(TTS)
    self.prot = pd.read_table(path_init+'/prot.dat', header=0)
    self.env = pd.read_table(path_init+'/environment.dat', header=None)
    self.path_init = path_init
    self.path_output = path_output
    self.gene_out = self.path_output+'/gene_expr.csv'

    if not os.path.exists(self.path_output):
      os.makedirs(self.path_output)

    self.NO_PLOT = NO_PLOT

    self.fit_bygeneration = []
    self.T0 = T0 

    copy2(self.path_init+'/params.ini','../sim_files/current/params.ini')

    self.write_sim_files("../sim_files/current")

    # Reset the simulation data file if it exists.
    try:
      os.remove(self.gene_out)
    except OSError:
      pass
    open(self.gene_out, 'w').close()
    self.gene_expr_history("fit", 
                           ["g0","g1","g2", "g3", "g4","g5", "g6", "g7", "g8", "g9"], 
                           "event", "keep", g="gen")
    
    ## Run simulation 0
    sim.start_transcribing('../sim_files/current/params.ini', "../sim_files/future")

    fit, fut = self.fitness("../sim_files/future/save_tr_nbr.csv")
    self.fit_bygeneration.append(fit)

    self.gene_expr_history(fit, fut, "NoE", "NA")

    # Keep track of best genome
    self.best = [self.gene_info, self.prot]
    self.best_fit = fit

    if not NO_PLOT:
      # Plot the genome using our draw_genome.py
      self.plot_current_genome()


  def __str__(self):
    """
    Allows for easy printing of genome elements.
    """
    s = "Genome Info :\n--------------\n"
    s += "Size : "+str(self.Size)+"\n"
    s += "Gene Info : \n"+str(self.gene_info)+"\n\n"
    s += "Barriers : \n"+str(self.prot)+"\n\n"  
    return s


  def nearest_obj_distance(self, p):
    """ 
    Takes a position p, returns nearest distance with nearest gene and barrier.
    This is needed to ensure that genes and barriers don't become less that sec_dist apart.
    """

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

    
    ## for extremities of genome cases :
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
    """
    Takes a start and end position (s, e) and checks if they're ok
    for inversion.
    """
    
    # check if pos is in gene
    if self.pos_in_gene(s) or self.pos_in_gene(e) :
      return False 

    # check the distance between barriers and gene
    s_gd, s_bd = self.nearest_obj_distance(s)
    e_gd, e_bd = self.nearest_obj_distance(e)

    if s_gd+e_bd < self.sec_dist or s_bd+e_gd < self.sec_dist :
      return False

    if s_gd+e_gd < self.sec_dist or s_bd+e_bd < self.sec_dist :
      return False
    
    return True

  def pos_in_gene(self, p):
    """
    Returns True if position p is in a gene.
    """

    for i, row in self.gene_info.iterrows():
      if row["TUorient"] == "+":
        if p>=row['TSS_pos'] and p<=row['TTS_pos']:
          return True
      else:
        if p>=row['TTS_pos'] and p<=row['TSS_pos']:
          return True
    return False


  def get_inv_pos(self):
    """
    Returns 2 good positions for inversion.
    """

    s,e = np.sort(np.random.randint(0, self.Size, size=2))
    i=0
    while not self.good_inv_pos(s,e):
      s,e = np.sort(np.random.randint(0, self.Size, size=2))
      i+=1
    
    return s,e
  

  def frag_length(self):
    """
    Returns size of fragment for insertion or deletion 
    (with respect to indel_var).
    """
        
    var = 0 if self.indel_var==0 else np.random.randint(-self.indel_var,self.indel_var)
    return abs(self.sec_dist + self.indel_var + var)


  def insertion(self):
    """
    Insertion event in the genome.
    """
        
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


  def barr_between_pos(self, s, e):
    """
    Returns True if there is a barrier between argument positions.
    This is needed to ensure deletions don't remove barriers.
    """
        
    for i, row in self.prot.iterrows():
      if row['prot_pos'] > s and row['prot_pos'] <= e:
        return True
    return False
  
  def get_barr_between_pos(self, s, e):
    """
    Returns indices of barriers between argument positions.
    This is useful during inversion
    """
        
    barrs = []
    for i, row in self.prot.iterrows():
      if row['prot_pos'] > s and row['prot_pos'] <= e:
        barrs.append(i)
    return barrs

  def deletion(self):
    """
    Deletion event in the genome.
    """
        
    l = self.frag_length()

    # find good del position
    p1 = np.random.randint(0, self.Size)
    p2 = p1 + np.random.choice([-1,1])*l
    s, e = np.sort([p1,p2]) 

    while not self.good_inv_pos(s,e) or self.barr_between_pos(s,e) or s <= 0 or e > self.Size:
      l = self.frag_length()
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

  def find_genes(self,s,e):
    """
    Takes a start and end position (s, e), returns a list of all genes in between.
    """
    
    genes_list = []
    for index,row in self.gene_info.iterrows(): #iterate over all start positions
      start_gene = row["TSS_pos"]
      if start_gene >= s and start_gene <= e:
        genes_list.append(index)
    # Note : no need to check for end of genes because find_genes is used only
    #        for inversions and inversions will never cut a gene
    return genes_list
          
  
  def inversion(self):
    """
    Inversion event in the genome.
    """
    
    start_inv,end_inv = self.get_inv_pos()

    genes_to_inv = self.find_genes(start_inv,end_inv)
    genes_copy = cp.deepcopy(self.gene_info)
    if len(genes_to_inv) :
      for ind in genes_to_inv:
        gene = cp.deepcopy(self.gene_info.iloc[ind,:])
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
    """
    Takes the path to simulation folder and writes the current genome to
    TSS, TTS, prot and gff files.
    """
    self.gene_info.to_csv(path_to_sim+'/TSS.dat', sep="\t", columns=['TUindex', 'TUorient', 'TSS_pos', 'TSS_strength'], index=False)
    self.gene_info.to_csv(path_to_sim+'/TTS.dat', sep="\t", columns=['TUindex', 'TUorient', 'TTS_pos', 'TTS_proba_off'], index=False)
    self.prot.to_csv(path_to_sim+'/prot.dat', sep="\t", index=False)

    with open(self.path_init+'/tousgenesidentiques.gff', 'r') as f:
      gff = f.readlines()

    gff[3] = '##sequence-region\ttousgenesidentiques\t1\t'+str(self.Size)+'\n'
    gff[4] = 'tousgenesidentiques\tRefSeq\tregion\t1\t'+str(self.Size)+'\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n'


    with open(path_to_sim+'/tousgenesidentiques.gff', 'w') as f:
      f.writelines(gff)



  def run_generation(self, path_to_sim, plot_it=True, plot_gen=10):
    """
    Takes the path to simulation folder and plot boolean arguments.
    Runs one generation. One generation implies :
      - Do one of the 3 events according to their probability
      - Write new genome to simulation files
      - Run the simulation 
      - Check output and compute fitness
      - Decide if new genome is kept or not.
    """
  
    
    print("\n\n----\n")

    #keep old info in case mutation is BAAAAAD.
    old_gene_info = cp.deepcopy(self.gene_info)
    old_prot = cp.deepcopy(self.prot)

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
    fit, fut = self.fitness("../sim_files/future/save_tr_nbr.csv")

    # keep or not
    # if new fitness is inferior to old fit. 
    self.fit_bygeneration.append(fit)
    keep = True

    deltaU = self.fit_bygeneration[self.generation-1] - fit

    if deltaU > 0:
      p = np.exp(-deltaU/self.T0)
      # with proba p, keep it. So if rd>p, reset
      if np.random.random() > p :
        self.gene_info = cp.deepcopy(old_gene_info)
        self.prot = cp.deepcopy(old_prot)
        self.fit_bygeneration[self.generation] = self.fit_bygeneration[self.generation-1]
        keep = False
      
    self.gene_expr_history(fit, fut, event, keep)

    # update best fitness  
    if fit > self.best_fit :
      self.best = [self.gene_info, self.prot]
      self.best_fit = fit


    # plot fitness
    if not self.NO_PLOT:
      if plot_it :
        self.plot_sim()
        plt.pause(0.001)

      # plot genome
      if (self.generation % plot_gen) == 0:
        self.plot_current_genome()

  
  def run_generation_no_events(self, path_to_sim, plot_it=True):
    """
    See documentation for run_generation() function.
    This essentially does the same thing, save for the mutation event on the genome.
    Allows to check for noise in transcription profiles.
    """
    print("\n\n----\n")

    self.evs.append('no')
    self.write_sim_files(path_to_sim)
    
    # update gen counter
    self.generation += 1
    print("Now at generation", self.generation)

    # run simulation
    sim.start_transcribing('../sim_files/current/params.ini', "../sim_files/future")
    fit, fut = self.fitness("../sim_files/future/save_tr_nbr.csv")

    # keep or not
    # if new fitness is inferior to old fit. 
    self.fit_bygeneration.append(fit) 
    self.gene_expr_history(fit, fut, 'NoE', 'NA')


    # plot fitness
    if not self.NO_PLOT:
      if plot_it :
        self.plot_sim()
        plt.pause(0.001)


  def gene_expr_history(self, fit, fut, ev, keep, g=None):
    """
    Takes current computed fitness (fit), relative gene expression (fut), last event (ev)
    and wether the mutation is kept or not (keep).
    Writes this info in gene_expr.csv file.
    """
    
    g = self.generation if g is None else g
    line = [g]+[ev]+[fit]+[keep]+list(fut)

    with open(self.gene_out, "a") as fp:
      wr = csv.writer(fp, dialect='excel')
      wr.writerow(line)



  def plot_sim(self):
    """
    Plots the fitness by generation in real time, indicates mutation events on plot.
    """
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
    """
    Takes path to simulation folder.
    Computes fitness of modified genome from simulation output.
    Returns fitness and relative gene expression levels (ordered g0->g9).
    """
    
    #takes expression profile stored after simulation at path_to_sim
    future = pd.read_table(path_to_sim, header=None).iloc[:,0]
    future = pd.concat([future, self.gene_info['TUindex']], axis=1)
    future.sort_values(by=['TUindex'])
    future = 1.0*future.iloc[:,0]
    ideal = self.env.iloc[:,1]
    future = future/sum(future)
    # Old fitness function := np.exp(-sum((future-ideal)/ideal))
    fit = np.exp(-sum(abs(np.log(future/ideal))))

    return fit, future


  def write_gb_file(self, path):
    """
    Takes a path to the plot directory.
    Generates a genbank file from current genome, used by our plot_genome function.
    """

    s5 = '     ' # before feature name
    s7 = '       '
    s9 = '         ' # after barrier
    s12 = '            ' # after gene
    gbString = "LOCUS"+s7+"TOUSGENESIDENTIQUES"+s7+str(self.Size)+" bp    DNA"+s5+"PLN"+s5+"DEC-2018"
    gbString += "\n\nFEATURES"+s12+" Location/Qualifiers"

    for i, row in self.gene_info.iterrows():
      gbString += "\n"+s5+"gene"+s12
      if row['TUorient'] == '+':
        gbString += str(row['TSS_pos'])+'..'+str(row['TTS_pos'])
      else :
        gbString += "complement("+str(row['TTS_pos'])+'..'+str(row['TSS_pos'])+")"

      gbString += "\n"+s5+s12+"    /gene=\"g"+str(row['TUindex'])+"\""

    for i, row in self.prot.iterrows():
      gbString += "\n"+s5+"barrier"+s9+str(row['prot_pos'])+'..'+str(row['prot_pos']+20)

    gbString += "\nORIGIN\n//"
    with open(path, 'w') as f:
      f.write(gbString)




  def plot_current_genome(self, title="genome_plot"):
    """
    Plots current genome.
    """
    path = "../plotting/"+title+"_gen_"+str(self.generation)+'.gb'
    self.write_gb_file(path)
    plot_genome(path)
		


  def change_RNAPS(self, nb_pol, sec_dist):
    """
    Takes a number of RNA-pols and a transcsription unit size.
    Modifies their associated values in param.ini for simulation.
    """
    f = open(r'../sim_files/init/params.ini', 'r')
    lines = f.readlines()
    lines[26] = "RNAPS_NB = "+str(nb_pol)+'\n'
    lines[20] = "DELTA_X = "+str(sec_dist)+'\n'
    f.close()
    f = open(r'../sim_files/init/params.ini', 'w')
    f.writelines(lines)
    f.close()


    
