#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome
import matplotlib.pyplot as plt

plan = pd.read_csv('../sim_files/mat_plan.csv', header = 0, sep=",")
print(plan)

len_noise = 20 # iterations without mutations to se the noise due to transcription
len_simu = 100 # length of real simulation including mutations

sec_dist=[40,80] #60
indel_var=[0,20] #0
p_inv=[0.05, 0.2] #0.1
T0=[0.001, 0.01] # 0.001
nb_pol = [5, 9] #6

#G = Genome(30000, '../sim_files/init')

for i, row in plan.iterrows():
  G = Genome(30000, '../sim_files/init', sec_dist=sec_dist[row['sec_dist']], indel_var=indel_var[row['indel_var']], 
	p_inv=p_inv[row['p_inv']], T0=T0[row['T0']], nb_pol=nb_pol[row['nb_pol']], path_output = '../sim_files/output/'+row['alias'], NO_PLOT = True)

  [G.run_generation_no_events('../sim_files/current') for i in range(len_noise)]
  [G.run_generation('../sim_files/current') for i in range(len_simu)]

  # to see how look like the best genome we made
  l =['Best genome found\n\n', str(G.best[0]), 
      '\n\nBarriers positions\n\n', str(G.best[1]), 
      '\n\nBest fitness found\n\n', str(G.best_fit)]
  with open("../sim_files/output/"+row["alias"]+'/best_gene.txt', "w") as output:
    for el in l:
      output.write(el)


    





