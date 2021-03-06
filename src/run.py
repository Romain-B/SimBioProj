#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome

# Read the experimental plan matrix
plan = pd.read_csv('../sim_files/mat_plan.csv', header = 0, sep=",")
print(plan)

len_noise = 20 # iterations without mutations to se the noise due to transcription

# Levels of the various factors in the experimental plan :
indel_var=[50, 500] #0
p_inv=[0.05, 0.2] #0.1
T0=[0.001, 0.01] # 0.001
nb_pol = [5, 9] #6
nb_iter = [250, 500] # length of 'real' simulation (with mutation events)

nfails = 0

# Proceed to plan :
for i, row in plan.iterrows():

  G = Genome(30000, '../sim_files/init', indel_var=indel_var[row['indel_var']], 
  	  p_inv=p_inv[row['p_inv']], T0=T0[row['T0']], nb_pol=nb_pol[row['nb_pol']], path_output = '../sim_files/output/'+row['alias'], NO_PLOT = True)
  [G.run_generation_no_events('../sim_files/current') for i in range(len_noise)]

  try :
    [G.run_generation('../sim_files/current') for i in range(nb_iter[row['nb_iter']])]
    
    # Write the best genome to gb (for later plotting)
    G.write_gb_file("../sim_files/output/"+row['alias']+"/best_genome.gb")
    l =['Best genome found\n\n', str(G.best[0]), 
        '\n\nBarriers positions\n\n', str(G.best[1]), 
        '\n\nBest fitness found\n\n', str(G.best_fit)]

    with open("../sim_files/output/"+row["alias"]+'/best_gene.txt', "w") as output:
      for el in l:
        output.write(el)
  
  # keep info about the failing genome
  except :
    # signal that we enter except :
    print('------------ in except ------------')
    # signal despair of the programmers :
    print('NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')
    # keep track of number of fails :
    nfails += 1
    G.write_gb_file("../sim_files/output/"+row['alias']+"/bad_genome.gb")
    l =['Failling genome\n\n', str(G.best[0]), 
      '\n\nBarriers positions\n\n', str(G.best[1]), 
      '\n\nFitness found\n\n', str(G.best_fit)]
    with open("../sim_files/output/"+row['alias']+"/bad_genome.txt", "w") as output:
      for el in l:
        output.write(el)

# friendly reminder of the simulation failure rate.           
print("fail ", nfails)


    





