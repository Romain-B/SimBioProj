#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome
import matplotlib.pyplot as plt


#G = Genome(30000, '../sim_files/init', NO_PLOT=True)
G = Genome(30000, '../sim_files/init', indel_var=500, p_inv=0.1, NO_PLOT=True)

best_genome = False

if best_genome :
  [G.run_generation_no_events('../sim_files/current') for i in range(5)]
  [G.run_generation('../sim_files/current') for i in range(5)]

  # to see how look like the best genome we made
  l =['Best genome found\n\n', str(G.best[0]), 
      '\n\nBarriers positions\n\n', str(G.best[1]), 
      '\n\nBest fitness found\n\n', str(G.best_fit)]
  with open("../sim_files/current/best_genome.txt", "w") as output:
    for el in l:
      output.write(el)


  # to keep information about how fitness evolve with events
  with open("../sim_files/current/fitness_by_event.txt", "w") as output:
    output.write(" ".join(map(str,G.evs)))
    output.write("\n")
    output.write(" ".join(map(str,G.fit_bygeneration)))
    output.write("\n")
    
  #just to have a first idea about how fitness evolves with each event
  diff_fit = np.array(G.fit_bygeneration[:-1]) - np.array(G.fit_bygeneration[1:])
  fits = pd.DataFrame({'events':G.evs[:-1], 'fit':diff_fit})
  print("\n----\n", fits.groupby('events').mean())


  plt.show(block=True)








else :
  #[G.run_generation_no_events('../sim_files/current') for i in range(50)]
  [G.run_generation('../sim_files/current') for i in range(100)]
  #plt.savefig('../sim_files/fitness_by_gen.png')
  plt.show(block=True)




