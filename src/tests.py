#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome

if False:
  path = '../sim_files/init/TSS.dat'

  TSS = pd.read_table(path, header=0)
  # print(TSS)


  path = '../sim_files/init/TTS.dat'

  TTS = pd.read_table(path, header=0)
  # print(TTS)

  gene_info = TSS.merge(TTS)
  
  

  
      


  

G = Genome(30000, '../sim_files/init')
#print(G)

#print(G.good_inv_pos(2500, 3500))
#print(G.good_inv_pos(2001, 2990))

G.inversion()
G.inversion()