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

##Testing inversion:
#print(G)
#G.inversion()
#print(G)
#G.inversion()

# print(G.good_inv_pos(2500, 3500))
# print(G.good_inv_pos(2001, 2990))

# print(G.frag_length())
# print(G.frag_length())
# print(G.frag_length())
print("######################################\n")
print("INSERTION\n")

G.insertion()
# print(G)

print("######################################\n")
print("INSERTION\n")

G.insertion()
# print(G)

print("######################################\n")
print("DELETION\n")


G.deletion()
# print(G)

print("######################################\n")
print("DELETION\n")

G.deletion()
G.deletion()
# print(G)


for i  in range(20):
  G.inversion()
  print(G)



G.write_sim_files('../sim_files/current')

