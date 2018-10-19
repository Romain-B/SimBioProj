#/usr/bin/env python
#coding:utf-8

import pandas as pd
import numpy as np


class Genome(object):

	def __init__(self, Size, path_init):
		self.Size = Size
		self.sec_dist = 60 # security distance of prot and genes

		TSS = pd.read_table(path_init+'/TSS.dat', header=0)
		TTS = pd.read_table(path_init+'/TTS.dat', header=0)
		self.gene_info = TSS.merge(TTS)
		self.prot = pd.read_table(path_init+'/prot.dat', header=0)['prot_pos']
		self.env = pd.read_table(path_init+'/environment.dat')

	def __str__(self):
		s = "Genome Info :\n--------------\n"
		s += "Size : "+str(self.Size)+"\n"
		s += "Gene Info : \n"+str(self.gene_info)+"\n\n"
		s += "Barriers : \n"+str(self.prot)+"\n\n"	
		return s



		
		
# class Gene:
# 	"""docstring for Genome"""

# 	def __init__(self, s,e,ori,prom_str):
# 		self.s = s
# 		self.e = e
# 		self.ori = ori
# 		self.prom_str = prom_str