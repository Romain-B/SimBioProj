import pandas as pd
import numpy as np

class Gene:
	"""docstring for Genome"""

	def __init__(self, s,e,ori,prom_str):
		self.s = s
		self.e = e
		self.ori = ori
		self.prom_str = prom_str


class Genome(object):
	"""docstring for Genome
"""
	def __init__(self, Size,path_to_env,path_tss,path_tts,path_prot):
		self.Size = Size
		st_genes = pd.read_table(path_tss)
		#self.genes = 

		self.env = pd.read_table(path_to_env)
	
	def 
		
		
