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
		self.prot = pd.read_table(path_init+'/prot.dat', header=0)
		self.env = pd.read_table(path_init+'/environment.dat')

	def __str__(self):
		s = "Genome Info :\n--------------\n"
		s += "Size : "+str(self.Size)+"\n"
		s += "Gene Info : \n"+str(self.gene_info)+"\n\n"
		s += "Barriers : \n"+str(self.prot)+"\n\n"	
		return s


	def nearest_obj_distance(self, p):
		# Needed to ensure that genes and barriers don't become less that sec_dist apart.
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
		# checks is inversion positions are ok.
		
		# check if pos is in gene
		if self.pos_in_gene(s) or self.pos_in_gene(e) :
			return False 

		# check the distance between barriers and gene
		s_gd, s_bd = self.nearest_obj_distance(s)
		e_gd, e_bd = self.nearest_obj_distance(e)
		
		return True if s_gd+e_bd > self.sec_dist and s_bd+e_gd > self.sec_dist else False

	def pos_in_gene(self, p):
		# Returns True if position is in a gene

		for i, row in self.gene_info.iterrows():
			if p>row['TSS_pos'] and p<row['TTS_pos']:
				return True
		return False


	def get_inv_pos(self):
		# Returns 2 good positions for inversion

		s,e = np.sort(np.random.randint(0, self.Size, size=2))
		i=0
		while not self.good_inv_pos(s,e):
			s,e = np.sort(np.random.randint(0, self.Size, size=2))
			i+=1
		
		return s,e
	




		
		
# class Gene:
# 	"""docstring for Genome"""

# 	def __init__(self, s,e,ori,prom_str):
# 		self.s = s
# 		self.e = e
# 		self.ori = ori
# 		self.prom_str = prom_str