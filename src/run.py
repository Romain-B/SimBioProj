#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome
import matplotlib.pyplot as plt



change_RNAPS_NB(123)

def change_RNAPS_NB(n_poly):
	"""it modifies the number of polymerases in param.ini"""
	f = open(r'../sim_files/init/params.ini', 'r')
	lines = f.readlines()
	lines[26] = "RNAPS_NB = "+str(n_poly)+'\n'   
	f.close()
	f = open(r'../sim_files/init/params.ini', 'w')
	f.writelines(lines)
	f.close()