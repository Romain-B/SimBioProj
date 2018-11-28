#/usr/bin/env python
#coding:utf-8

import numpy as np
import pandas as pd
from Genome import Genome
import matplotlib.pyplot as plt


G = Genome(30000, '../sim_files/init')

[G.run_generation_no_events('../sim_files/current') for i in range(100)]
[G.run_generation('../sim_files/current') for i in range(100)]
plt.show(block=True)
