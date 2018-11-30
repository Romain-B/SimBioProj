# SimBioProj
Computational biology project


The aim of this project is to generate a genome evolution algorithm without mutation and use [https://github.com/bilalelhoudaigui/TCDS-v2.git] to compute its transcriptional profile (`simulation2.py` file is from the repo). 
This profile is then used to compute a fitness and, on the basis of a Monte Carlo Metropolis algorithm, tend to an optimal profile for the simulated environment.

Simulations plot fitness over time (generations) and the type of modification applied are indicated on the plot.
During the simulation, an output file (specified by `path_output`, `sim_files/gene_expr.csv` by default) is generated and updated at each iteration containing the generation, mutation event, fitness and expression level of each gene in the genome.

To allow for graphical representation of the genome, we generated dummy genebank files and use tools from the `reportlab` package, which you can install easily with `pip3` (using the `--user` option should you lack admin rights). 
The resulting images are in the `plotting/` folder. 
Please note that, as of now, the default setting for `run_generation()` is to plot the genome every 10 generation (which may result in quite a bit of outputs for larger simulations), tweak the  `plot_gen` parameter if needed.

