#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from scipy.stats import norm
from Quantitative import *
import argparse
from argparse import RawTextHelpFormatter 
  
  
class TrajectoriesForSNP:
  def __init__(self,snpinfo,trajcount):
    trajectories=[]
    for i in range(0,trajcount):
      trajectories.append([])
    self.__traj=trajectories
    self.__snpinfo=snpinfo
    self.__trajcount=trajcount

  def appendFreq(self,repeatnumber,frequency):
    self.__traj[repeatnumber].append(frequency)
  
  def generations(self):
    return len(self.__traj[0])
  
  def snpinfo(self):
    return self.__snpinfo
  
  def repeats(self):
    return self.__trajcount
  
  def trajectory(self,i):
    return self.__traj[i]

#----------

def print_phenotypes(rep,phenotypes,snapshots,f_out):

  snapshots=snapshots|set([1])

  ngen=len(phenotypes)
  nind=len(phenotypes[0])

  for ng in snapshots:#range(1,ngen+1):
    for ni in range(1,nind+1):
      f_out.write("\t".join(map(str,[rep,ng,phenotypes[ng-1][ni-1]]))+"\n")



#----------------------------------------


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Simulator of quantitative trait evolution to a new fitness optimum using forward simulations in diploids assuming random mating and free recombination between loci. 
Details can be found in the associated publication. In brief, a Gaussian fitness function can be specified for the phenotypic range between zero and one by the --selection argument, where the selection strength is defined by the "fitmin" and "fitmax" value, the new optimum by the "mean" and the width of the fitness function by the standard deviation "stdDev".

Two different types of output files are provided:

The outfile-freq is specifying the frequency at each locus (snp_id) for the specified generations (--snapshots) and the specified replicates (--repeat-simulations) in the following format, e.g.:
# snp_id	generation	frequency	replicate
1	0	0.05	1
2	0	0.05	1
3	0	0.05	1
4	0	0.05	1
5	0	0.05	1
1	10	0.1925	1
2	10	0.21	1
3	10	0.3575	1
4	10	0.2625	1
5	10	0.24	1
...

The outfile-pheno is specifying the phenotype values (in the range from zero t0 one) for each individual in the population (number of individuals specified by --Ne) for each genertion (--snapshots) and replicates (--repeat-simulations) the following format, e.g.:
# replicate	generation	phenotype
1	1	0.0
1	1	0.0
1	1	0.0
1	1	0.2
1	1	0.1
1	1	0.0
...

""")
parser.add_argument('--Ne', required=True, dest='ne', type=int, help="The number of diploid individuals in the population.")
parser.add_argument('--selection', required=True, dest='selection', type=str, help="Specifying the selection in the form fitmin:fitmax:mean:stdDev that descripbe the Gaussian fitness function. Note: the mean should be within the interval from zeo to one. See corresponding publication for futther details.")
parser.add_argument('--variant-contributions', required=True, dest='varcontri', type=str, help="A comma separated list specifiying the contribution of every variant to the phenotype. Note: the number of variant contributions has to be identicatical to the number of starting frequencies provided under the -p option.")
parser.add_argument('-p', required=True, dest='p', type=str, help="A comma separated list of starting allele frequencies. Note: the number of starting frequencies has to be identicatical to the number of variant contributions provided under the --variant-contributions option.")
parser.add_argument('--repeat-simulations', required=True, dest='repsim', type=int, help="The number of independent replicate populations to be simulated.")
parser.add_argument('--snapshots', required=True, dest='snapshots', type=str, help="The number of diploid individuals in the population.")
parser.add_argument('--outfile-freq', required=True, dest='outfreq', type=str, help="Name of the output file for frequency trajectory information.")
parser.add_argument('--outfile-pheno', required=True, dest='outpheno', type=str, help="Name of the output file for the phenotype information.")
options = parser.parse_args()


repsim = int(options.repsim)
twone = int(options.ne) * 2
p = map(float,options.p.split(","))                       # start allele frequencies
selection = map(float,options.selection.split(":"))        # selection coefficient, mean, standard deviation
phenocontri = map(float, options.varcontri.split(","))               # contribution of the variants to the phenotype
assert(len(phenocontri) == len(p))
snpcount=len(phenocontri)

startc = [int(twone*i) for i in p]
snapshots=set(map(int, options.snapshots.split(","))) 
maxgen = max(snapshots)

# initialize phenotypes


fitnesCalc=FitnessCalculator(selection[0],selection[1],selection[2],selection[3])

# initialize the trajectory array; I fucking hate 3D arrays...

# simulate and print frequency file
f_freq=open(options.outfreq,"w")
f_freq.write("# snp_id\tgeneration\tfrequency\treplicate\n")

f_out=open(options.outpheno,"w")
f_out.write("# replicate\tgeneration\tphenotype\n")
for rep in range(0,repsim):
  # initialize
  phenotypes=[]
  pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
  for i in range(0,snpcount): # 1st generation
    freq=pop.get_frequencyAt(i)
    f_freq.write("\t".join(map(str,[i+1,0,freq,rep+1]))+"\n")
  phenotypes.append(pop.get_relative_phenotypes())

  for g in range(0, maxgen): # remaining generations
    pop=pop.getNextGeneration(twone,fitnesCalc,phenocontri)
    phenotypes.append(pop.get_relative_phenotypes())
    if (g+1) in snapshots: 
      for i in range(0,snpcount):
        freq=pop.get_frequencyAt(i)
        f_freq.write("\t".join(map(str,[i+1,g+1,freq,rep+1]))+"\n")
  print_phenotypes(rep+1,phenotypes,snapshots,f_out)

f_out.close()        
f_freq.close()


# output phenotypes to file


sys.exit()
########
#phenotypes=[]
#for i in range(0,repsim):
#        phenotypes.append([])
#
#for rep in range(0,repsim):
#        # initialize
#        pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
#        for i in range(0,snpcount):
#                freq=pop.get_frequencyAt(i)
#                trajectories[i].appendFreq(rep,freq)
#        phenotypes[rep].append(pop.get_relative_phenotypes())
#
#        for g in range(0, maxgen):
#                pop=pop.getNextGeneration(twone,fitnesCalc,phenocontri)
#                phenotypes[rep].append(pop.get_relative_phenotypes())
#                for i in range(0,snpcount):
#                        freq=pop.get_frequencyAt(i)
#                        trajectories[i].appendFreq(rep,freq)







