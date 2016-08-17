#!/bin/csh
#Shell script to submit a set of HdpCluster MCMC runs to the CSC cluster (Francesca)
##
##
module load matlab
qsub -t 1-20 MDISim.pbs
qstat
