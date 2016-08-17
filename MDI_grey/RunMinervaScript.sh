#!/bin/csh
#Shell script to submit a set of HdpCluster MCMC runs to the CSC cluster (Minerva)
##
##
module unload msub
module load matlab
qsub -t 1-25  -v initialise=0,dataName=MetabricImage         RunMDI.pbs
qsub -t 1-25  -v initialise=0,dataName=MetabricImagePca      RunMDI.pbs
qsub -t 1-25  -v initialise=0,dataName=MetabricMoments       RunMDI.pbs
qstat
showq -r -u smsgal
