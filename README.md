### MFACT: An MPI Fast Application Classification Tool based on SST DUMPI Trace. 

If you use the resources available here in your work, please refer to 
this paper as the source.

  @INPROCEEDINGS{7516058, \
  author={Z. Tong and S. Pakin and M. Lang and X. Yuan}, \
  booktitle={2016 IEEE International Parallel and Distributed Processing Symposium (IPDPS)}, \
  title={Fast Classification of MPI Applications Using Lamport's Logical Clocks}, \
  year={2016}, \
  pages={618-627}, \
  doi={10.1109/IPDPS.2016.40}, \
  ISSN={1530-2075}, \
  month={May}}

Please refer to the SST-DUMPI-README for copyright information.

This README file is used for Fast Classification Tool 
and it is based on the sst-dumpi library from Sandia National Lab. 
This package is part of the Message Disembogulator Suite (MDS), 
known internally as LA-CC-13-036. It concurrently runs multiple copy
of updated dumpi2ascii to process dumpi traces and to model application 
performance over an arbitrary number of network configurations in one
trace run.

Before installation, be sure to install openmpi or mpich \
Files that have been updated for this tool include:

* ~/sst-dumpi/dumpi/bin/dumpi2ascii.c 
* ~/sst-dumpi/dumpi/bin/dumpi2ascii-callback.c 
* ~/sst-dumpi/dumpi/bin/dumpi2ascii-callback.h 

### To Configure 
  ./bootstrap.sh

### To configure
  ./configure CC=mpicc CXX=mpicxx --enable-libdumpi --prefix=$HOME/$DUMPI_PATH

### To Make and install
  make 
  
  make install

### Edit the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/sst-dumpi/dumpi/bin

### To run fast classification tool requires three parameters
*-n: number of ranks deployed need to the same as the number of dumpi traces (ranks)\
*-x: number of ranks per compute node (required to differentiate intra-node and inter-node communication \
*-p: file prefix of dumpi traces 

For instanace, to run a 64 ranks of NPB-BT with 1 rank per node: 

 mpirun -n 64 dumpi2ascii 1 tests/npb-bt-traces/dumpi-2015.07.22.20.22.28-

As a result, a file named "summary" will be generated to show the performance predictions over 147 pre-defined
network configurations. 

There are three sections in the perfomrance summary:
1. Summary of time counters by percentage and classification results 
2. Communication sensitivity to various network configurations
3. Performance summary of prediction results

### Simulation environment parameters:
In ~/sst-dumpi/dumpi/bin/dumpi2ascii-callback.h: 
* MAXNUMRECORD:  maximum number of MPI records in each trace 
* MAXNUMCONFIGS: maximum number of network configurations 
* MEMORY_BW: memory bandwidth for the target system
* MEMORY_LAT: memory latency for the target system  
* OVERLAP_CONTROL: overlap options, lat-first (default) or bw-first

Computation scaling factor is defined and can be adjusted in ~/sst-dumpi/dumpi/bin/dumpi2ascii.c
* config.comp_factor = 1.0; 


### Notes:
1. List of intercepted MPI Operations are shown in ~/sst-dumpi/dumpi/bin/dumpi2ascii-callback.h \
  This covers the majority of the frequently used MPI functions. 
  
2. User defined datatypes are sometimes not recorded properly in dumpi traces. \
  For more details, check trace_initTypes() in ~/sst-dumpi/dumpi/bin/dumpi2ascii.c

3. If you run larger traces on a single node, please check the file limits in /etc/security/limits.conf and update accordingly. \
 user  hard  nofile  500000 \
 user  soft  nofile  500000 
 
### Predictive Modeling
1. Please see how the predictive model is built in model.R.
