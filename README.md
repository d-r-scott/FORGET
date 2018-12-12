# grouping
Candidate grouping algorithm for CRAFT Real-Time FRB detection pipeline
Author: David Scott
Date: 2018-12-12

This is an implementation of the clustering algorithm used in Heimdall written by Ben Barsdell and Andrew Jameson.
  (https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu)
  
The main python file (grouping.py) is intended to be run similarly to FREDDA's Friends-of-Friends script written by Keith Bannister.

Usage is as follows:

  ./grouping.py filename.cand
  
Options available include:

  -h, --help: Show the help message
    
  -t TTOL, --ttol TTOL: Set the tolerance for the number of time samples events may be separated by in order to be considered coincident [Default: 3]
    
  -d DMTOL, --dmtol DMTOL: Set the tolerance for DM of candidates, in units of pc cm^-3 [Default: 2.0]
    
  -w WTOL, --wtol WTOL: Set the tolerance for boxcar width of candidates, in units of time samples [Default: 2]
    
  --tmin TMIN: Set the earliest time sample to consider [Default: 0]

  --tmax TMAX: Set the latest time sample to consider [Default: 9999999999] 
    
  --dmmin DMMIN: Set the minimum DM to consider (in pc cm^-3} [Default: 0]
    
  --wmax WMAX: Set the maximum width to consider (in time samples) [Default: 32]
    
  -snmin SNMIN: Set the minimum S/N to consider [Default: 0]

  -p, --plot: Show some plots of the candidates in time/DM/width space. Colour indicates S/N, with more green having higher S/N. Left: Candidates before grouping. Middle: Candidates before grouping with error bars indicating the tolerance ranges. Right: Candidates after grouping.

  -l --latency: Calculate the latency (run-time) of the grouping algorithm. Does not include the time taken for file writing and plot generation.
