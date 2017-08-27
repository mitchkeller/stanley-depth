#!/bin/bash
#PBS -lnodes=1:ppn=1,mem=1gb,walltime=48:00:00

/home/kellermt/partitionability/sdepth_down.py $n $k $e $part