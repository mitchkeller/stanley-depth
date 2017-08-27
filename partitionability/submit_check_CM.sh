#!/bin/bash
#PBS -lnodes=1:ppn=1,mem=1gb,walltime=48:00:00

/home/kellermt/partitionability/check_CM.sage $n $k $e $part