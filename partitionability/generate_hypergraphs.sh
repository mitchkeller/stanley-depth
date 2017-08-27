#!/bin/bash
#PBS -lnodes=1:ppn=1,mem=200mb,walltime=24:00:00
/home/kellermt/nauty25r9/genbgL -z -d0:$k -D$e:$k $n $e $res/$mod /home/kellermt/partitionability/$n-hypergraphs/$n-$k-$e-hypergraphs-$res.g6
#gzip -9 /home/kellermt/sdepth/$n-hypergraphs/$n-$k-$e-hypergraphs-$res.g6