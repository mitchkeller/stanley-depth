#!/usr/bin/env python3
# encoding: utf-8
"""
sdepth_down.py

Created by Mitch Keller on 2016-07-06.
Takes k-uniform hypergraphs (equivalently Cohen-Macaulay simplicial complexes of dimension k-1)
and determines if the Stanley depth of S/I (where I is the associated Stanley-Reisner ideal) is k.
If so, the hypergraph (in g6 format) is written out to the 'good' file. If not, it is written to the
'bad' file. Does complementation so that input is only required up to n/2 and k/2.

Copyright (c) 2014-2016 . All rights reserved.

Command line arguments (in order):
n - size of the universal set
k - size of sets in antichain
e - number of sets in antichain
part - Which batch of (n,k,e) data is this?
"""

from sdepth import *
from ast import literal_eval
import networkx as nx
from itertools import combinations
import sys
import pprint

def graph_to_antichain(g,n,e):
    x = list(map(set,list(map(nx.all_neighbors,[g]*(n+e),g.node))))
    return [s for s in x if (len(s) > 0 and max(s) < n)]

def process_antichains(graphs,antichains,n,k,e,part):
    filename_prefix = "/home/kellermt/partitionability/"+str(n)+"-hypergraphs/results/result-"+str(n) + "-" + str(k) + "-" + str(e) + "-" + part + "-"
    filename_suffix = ".g6"
    print("Processing antichains with n={}, k={}, e={}, and part={}".format(n,k,e,part))
    pprint.pprint(antichains)
    count = len(antichains)
    results = zip(map(sdepth,[[frozenset([])]]*count,antichains,[k]*count,[set(range(n))]*count,[None]*count),graphs)

    with open(filename_prefix+"good"+filename_suffix,'w') as good,\
         open(filename_prefix+"bad"+filename_suffix,'w') as bad:
        for (sd_result,graph) in results:
            g=nx.generate_graph6(graph,header=False)
            if (sd_result):
                good.write(g+"\n")
            else:
                print("Result: {0}, graph = {1}, antichain = {2}".format(sd_result,g,antichains[graphs.index(graph)]))
                bad.write(g+"\n")
    return

def main():
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    e = int(sys.argv[3])
    part = sys.argv[4]

    infile = "/home/kellermt/partitionability/"+str(n)+"-hypergraphs/"+str(n)+"-"+str(k)+"-"+str(e)+"-hypergraphs-"+part+".g6"

    graphs = nx.read_graph6(infile)
    if (e<=1 or e >= choose(n,k)-1):
        graphs = [graphs]
    n_list = [n]*len(graphs)
    e_list = [e]*len(graphs)
    antichains_n_k_e = list(map(graph_to_antichain,graphs,n_list,e_list))
    antichains_n_n_minus_k_e = [[set(range(n)) - x for x in a] for a in antichains_n_k_e]
    k_sets = set(map(frozenset,combinations(range(n),k)))
    antichains_n_k_n_choose_k_minus_e = [list(map(set,k_sets - set(map(frozenset,x)))) for x in antichains_n_k_e]
    antichains_n_n_minus_k_n_choose_k_minus_e = [[set(range(n)) - x for x in a] for a in antichains_n_k_n_choose_k_minus_e]

    process_antichains(graphs,list(map(stanley_reisner_generators,antichains_n_k_e,[set(range(n))]*len(graphs))),n,k,e,part)
    if k != n-k:
        process_antichains(graphs,list(map(stanley_reisner_generators,antichains_n_n_minus_k_e,[set(range(n))]*len(graphs))),n,n-k,e,part)
    if e != choose(n,k) - e:
        process_antichains(graphs,list(map(stanley_reisner_generators,antichains_n_k_n_choose_k_minus_e,[set(range(n))]*len(graphs))),n,k,choose(n,k)-e,part)
    if k != n-k and e != choose(n,k) - e:
        process_antichains(graphs,list(map(stanley_reisner_generators,antichains_n_n_minus_k_n_choose_k_minus_e,[set(range(n))]*len(graphs))),n,n-k,choose(n,k)-e,part)    







if __name__ == '__main__':
    main()
