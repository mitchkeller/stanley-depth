#!/usr/bin/env python3
# encoding: utf-8
"""
test_scc.py


Copyright (c) 2014-2017 Mitchel T. Keller and Stephen J. Young.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


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

def sat_strong_comb_crit_antichain(antichain):
    return (antichain,sat_strong_comb_crit(make_downset(antichain),4))

def process_input_to_sets(s):
    return set(map(frozenset,literal_eval(s)))

def graph_to_antichain(g,n,e):
    x = list(map(set,list(map(nx.all_neighbors,[g]*(n+e),g.node))))
    return [s for s in x if (len(s) > 0 and max(s) < n)]

def process_antichains(graphs,antichains,n,k,e,part):
    filename_prefix = "/home/kellermt/sdepth/"+str(n)+"-hypergraphs/results/result-"+str(n) + "-" + str(k) + "-" + str(e) + "-" + part + "-"
    filename_suffix = ".g6"

    count = len(antichains)
    results = zip(map(check_sdepth,antichains,[n]*count,[k]*count),graphs)

    with open(filename_prefix+"bad_degree"+filename_suffix,'w') as bad_degree,\
         open(filename_prefix+"fail_scc_down"+filename_suffix,'w') as fail_scc_down,\
         open(filename_prefix+"thin_layer_splits"+filename_suffix,'w') as thin_layer_splits,\
         open(filename_prefix+"splits"+filename_suffix,'w') as splits,\
         open(filename_prefix+"thin_layer_pass_sdepth_up"+filename_suffix,'w') as thin_layer_pass_sdepth_up,\
         open(filename_prefix+"thin_layer_bot_half"+filename_suffix,'w') as thin_layer_bot_half,\
         open(filename_prefix+"pass_sdepth_up"+filename_suffix,'w') as pass_sdepth_up,\
         open(filename_prefix+"thin_layer_fail_sdepth_up"+filename_suffix,'w') as thin_layer_fail_sdepth_up,\
         open(filename_prefix+"fail_sdepth_up"+filename_suffix,'w') as fail_sdepth_up:
        for (sd_result,graph) in results:
            g=nx.generate_graph6(graph,header=False)
            if (sd_result == "bad_degree"):
                bad_degree.write(g+"\n")
            elif (sd_result == "fail_scc_down"):
                fail_scc_down.write(g+"\n")
            elif (sd_result == "thin_layer_splits"):
                thin_layer_splits.write(g+"\n")
            elif (sd_result == "thin_layer_bot_half"):
                thin_layer_bot_half.write(g+"\n")
            elif (sd_result == "splits"):
                splits.write(g+"\n")
            elif (sd_result == "thin_layer_pass_sdepth_up"):
                thin_layer_pass_sdepth_up.write(g+"\n")
            elif (sd_result == "pass_sdepth_up"):
                pass_sdepth_up.write(g+"\n")
            elif (sd_result == "thin_layer_fail_sdepth_up"):
                thin_layer_fail_sdepth_up.write(g+"\n")
            elif (sd_result == "faail_sdepth_up"):
                fail_sdepth_up.write(g+"\n")
    return

def main():
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    e = int(sys.argv[3])
    part = sys.argv[4]

    infile = "/home/kellermt/sdepth/"+str(n)+"-hypergraphs/"+str(n)+"-"+str(k)+"-"+str(e)+"-hypergraphs-"+part+".g6"

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
    
    process_antichains(graphs,antichains_n_k_e,n,k,e,part)
    process_antichains(graphs,antichains_n_n_minus_k_e,n,n-k,e,part)
    process_antichains(graphs,antichains_n_k_n_choose_k_minus_e,n,k,choose(n,k)-e,part)
    process_antichains(graphs,antichains_n_n_minus_k_n_choose_k_minus_e,n,n-k,choose(n,k)-e,part)    







if __name__ == '__main__':
    main()
