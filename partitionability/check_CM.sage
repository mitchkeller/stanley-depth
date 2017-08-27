#!/usr/bin/env sage

"""
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
"""

import networkx as nx
import sys

n = int(sys.argv[1])
k = int(sys.argv[2])
e = int(sys.argv[3])
part = sys.argv[4]

infile = "/home/kellermt/partitionability/{0}-hypergraphs/results/result-{0}-{1}-{2}-{3}-bad.g6".format(n,k,e,part)
outfile = "/home/kellermt/partitionability/{0}-hypergraphs/results/result-{0}-{1}-{2}-{3}-bad-CM.g6".format(n,k,e,part)
graphs = nx.read_graph6(infile)

with open(infile) as f:
    i = 0
    for l in f:
        if i > 1:
            break
        i += 1

if (i <= 1):
    graphs = [graphs]
    
if e > binomial(n,k)/2:
    num_hypergraph_edges = binomial(n,k) - e
else:
    num_hypergraph_edges = e

not_CM = 0

with open(outfile,'w') as out:
    for g in graphs:
        temp = list(map(set,list(map(nx.all_neighbors,[g]*(n+num_hypergraph_edges),g.node))))
        hypergraph_antichain = [s for s in temp if (len(s) > 0 and max(s) < n)]
        if e != num_hypergraph_edges:
            if k != len(hypergraph_antichain[0]):
                antichain = set(map(frozenset,Subsets(range(n),k))) - set(map(frozenset,[set(range(n)) - x for x in hypergraph_antichain]))
                #print(set(map(frozenset,Subsets(range(n),n-k))))
                #print(set(map(frozenset,[set(range(n)) - x for x in hypergraph_antichain])))
                #print(antichain)
                #[set(range(n)) - x for x in set(map(frozenset,Subsets(range(n),k))) - set(map(frozenset,hypergraph_antichain))]
            else:
                antichain = set(map(frozenset,Subsets(range(n),k))) - set(map(frozenset,hypergraph_antichain))
        else:
            if k != len(hypergraph_antichain[0]):
                antichain = [set(range(n)) - x for x in hypergraph_antichain]
            else:
                antichain = hypergraph_antichain

        c = SimplicialComplex(antichain)

        if c.is_cohen_macaulay(1):
            #print(antichain)
            out.write(nx.generate_graph6(g,header=False)+"\n")
        else:
            not_CM += 1

print("Processed {} antichains that were not Cohen-Macaulay.".format(not_CM))
    
                            
