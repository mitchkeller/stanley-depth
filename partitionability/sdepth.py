#!/usr/bin/env python3
# encoding: utf-8
"""
sdepth.py

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
import pprint
from itertools import combinations

def sdepth(I,J,k,max_set=None,witness = None):
    # Determine whether the Stanley depth of I/J is at least k.
    # I and J are assumed to be provided as a lists of frozensets
    # If max_set is not provided then max_set is the union of all the elements of I and J
    #print(I,J,k)

    if max_set == None:
        max_set = set([]).union(*I).union(*J)
    min_set_size = min([len(i) for i in I])
    remains = {r : set([frozenset(c) for c in combinations(max_set,r) if any([frozenset(c) >= i for i in I]) and not any([frozenset(c) >= j for j in J])]) for r in range(min_set_size,k+1)}
    #pprint.pprint(remains)
    if not combinatorial_criterion(remains):
        return False
    if any((not combinatorial_criterion(remains,[e]))  for r in remains for e in remains[r]):
        return False
    
    if witness == None:
        return sdepth_internal([],remains,max_set,k)
    else:
        #pprint.pprint(remains)
        return sdepth_internal(witness,remains,max_set,k)

def find_sdepth(I,J,max_set = None,witness=None):
    # Return the value of the Stanley depth of I/J
    # I and J are assumed to be provided as lists of frozensets
    # If max_set is not provided then max_set is the union of all the elements of I and J
    if max_set == None:
        max_set = set([]).union(*I).union(*J)
    k = min([len(a) for a in I])
  
    last_witness = [(i,i) for i in I if len(i) == k]
    while True:
        #print('Stanley Depth at least ' + repr(k))
        current_witness = []
        if not sdepth(I,J,k+1,max_set,current_witness):
            if not witness == None:
                for w in last_witness:
                    witness.append(w)
            break
        else:
            k += 1
            last_witness = current_witness
    return k


# Note that remains are destroyed in this implimentation
def sdepth_internal(fixed_intervals, remains, max_set, k):
    # Internal function controlling the iteration of backtracking
    # fixed_intervals is the set of intervals that have been decided in the backtraking code, as a side effect gives the witness
    # remains is the remainder of the upset given those fixed intervals
    # remains is stored as ranked dictionary of frozensets
    #pprint.pprint(fixed_intervals)
    if len(remains) == 0:
        return True
#    if len(remains[k]) == 2:
#        pprint.pprint(fixed_intervals)
    minimum_rank = min(remains.keys())
    if len(remains[minimum_rank]) == 0:
        del remains[minimum_rank]
        recursive_result = sdepth_internal(fixed_intervals, remains, max_set,k)
        if recursive_result:
            return True
        else:
            remains[minimum_rank] = set([])
        return False

    # All that is left are rank k elements, add them to fixed intervals and return
    if minimum_rank == k:
        fixed_intervals.extend([(v,v) for v in remains[k]])
        return True

    # Check if all that remains is to do a bipartite matching
    if minimum_rank == k-1:
        graph = { v : [u for u in remains[k] if v <= u] for v in remains[k-1]}
        matching = bipartiteMatch(graph)
        if (len(matching[0]) < len(remains[k-1])): # matching wasn't complete
            return False
        else:
            # MTK also changed this. Commented out line 91 and indented lines 92 and 93
            #fixed_intervals.extend([interval for interval in fixed_intervals]) # add the fixed intervals
            fixed_intervals.extend([(matching[0][m],m) for m in matching[0]]) # add the matching intervals
            fixed_intervals.extend([(b,b) for b in remains[k] - set(matching[0])]) # add the unused maximal elements
            return True
   
    cube_size = k - minimum_rank
    # Create the containment graph between minimals and rank k items
    # Remove any edges that intersect with fixed intervals
    graph = {lower : [(frozenset(upper) | lower) for upper in combinations(max_set-lower,cube_size) if (frozenset(upper) | lower) in remains[k] and not any([interval_intersection(f,(lower,(frozenset(upper)|lower))) for f in fixed_intervals])] for lower in remains[minimum_rank]}
    # Good candidates is an iterator over non-intersecting perfect matchings in the graph

    num_mins = len(graph)
    for candidate in good_candidates(graph):
        #pprint.pprint(candidate)
        fixed_intervals.extend(candidate)
        new_covered = build_intervals(candidate)
        #pprint.pprint(new_covered)
        for j in range(minimum_rank,k+1):
            remains[j] -= new_covered[j]
        recursive_result = sdepth_internal(fixed_intervals,remains,max_set,k)
        if recursive_result:
            return True
        else:
            # Since the recursion failed we need to remove the candidates from the fixed list and restore the intervals we removed
            #pprint.pprint(fixed_intervals)
            # This is where MTK had to make a change, since something about the way the recursion was being done was not correct.
            #fixed_intervals = fixed_intervals[:-num_mins]
            #fixed_intervals = [x for x in fixed_intervals if x[0] not in new_covered[min(new_covered.keys())]]
            #fixed_intervals = [x for x in fixed_intervals if x[0] not in new_covered[minimum_rank]]
            del fixed_intervals[-num_mins:]
            #print("Deleted from fixed_intervals")
            #pprint.pprint(new_covered[min(new_covered.keys())])
            #pprint.pprint([x for x in fixed_intervals if x[0] not in new_covered[min(new_covered.keys())]])

            for j in range(minimum_rank,k+1):
                remains[j] |= new_covered[j]
    return False

# Build a ranked dictionary that contains all the sets covered by the given intervals
def build_intervals(intervals):
    lower_sizes = [len(interval[0]) for interval in intervals]
    upper_sizes = [len(interval[1]) for interval in intervals]
    delta_sets = [interval[1] - interval[0] for interval in intervals]
    num_intervals = len(intervals)
    covered = {j : set([(intervals[i][0] | frozenset(c)) for i in range(num_intervals) for c in combinations(delta_sets[i],j-lower_sizes[i]) if lower_sizes[i] <= j <= upper_sizes[i]]) for j in range(min(lower_sizes),max(upper_sizes)+1)}
    return covered

# Return lists of candidate intervals that do not intersect
# graph is a bipartite graph edge list with minimals and allowable maxes
def good_candidates(graph):
    if len(graph) > len( bipartiteMatch(graph)[0]):
        return []

    minimals = [x for x in graph]

    
    # Now set up the rest of the data structure
    # By "sorting" minimals we can change the order of generation
    pairing_length = [len(graph[x]) for x in minimals]
    current_state = [-1 for x in minimals]
    idx = 0
    maximal_size = len(minimals)

    while True:
        if current_state[idx] == -1:
            current_state[idx] = 0
        else:
            # Back up the current state until we reach a valid change
            #print("Backtracking")
            #pprint.pprint(graph)
            while current_state[idx] == pairing_length[idx] - 1:
                current_state[idx] = -1
                idx -= 1
            # If the rollback takes us all the way back, we have generated everything 
            if idx < 0:
                break
            else:
                current_state[idx] += 1
        # Check whether this is a good state
        bottom = minimals[idx]
        top = graph[bottom][current_state[idx]]

        if not any((minimals[j] <= top) and (bottom <= graph[minimals[j]][current_state[j]]) for j in range(idx)):
            if idx == maximal_size - 1:
                yield [(minimals[i],graph[minimals[i]][current_state[i]]) for i in range(maximal_size)]
            else:
                idx += 1


def interval_intersection(A,B):
    # Returns True if the interval (A[0],A[1]) and the interval (B[0],B[1]) intersect
    # return (A[0] | B[0]) <= (A[1] & B[1])
    return A[0] <= B[1] and B[0] <= A[1]
    
def combinatorial_criterion(poset,antichain= None):
    # Tests whether the upset of the poset from the antichain will have enough space to go to the maximal elements
    # Changed the default argument to None because you never ever ever want a mutable as a default argument valut in Python.
    if antichain == None:
        antichain = [frozenset([])]
    k = min(poset.keys())
    max_height = max(poset.keys())
    counts = { r : len([x for x in poset[r] if any([x >= A for A in antichain])]) for r in poset}
    # The control on this while loop was <, which meant that you never checked to see if you ran out of k-sets for the final matching.
    while k <= max_height:
        if counts[k] == 0:
            del counts[k]
        elif counts[k] < 0:
            return False
        else:
            counts = { r : counts[r] - counts[k]*choose(max_height-k,r-k) for r in counts}
            del counts[k]
        k += 1
    return True
    
def stanley_reisner_generators(antichain,max_set = None):
    '''Given the facets of a simplicial complex (antichain of maximal elements 
    of a downset) find the generators of the associated Stanley-Reisner ideal (
    minimal elements of the complementary upset)'''
    if max_set == None:
        max_set = set().union(*antichain)
    maximum_size = max(len(A) for A in antichain)+1
    generators = []
    for r in range(maximum_size+1):
        generators.extend([frozenset(c) for c in combinations(max_set,r) if not any(frozenset(c) >= g for g in generators) and not any(frozenset(c) <= A for A in antichain)])
    return generators

    
# Hopcroft-Karp bipartite max-cardinality matching and max independent set
# David Eppstein, UC Irvine, 27 Apr 2002
# http://code.activestate.com/recipes/123641-hopcroft-karp-bipartite-matching/

def bipartiteMatch(graph):
	'''Find maximum cardinality matching of a bipartite graph (U,V,E).
	The input format is a dictionary mapping members of U to a list
	of their neighbors in V.  The output is a triple (M,A,B) where M is a
	dictionary mapping members of V to their matches in U, A is the part
	of the maximum independent set in U, and B is the part of the MIS in V.
	The same object may occur in both U and V, and is treated as two
	distinct vertices if this happens.'''
	
	# initialize greedy matching (redundant, but faster than full search)
	matching = {}
	for u in graph:
		for v in graph[u]:
			if v not in matching:
				matching[v] = u
				break
	
	while 1:
		# structure residual graph into layers
		# pred[u] gives the neighbor in the previous layer for u in U
		# preds[v] gives a list of neighbors in the previous layer for v in V
		# unmatched gives a list of unmatched vertices in final layer of V,
		# and is also used as a flag value for pred[u] when u is in the first layer
		preds = {}
		unmatched = []
		pred = dict([(u,unmatched) for u in graph])
		for v in matching:
			del pred[matching[v]]
		layer = list(pred)
		
		# repeatedly extend layering structure by another pair of layers
		while layer and not unmatched:
			newLayer = {}
			for u in layer:
				for v in graph[u]:
					if v not in preds:
						newLayer.setdefault(v,[]).append(u)
			layer = []
			for v in newLayer:
				preds[v] = newLayer[v]
				if v in matching:
					layer.append(matching[v])
					pred[matching[v]] = v
				else:
					unmatched.append(v)
		
		# did we finish layering without finding any alternating paths?
		if not unmatched:
			unlayered = {}
			for u in graph:
				for v in graph[u]:
					if v not in preds:
						unlayered[v] = None
			return (matching,list(pred),list(unlayered))

		# recursively search backward through layers to find alternating paths
		# recursion returns true if found path, false otherwise
		def recurse(v):
			if v in preds:
				L = preds[v]
				del preds[v]
				for u in L:
					if u in pred:
						pu = pred[u]
						del pred[u]
						if pu is unmatched or recurse(pu):
							matching[v] = u
							return 1
			return 0

		for v in unmatched: recurse(v)


# From http://stackoverflow.com/questions/3025162/statistics-combinations-in-python
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

    
	      
