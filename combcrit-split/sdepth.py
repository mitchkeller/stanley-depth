#!/usr/bin/env python
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

#import numpy
#import itertools
from functools import reduce
from collections import Counter
from itertools import combinations
from itertools import product
from itertools import chain
from itertools import dropwhile
from copy import copy
from math import ceil
import sys


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



def build_interval(min,max):
    """Returns a set of frozensets that is the interval of subsets of max that contain min"""
    if (not (min <= max)):
        raise ValueError("min is not a subset of max")
    interval = powerset(max-min)
    return {s | min for s in interval}

def make_downset(antichain):
    """Given an antichain in a subset lattice, returns its downset"""
    downset = set(map(frozenset,antichain))
    for s in antichain:
        downset = downset | powerset(s)
    return downset

def count_by_size(set):
    """Returns an array giving the number of subsets in set of each size"""
    return Counter(map(len,set))

def sat_comb_crit(c,k):
    """Given a Counter c storing the number of subsets of each size in a poset, determine if that poset satisfies the combinatorial criterion for size k"""
    for i in range(k+1):
        new_intervals = c[i]
        for j in range(i,k+1):
            c[j]-=(new_intervals * choose(k-i,j-i))
        if (c.most_common()[:-2:-1][0][1] < 0):
            return False
    return True


def sat_strong_comb_crit(poset,k):
    """
    Returns True if the poset satisfies the strong
    combinatorial criterion for size k. Returns False otherwise.
    """
    for s in poset:
        #print("Checking set "+str(s))
        subposet = {x - s for x in poset if s <= x}
        if (not sat_comb_crit(count_by_size(subposet),k-len(s))):
            #print("Fail!")
            return(False)
    return(True)

def splits_on_x(antichain,x):
    """
    Determines if the down set of antichain can be split into two down
    sets, one where every set contains the element x and the other where
    no set contains the element x.
    """
    have_x = { s for s in antichain if x in s }
    no_have_x = { s for s in antichain if x not in s }
    down_no_have_x = make_downset(no_have_x)
    have_x_without_x = { s - {x} for s in have_x }
    down_have_x_without_x = make_downset(have_x_without_x)
    return (down_have_x_without_x <= down_no_have_x)

def find_splits(antichain):
    """
    Returns a dict with key 'good' containing a list of the integers x in
    set.union(*antichain) for which antichain splits on x and both parts
    of the split satisfy the strong combinatorial criterion. The 'bad' key
    contains a list of the integers x in set.union(*antichain) for which
    antichain splits on x but one part of the split fails the strong
    combinatorial criterion.
    """
    domain = set.union(*map(set,antichain))
    k = max(map(len,antichain))
    d=dict(good=[], bad= [])
    for x in domain:
        have_x = { s for s in antichain if x in s }
        no_have_x = { s for s in antichain if x not in s }
        down_no_have_x = make_downset(no_have_x)
        have_x_without_x = { s - {x} for s in have_x }
        down_have_x_without_x = make_downset(have_x_without_x)
        if (down_have_x_without_x <= down_no_have_x):
            if (sat_strong_comb_crit(down_no_have_x,k) and sat_strong_comb_crit(down_have_x_without_x,k-1)):
                d['good'].append(x)
            else:
                d['bad'].append(x)
    return d

def splits(antichain,n,k):
    """
    Modificatin of find_splits that returns a boolean.  Modified to be more efficient for use in check_sdepth
    """
    domain = set(range(n))
    for x in domain:
        have_x = { s for s in antichain if x in s }
        no_have_x = { s for s in antichain if x not in s }
        have_x_without_x = { s - {x} for s in have_x }
        for s in have_x_without_x:
            for t in no_have_x:
                if(s <= t):
                    break
            else:
                break
        else:
            down_no_have_x = make_downset(no_have_x)
            down_have_x_without_x = make_downset(have_x_without_x)
            if (sat_strong_comb_crit(down_no_have_x,k) and sat_strong_comb_crit(down_have_x_without_x,k-1)):
                return(True)
    return(False)

def strong_splits(antichain,n,k):
    """
    Checks whether there is a split that such that all the k-1 elements sets are covered in the cube containing a splitting variable.
    """
    test_size = choose(n-1,k-1)
    domain = set(range(n))
    for x in domain:
        #print(x)
        no_have_x = {s for s in antichain if x not in s}
        #print(no_have_x)
        down_no_have_x = make_downset(no_have_x)
        #print({s for s in down_no_have_x if len(s) == (k-1)})
        if ( test_size == len({s for s in down_no_have_x if len(s) == (k-1)})):
            return(True)
    return(False)

            
    
# Next two functions from http://rosettacode.org/wiki/Power_set#Python
def list_powerset(lst):
    return reduce(lambda result, x: result + [subset + [x] for subset
                                              in result],lst, [[]])
 
def powerset(s):
    return frozenset(map(frozenset, list_powerset(list(s))))

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

def comb_crit(c,k):
    """Given a Counter c storing the number of subsets of each size in a pooset which satisfies the combintorial critereon for level k, return an array of the number of intervals starting at a given level"""
    for i in range(k+1):
        new_intervals = c[i]
        for j in range(i+1,k+1):
            c[j]-=(new_intervals * choose(k-i,j-i))
    return c

def check_upset_slow(antichain):
   """Checks whether the there is a Stanley partition of the complement of the down set of antichain that reaches the k+1 element sets"""
   max_set = set.union(*map(set,antichain))
   k = max(map(len,antichain))
   downset = make_downset(antichain)
   upset = {s for s in (powerset(max_set) - downset)  if len(s) <= k+1} 
   if not(sat_strong_comb_crit(upset,k+1)):
       return False
   intervals = comb_crit(count_by_size(upset),k+1)
   minimals = [];
   # Do not need to go up to rank k+1 because those are all trivial
   for i in range(k+1):
       rank = {s for s in upset if len(s) == i}
       if (( len(rank) >= intervals[i]) and (intervals[i]>0 )):
           minimals.append( combinations(rank,intervals[i]))
   minimal_iter = product(*minimals)
   print("Built Minimal Iterator")
   count = 0
   for bottoms in minimal_iter:
       print("--- Starting Maximal Iterator: " + str(count))
       maximals = []
       flat_bottoms = list(chain.from_iterable(bottoms))
       for x in flat_bottoms:
           maximals.append(combinations(max_set-x, k+1-len(x)))
       maximal_iter = product(*maximals)
       print("--- Built Maximal Iterator: " + str(count))
       count +=  1
       for tops in maximal_iter:
           if not(collison(flat_bottoms,tops)):
               print(flat_bottoms)
               print(tops)
               return True
   return False
   
def collison(bottoms,tops):
    """Checks whether the intervalls [B_i,T_i] intersect, returns True if there is a collision"""
    for i in range(len(bottoms)):
       for j in range(i+1,len(bottoms)):
           min_inter = (bottoms[i] | bottoms[j])
           if ( (min_inter <= (frozenset(tops[i]) | bottoms[i])) and (min_inter <= (frozenset(tops[j]) | bottoms[j])) ):
               return True
    return False

def check_upset(antichain):
    """Checks whether the there is a Stanley partition of the complement of the down set of antichain that reaches the k+1 element sets"""
    witness = []
    max_set = set.union(*map(set,antichain))
    k = max(map(len,antichain))
    downset = make_downset(antichain)
    upset = {s for s in (powerset(max_set) - downset)  if len(s) <= k+1} 
    if not(sat_strong_comb_crit(upset,k+1)):
        return [antichain, False, witness]
    upset = {s for s in upset if len(s) <= k}
    return [antichain, check_upset_internal([],upset,max_set,k,witness), witness]

def check_upset_internal(fixed_intervals, remains, max_set, k,witness):
    """Internal code for check_upset.  fixed_intervals is the intervals that have already beed decided in the backtracking code, remains is the upset remove those intervals, max_set is the maximum set of the cube we are working in, and k is the rank of the antichain (we are trying to get to level k+1 in the upset)"""
    #print(fixed_intervals)
    if len(remains) == 0:
        for interval in fixed_intervals:
            witness.append(interval)
        return True
    minimal_rank = min(map(len,remains))

    # All that remains is to match the k-sets up to the (k+1)-sets. Rather
    # than backtracking, let's do this by bipartite matching
    if minimal_rank == k:
        #print("New code!")
        remains_k = remains # We only have k-sets

        # Get all the (k+1)-sets and then remove those that are already
        # maximal elements of fixed intervals
        remains_k_plus_1 = set(map(frozenset,combinations(max_set,k+1)))
        for interval in fixed_intervals:
            remains_k_plus_1 = remains_k_plus_1 - {interval[1]}

        # Form the required type of graph
        graph = dict()
        for x in remains_k_plus_1:
            graph[x] = list(y for y in remains_k if y <= x)
        matching = bipartiteMatch(graph)
        if len(matching[0]) < len(remains_k): # Matching wasn't complete
            return False
        else:
            for interval in fixed_intervals: # Fixed intervals work
                witness.append(interval)
            for m in matching[0]: # And need the matching
                witness.append([m,matching[0][m]])
            return True

    cube_size = k+1-minimal_rank
    minimals = [s for s in remains if len(s) == minimal_rank]
    maximals = []
    for x in minimals:
        valid_covers = []
        for c in combinations(max_set -x, cube_size):
            if not(interval_intersection(fixed_intervals,x,frozenset(c))):
                valid_covers.append(frozenset(c))
        maximals.append(valid_covers)

    maximal_iterator = product(*maximals)
    #print(minimals)
    #print(num_mins)
    for candidate in good_candidates(maximal_iterator,minimals):
        #print("candidate is " + repr(candidate))
        new_fixed = [];
        new_removed = set();
        for old in fixed_intervals:
            new_fixed.append(old)
        for i,x in enumerate(minimals):
            top = frozenset(x | candidate[i])
            new_fixed.append([x, top])
            new_removed |= build_interval(x,top)
        recursive_result = check_upset_internal(new_fixed, remains - new_removed, max_set, k,witness)
        if recursive_result:
            #print(witness)
            return True
    return False

def interval_intersection(intervals,A,B):
    """Check whether the interval [A, A \cup B] intersects with any of the intervals in intervals"""
    for x in intervals:
        if (( A <= x[1]) and (x[0]-A <= B)):
            return True
    return False





def good_candidates(maximal_iterator,minimals):
    reject_candidate = False
    for candidate in maximal_iterator:
        if (reject_candidate and candidate[bad_i] == bad_cand_i and candidate[bad_j] == bad_cand_j):
            continue
        else:
            reject_candidate = False
            bad_i = None
            bad_j = None
            bad_cand_i = None
            bad_cand_j = None
        for j in range(1,len(minimals)):
            if reject_candidate:
                break
            for i in range(j):
                if (( minimals[i] - minimals[j] <= candidate[j]) and (minimals[j] - minimals[i]) <= candidate[i]):
                    reject_candidate = True
                    bad_i = i
                    bad_j = j
                    bad_cand_i = candidate[i]
                    bad_cand_j = candidate[j]
                    #print("We need to reject this sucker!")
                    break
        if not reject_candidate:
            yield(candidate)


def check_sdepth(antichain,n,k):
    """Given an antichain check whether it satisfies the degree conditions, passes the downset passes the strong combinatorial critereon, if it splits, if the sdepth of the downset is k, and the sdepth of the upset is k+1"""
    #print(antichain,file=sys.stderr)
    max_set = set(range(n))
    # Check whether every element appears in an element of the antichain
    for i in max_set:
        for x in antichain:
            if i in x:
                break
        else:
            return("bad_degree")
    
    # Check whether no elment appears in every element of the antichain
    for i in max_set:
        for x in antichain:
            if not(i in x):
                break
        else:
            return("bad_degree")

    
    # We now want to check whether the downset of the antichain passes the SCC
    downset = make_downset(antichain)
    if not(sat_strong_comb_crit(downset,k)):
        return("fail_scc_down")

    splits_res = False
    thin_layer_res = False
    
    # At this point, we know that the antichain passes the degree test and SCC
    # We'll now look to see if the (k-1)-sets are covered.

    down_level = {s for s in downset if len(s) == k-1}
    if(len(down_level) == choose(n,k-1)):
        thin_layer_res = True


    antichain = set(map(frozenset,antichain))

        
    # Now we check whether the downset splits

    if splits(antichain,n,k):
        #print("\t" + find_splits(antichain)['good'])
        #return("splits" + str(find_splits(antichain)['good']))
        splits_res = True

    if (thin_layer_res and splits_res):
        return("thin_layer_splits")
    elif (thin_layer_res and k < ceil(n/2)):
        return("thin_layer_bot_half")
    elif splits_res:
        return("splits")

    # Verify that the SCC holds on the upset
    #upset = {s for s in (powerset(max_set) - downset)  if len(s) <= k+1}
    #if not(sat_strong_comb_crit(upset,k+1)):
    #    return("fail_scc_up")

    # Determine whether the sdepth of the upset is >= k+1 or not
    upset = {s for s in (powerset(max_set) - downset)  if len(s) <= k}
    #upset = {s for s in upset if len(s) <= k}
    witness = []
    if(check_upset_internal([],upset,max_set,k,witness)):
        if(thin_layer_res):
            #print(witness)
            return("thin_layer_pass_sdepth_up")
        else:
            #print(witness)            
            return("pass_sdepth_up")
    else:
        if(thin_layer_res):
            return("thin_layer_fail_sdepth_up")
        else:
            return("fail_sdepth_up")

def main():
    pass

if __name__ == '__main__':
    main()
