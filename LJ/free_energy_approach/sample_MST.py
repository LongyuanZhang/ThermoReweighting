import sys
from enum import Enum
from math import isnan
from operator import itemgetter
import copy
from random import shuffle

import networkx as nx
from networkx.utils import UnionFind, not_implemented_for, py_random_state

N_ATOMS=10

class EdgePartition(Enum):
    """
    An enum to store the state of an edge partition. The enum is written to the
    edges of a graph before being pasted to `kruskal_mst_edges`. Options are:

    - EdgePartition.OPEN
    - EdgePartition.INCLUDED
    - EdgePartition.EXCLUDED
    """

    OPEN = 0
    INCLUDED = 1
    EXCLUDED = 2

def my_kruskal_mst_edges(
    G, minimum, weight="weight", keys=True, data=True, ignore_nan=False, partition=None
):
    subtrees = UnionFind()
    edges = G.edges(data=True)

    """
    Sort the edges of the graph with respect to the partition data. 
    Edges are returned in the following order:

    * Included edges
    * Open edges from smallest to largest weight
    * Excluded edges
    """
    included_edges = []
    open_edges = []
    for e in edges:
        d = e[-1]
        wt = d.get(weight, 1)
        if isnan(wt):
            if ignore_nan:
                continue
            raise ValueError(f"NaN found as an edge weight. Edge {e}")

        edge = (wt,) + e
        if d.get(partition) == EdgePartition.INCLUDED:
            included_edges.append(edge)
        elif d.get(partition) == EdgePartition.EXCLUDED:
            continue
        else:
            open_edges.append(edge)

    if minimum:
        sorted_open_edges = sorted(open_edges, key=itemgetter(0))
    else:
        sorted_open_edges = sorted(open_edges, key=itemgetter(0), reverse=True)

    # Condense the lists into one
    included_edges.extend(sorted_open_edges)
    random_edges = included_edges
    del open_edges, sorted_open_edges, included_edges
    shuffle(random_edges)
    
    for wt, u, v, d in random_edges:
        if subtrees[u] != subtrees[v]:
            if data:
                yield u, v, d
            else:
                yield u, v
            subtrees.union(u, v)

def my_kruskal_mst(G, weight="weight", algorithm="kruskal", ignore_nan=False):
    edges = my_kruskal_mst_edges(
        G, algorithm, weight, keys=True, data=True, ignore_nan=ignore_nan
    )
    T = G.__class__()  # Same graph class as G
    T.graph.update(G.graph)
    T.add_nodes_from(G.nodes.items())
    T.add_edges_from(edges)
    return T

G=nx.complete_graph(N_ATOMS)

for i in range(250):
    sys.stdout = open("edges%s"%i, "w")
    MST=my_kruskal_mst(G)
    edges_list = [e for e in MST.edges]
    for nodes in edges_list:
        print("%s,%s"%(nodes[0]+1,nodes[1]+1))

sys.stdout.close()
