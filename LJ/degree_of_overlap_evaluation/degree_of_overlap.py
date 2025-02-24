import mdtraj as md
import numpy as np
import networkx as nx
from sympy import *
from collections import Counter

cutoff = 0.15

# construct current graph
def construct_current_graph(traj_frame):
    G = nx.Graph()
    G.add_nodes_from(list(range(traj_frame.n_atoms)))
    cluster = False
    for i in range(traj_frame.n_atoms):
        for j in range(i+1, traj_frame.n_atoms):
            distance = md.compute_distances(traj_frame, [(i, j)])[0]
            if distance < cutoff:
                G.add_edge(i, j)
    cluster = nx.is_connected(G)
    return G

def calc_overrepresented(current_graph):
    num_nodes = current_graph.number_of_nodes()
    N = 0
    L1 = nx.laplacian_matrix(current_graph)
    L2 = L1.todense()
    L2sub = L2[1:, 1:]
    N = np.linalg.det(L2sub)

    return N

gro = md.load("MD1.gro")
current_graph = construct_current_graph(gro) 
N = calc_overrepresented(current_graph)
print("degree of overlap = %.d"%np.round(N))
