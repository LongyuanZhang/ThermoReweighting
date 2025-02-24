import mdtraj as md
import numpy as np
import networkx as nx
from sympy import *
from collections import Counter

cutoff = 0.15

# load edges
def load_MST_edges(filename):
    G = nx.Graph()
    with open(filename) as f:
        lines = f.readlines()
    for line in lines:
        nodes = line.strip().split(",")
        G.add_edge(int(nodes[0])-1, int(nodes[1])-1)
    return G

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

def check_current_MST(current_graph_edges, traj_index):
    is_current_MST_connected = True
    current_MST = load_MST_edges("../edges%s"%traj_index)
    current_MST_edges = [set(e) for e in current_MST.edges]
    for e in current_MST_edges:
        if e not in current_graph_edges:
            is_current_MST_connected = False
            break

    return is_current_MST_connected

def calc_overrepresented(current_graph):
    num_nodes = current_graph.number_of_nodes()
    N = 0
    L1 = nx.laplacian_matrix(current_graph)
    L2 = L1.todense()
    L2sub = L2[1:, 1:]
    N = np.linalg.det(L2sub)

    return N

def write_list(filename, List):
    with open(filename, "w") as f:
        for element in List:
            f.write("%s\n"%element)

# list of sampled graphs
N_traj_sampled = 100
# list to record N, which is the number of graphs that can be the subgraphs of graph represented by current coordinates
overrepresented_factors = []
PE = np.array([])
for i in range(1, N_traj_sampled+1):
    PE_MST = np.loadtxt("../MD%s.xvg"%i, comments=["@","#"], usecols=(1,))
    PE = np.concatenate((PE, PE_MST))
write_list("PE_100traj.txt", PE)
count_list=[]
for i in range(1, N_traj_sampled+1):
    traj = md.load("../MD%s.xtc"%i, top="../MD1.gro")
    for j in range(traj.n_frames):
        traj_frame = traj[j]
        current_graph = construct_current_graph(traj_frame) 
        current_graph_edges = [set(e) for e in current_graph.edges]
        if check_current_MST(current_graph_edges, i):
            N = calc_overrepresented(current_graph)
            overrepresented_factors.append(N)
            count=True
        else:
            count=False
        count_list.append(count)
PE_counted = PE[count_list]
PE_after_correction = PE_counted / overrepresented_factors
write_list("overrepresented_factor.txt", overrepresented_factors)
write_list("PE_100traj_counted.txt", PE_counted)
write_list("is_frame_counted.txt", count_list)
print(np.mean(overrepresented_factors) * np.mean(PE_after_correction))
write_list("PE_100traj_corrected.txt", PE_after_correction)

max_overlap = np.max(overrepresented_factors)
PE_byN = [[]] * max_overlap
N_byN = [0] * max_overlap

c = Counter(overrepresented_factors)
for i in range(1, max_overlap+1):
    N_byN[i-1] = c[i]

N_byN = np.array(N_byN)
N_reweighted = (np.sum(N_byN)) / (np.sum(N_byN / np.array(list(range(1, max_overlap+1)))))
print(N_reweighted)
print(N_reweighted * np.mean(PE_after_correction))
