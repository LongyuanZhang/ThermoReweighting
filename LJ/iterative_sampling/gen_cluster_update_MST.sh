#!/bin/bash
set -e
generation=100
LenPerRun=200
AT=0.15
KAPPA=885375
R_0=0.15
N_LJ=10

WDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
MDP=/home/lzhang657/Schmidt_work/NUCL/thermo_reweight/pub_data/LJ/MD_update_random_MST/mdp
ITP=/home/lzhang657/Schmidt_work/NUCL/thermo_reweight/pub_data/LJ/MD_update_random_MST/itp

cat << EOF > contact.py
import sys
from enum import Enum
from math import isnan
from operator import itemgetter
import copy
from random import shuffle

import networkx as nx
from networkx.utils import UnionFind, not_implemented_for, py_random_state

class EdgePartition(Enum):
    """
    An enum to store the state of an edge partition. The enum is written to the
    edges of a graph before being pasted to "kruskal_mst_edges". Options are:

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

# Read the contact matrix as a graph (This includes ++ and -- contacts too)
G=nx.Graph(nx.drawing.nx_pydot.read_dot('contact.dot'))

print(int(nx.algorithms.components.is_connected(G)))

sys.stdout = open("mst", "w")
if nx.algorithms.components.is_connected(G):
  MST=my_kruskal_mst(G)
  for line in nx.generate_edgelist(MST, data=False):print(line)

sys.stdout.close()
EOF

cp $ITP/template.top ./topol.top
for i in `seq 1 $N_LJ`
do
cat << EOF >> ./topol.top
H                  1
EOF
done

for j in {0..4}
do

var=$((RANDOM%100))
mkdir -p run${j}
cp ../../MD_start_from_MC/N${N_LJ}/gro_gen_MC/gro0/MC${var}.gro run${j}/start.gro
cp topol.top run${j}
cp contact.py run${j}

cd run${j}

echo 'q' | gmx make_ndx -f start.gro &> /dev/null

echo "System: GROUP NDX_FILE=index.ndx NDX_GROUP=System" > mst.dat
echo "WHOLEMOLECULES ENTITY0=System" >> mst.dat
echo "mat: CONTACT_MATRIX ATOMS=System SWITCH={RATIONAL R_0=${R_0} NN=10000}" >> mst.dat
echo "dfs: DFSCLUSTERING MATRIX=mat" >> mst.dat
echo "nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1" >> mst.dat
echo "PRINT ARG=nat FILE=NAT" >> mst.dat
echo "DUMPGRAPH MATRIX=mat FILE=contact.dot" >> mst.dat
plumed driver --igro start.gro --plumed mst.dat
python contact.py
awk '{print $1+1","$2+1}' mst > edges1

# restraints setup via PLUMED
echo "System: GROUP NDX_FILE=index.ndx NDX_GROUP=System" > plumed.dat
echo "WRAPAROUND ATOMS=System AROUND=1" >> plumed.dat

n=1
while read p
do
echo "#D($p)" >> plumed.dat
echo "DISTANCE ATOMS=$p LABEL=d$n" >> plumed.dat
#echo "RESTRAINT ARG=d$n AT=$AT KAPPA=250.0 LABEL=uwall$n" >> plumed.dat
echo "UPPER_WALLS ARG=d$n AT=$AT KAPPA=$KAPPA EXP=2 EPS=1 OFFSET=0 LABEL=uwall$n" >> plumed.dat
echo "PRINT ARG=d$n,uwall${n}.bias STRIDE=500 FILE=uwall_$n" >> plumed.dat
echo " " >> plumed.dat
n=$((n+1))
done < edges1

#em
#cp $MDP/min.mdp .
#gmx grompp -f min.mdp -c start.gro -p topol.top -o min.tpr
#gmx mdrun -ntmpi 1 -ntomp 16 -deffnm min

#nvt
cp $MDP/nvt.mdp .
gmx grompp -f nvt.mdp -c start.gro -p topol.top -o nvt.tpr
gmx mdrun -ntmpi 1 -ntomp 16 -deffnm nvt --plumed plumed.dat

#MD1
cp $MDP/md.mdp .
gmx grompp -f md.mdp -p topol.top -c nvt.gro -t nvt.cpt -o MD1.tpr
gmx mdrun -ntmpi 1 -ntomp 16 -deffnm MD1 --plumed plumed.dat

#
for i in `seq 2 100`
do

echo 'q' | gmx make_ndx -f MD$(($i-1)).gro -o index.ndx

R_0=0.15
connectivity='0'

while [[ "$connectivity" != '1' ]]
do
        echo "System: GROUP NDX_FILE=index.ndx NDX_GROUP=System" > mst.dat
        echo "WHOLEMOLECULES ENTITY0=System" >> mst.dat
        echo "mat: CONTACT_MATRIX ATOMS=System SWITCH={RATIONAL R_0=${R_0} NN=10000}" >> mst.dat
        echo "dfs: DFSCLUSTERING MATRIX=mat" >> mst.dat
        echo "nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1" >> mst.dat
        echo "PRINT ARG=nat FILE=NAT" >> mst.dat
        echo "DUMPGRAPH MATRIX=mat FILE=contact.dot" >> mst.dat
        plumed driver --igro MD$((i-1)).gro --plumed mst.dat
        connectivity=$(python contact.py)
        if [[ "$connectivity" == '1' ]]
        then
                awk '{print $1+1","$2+1}' mst > edges${i}
                break
        fi
        R_0=$(awk -v var="$R_0" 'BEGIN{ printf "%.3f", var + 0.01 }')
done

# restraints setup via PLUMED
echo "System: GROUP NDX_FILE=index.ndx NDX_GROUP=System" > plumed.dat
echo "WRAPAROUND ATOMS=System AROUND=1" >> plumed.dat

n=1
while read p
do
echo "#D($p)" >> plumed.dat
echo "DISTANCE ATOMS=$p LABEL=d$n" >> plumed.dat
#echo "RESTRAINT ARG=d$n AT=$AT KAPPA=250.0 LABEL=uwall$n" >> plumed.dat
echo "UPPER_WALLS ARG=d$n AT=$AT KAPPA=$KAPPA EXP=2 EPS=1 OFFSET=0 LABEL=uwall$n" >> plumed.dat
echo "PRINT ARG=d$n,uwall${n}.bias STRIDE=500 FILE=uwall_$n" >> plumed.dat
echo " " >> plumed.dat
n=$((n+1))
done < edges${i}

gmx grompp -f md.mdp -c MD$(($i-1)).gro -p topol.top -o MD${i}.tpr
gmx mdrun -ntmpi 1 -ntomp 16 -deffnm MD${i} --plumed plumed.dat

done
cd $WDIR

done
