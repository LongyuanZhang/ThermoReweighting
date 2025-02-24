import numpy as np
from collections import Counter

kT = 2.479

def write_list(filename, List):
    with open(filename, "w") as f:
        for element in List:
            f.write("%s\n"%element)

P_tilde_i = np.loadtxt("P_tilde_i.txt")

PE = np.loadtxt("PE_250traj.txt")

PE_corrected_all = np.loadtxt("PE_corrected_bytraj.txt")
PE_counted_all = np.loadtxt("PE_counted_bytraj.txt")
P_tilde_i_times_PE_corrected = np.sum(P_tilde_i * np.array(PE_corrected_all))

overrepresented_factors = np.array([], dtype="int")
for i in range(1,11):
    factors_i = np.loadtxt("./MD/evaluation_alltree/overrepresented_factor_%s.txt"%i, dtype="int")
    overrepresented_factors = np.concatenate([overrepresented_factors, factors_i])

max_overlap = np.max(overrepresented_factors)
N_byN = [0] * max_overlap

c = Counter(overrepresented_factors)
for i in range(1, max_overlap+1):
    N_byN[i-1] = c[i]

N_byN = np.array(N_byN)
N_reweighted = (np.sum(N_byN)) / (np.sum(N_byN / np.array(list(range(1, max_overlap+1)))))
print(N_reweighted)
print(N_reweighted * P_tilde_i_times_PE_corrected)
