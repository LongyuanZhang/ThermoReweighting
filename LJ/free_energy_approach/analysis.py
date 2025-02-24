import numpy as np
from collections import Counter

kT = 2.479

def write_list(filename, List):
    with open(filename, "w") as f:
        for element in List:
            f.write("%s\n"%element)

# calculate of \tilde{P}_i
FE_kJPerMol = np.array([])
for i in range(250):
    FE = np.loadtxt("FE/TI_graph_%s.txt"%i)
    FE_kJPerMol = np.append(FE_kJPerMol, np.mean(FE))
FE_kT = FE_kJPerMol / kT
exp_FE = np.exp(-1 * FE_kT)
P_tilde_i = exp_FE / np.sum(exp_FE)
write_list("P_tilde_i.txt", P_tilde_i)

PE = np.array([])
for i in range(250):
    PE_MST = np.loadtxt("MD/MD%s.xvg"%i, comments=["@","#"], usecols=(1,))
    PE = np.concatenate((PE, PE_MST))
write_list("PE_250traj.txt", PE)
print(np.mean(PE))

P_tilde_i_times_PE_corrected = 0
PE_corrected_all = []
PE_counted_all = []
for i in range(250):
    PE_corrected_i = np.loadtxt("./MD/evaluation_alltree/PE_corrected_%s.txt"%i)
    PE_counted_i = np.loadtxt("./MD/evaluation_alltree/PE_counted_%s.txt"%i)
    PE_corrected_all.append(np.mean(PE_corrected_i))
    PE_counted_all.append(np.mean(PE_counted_i))
P_tilde_i_times_PE_corrected = np.sum(P_tilde_i * np.array(PE_corrected_all))
print(P_tilde_i_times_PE_corrected)

write_list("PE_corrected_bytraj.txt", PE_corrected_all)
write_list("PE_counted_bytraj.txt", PE_counted_all)

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
