import numpy as np

cols = "LigIdx PrtIdx v1 v2 v3 v4 v5 v6"
num_perturbations = 6
tra = 0.2
# -r 0.0872 <=> 5 deg
rot = 0.0872

print(cols)
for ceil in range(0, num_perturbations):
    tra_vec = np.random.uniform(tra * ceil, tra * ceil + tra, 3)
    rot_vec = np.random.uniform(rot * ceil, rot * ceil + rot, 3)
    print("0 0 " + " ".join(map(str, tra_vec)) + " " + " ".join(map(str, rot_vec)))
