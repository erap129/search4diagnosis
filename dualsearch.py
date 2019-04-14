import logiccircuit as lc
import random
from copy import deepcopy
from itertools import product

in1 = lc.Variable()
in2 = lc.Variable()
in3 = lc.Variable()
in4 = lc.Variable()
in5 = lc.Variable()
g10 = lc.Nand(in1, in3)
g11 = lc.Nand(in3, in4)
g16 = lc.Nand(in2, g11)
g19 = lc.Nand(in5, g11)
o1 = lc.Nand(g10, g16)
o2 = lc.Nand(g16, g19)

bool = lc.partial_eval([g10,g16,o1], input_dict={in1: True, in2: True, in3: True, in4: True, in5:True}, return_result=True)

def get_hitting_set(Lambda):
    all_sets =  list(product(*Lambda))
    for i in range(len(all_sets)):
        all_sets[i] = list(set(all_sets[i]))
    return all_sets

def PDDS(DP, Omega, Lambda, state):
    K = 1000
    for i in range(K):
        hitting_sets = get_hitting_set(Lambda)
        found = False
        for hitting_set in hitting_sets:
            if hitting_set not in Omega:
                if state == 'diagnosis':
                    if consistency_diagnosis(DP, hitting_set):
                        Omega.append(hitting_set)
                    else:
                        complement = [item for item in DP.elem if item not in hitting_set]


                found = True
                break
        if not found:
            return True

Lambda = [[g10, g11, g16, o1, o2]]
Omega = []



