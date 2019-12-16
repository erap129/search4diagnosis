import itertools
import os
import re

from pysat.solvers import Glucose3, Minisat22, Lingeling, MapleCM, Minicard, Cadical

import logiccircuit as lc
import random
from copy import deepcopy
from itertools import product, chain, combinations
import pycosat
import numpy as np
import threading
import ctypes
import time
import pandas as pd
import time
import pysat
timer = True
comp_counter = 1
obs_folder = 'data/observations'
sys_folder = 'data/system_descriptions'
Lambda = None
Omega = None

def generate_nand_clause(clause, output, input1, input2):
    return [[-clause, -input1, -input2, -output], [-clause, input1, output], [-clause, input2, output]]

def generate_and_clause(clause, output, input1, input2):
    return [[-clause, -input1, -input2, output], [-clause, input1, -output], [-clause, input2, -output]]

def generate_or_clause(clause, output, input1, input2):
    return [[-clause, input1, input2, -output], [-clause, -input1, output], [-clause, -input2, output]]

def generate_nor_clause(clause, output, input1, input2):
    return [[-clause, input1, input2, output], [-clause, -input1, -output], [-clause, -input2, -output]]

def generate_not_clause(clause, output, input1):
    return [[-clause, -input1, -output], [-clause, input1, output]]

def generate_xor_clause(clause, output, input1, input2):
    return [[-clause, -input1, -input2, -output], [-clause, input1, input2, -output], [-clause, input1, -input2, output], [-clause, -input1, input2, output]]

def generate_buffer_clause(clause, output, input1):
    return [[-clause, -input1, output], [-clause, input1, -output]]

def generate_xnor_clause(clause, output, input1, input2):
    global comp_counter
    temp_var = comp_counter
    comp_counter += 1
    all_clauses = []
    first_clause = generate_xor_clause(clause, temp_var, input1, input2)
    second_clause = generate_not_clause(clause, output, temp_var)
    all_clauses.extend(first_clause)
    all_clauses.extend(second_clause)
    return all_clauses

def generate_multiple_clause(clause_type, clause, output, inputs):
    global comp_counter
    if len(inputs) == 1:
        return str_to_clause[clause_type](clause, output, inputs[0])
    all_clauses = []
    prev_temp = inputs[0]
    for i in range(len(inputs)-1):
        temp_var = comp_counter
        comp_counter += 1
        if i == len(inputs)-2:
            clause_expr = str_to_clause[clause_type](clause, output, prev_temp, inputs[i + 1])
        elif i == 0:
            clause_expr = str_to_clause[clause_type](clause, temp_var, inputs[i], inputs[i+1])
        else:
            clause_expr = str_to_clause[clause_type](clause, temp_var, prev_temp, inputs[i + 1])
        prev_temp = temp_var
        all_clauses.extend(clause_expr)
    return all_clauses

def split_str_number(string):
    match = re.match(r"([a-z]+)([0-9]+)", string, re.I)
    if match:
        items = match.groups()
    else:
        return [string]
    return items

str_to_clause = {'nand': generate_nand_clause, 'and': generate_and_clause, 'xor': generate_xor_clause,
                 'or': generate_or_clause, 'nor': generate_nor_clause, 'inverter': generate_not_clause,
                 'xnor': generate_xnor_clause, 'buffer': generate_buffer_clause}

class DP:
    def __init__(self, sd_path, obs_path, obs_idx):
        self.name = None
        self.comps = None
        self.inputs = None
        self.outputs = None
        self.observations = None
        self.map_comp_to_int = {}
        self.sd_path = sd_path
        self.obs_path = obs_path
        self.parse_sd()
        self.parse_obs(obs_idx)

    def parse_sd(self):
        print('parsing system...')
        global comp_counter
        with open(self.sd_path) as fp:
            sd_items = fp.read().replace('\n', '').split('.')
            self.name = sd_items[0]
            self.healthy_clauses = {}
            input_strs = sd_items[1][1:-1].split(',')
            output_strs = sd_items[2][1:-1].split(',')
            for str in input_strs:
                self.map_comp_to_int[str] = comp_counter
                comp_counter += 1
            for str in output_strs:
                self.map_comp_to_int[str] = comp_counter
                comp_counter += 1
            clause_strs = sd_items[3][2:-2].split('],[')
            for clause_str in clause_strs:
                clause_elems = clause_str.split(',')
                clause_type = split_str_number(clause_elems[0])[0]
                clause_name = clause_elems[1]
                clause_output = clause_elems[2]
                clause_inputs = clause_elems[3:]
                if clause_output not in self.map_comp_to_int.keys():
                    self.map_comp_to_int[clause_output] = comp_counter
                    comp_counter += 1
                if clause_name not in self.map_comp_to_int.keys():
                    self.map_comp_to_int[clause_name] = comp_counter
                    comp_counter += 1
                self.healthy_clauses[clause_name] = generate_multiple_clause(clause_type, self.map_comp_to_int[clause_name],
                    self.map_comp_to_int[clause_output], [self.map_comp_to_int[ci] for ci in clause_inputs])


    def parse_obs(self, obs_idx):
        print(f'parsing observation {obs_idx}')
        self.obs_clauses = []
        with open(self.obs_path) as fp:
            obs_items = list(filter(lambda l: len(l)>0, fp.read().replace('\n', '').split('.')))
            obs_item = obs_items[obs_idx]
            obs_item_split = obs_item.split(',[')
            self.parse_single_obs(obs_item_split[1][:-2].split(','))

    def parse_single_obs(self, obs):
        for elem in obs:
            if '-' in elem:
                self.obs_clauses.append([-self.map_comp_to_int[elem[1:]]])
            else:
                self.obs_clauses.append([self.map_comp_to_int[elem]])

def get_hitting_set(Lambda):
    all_sets =  list(product(*Lambda))
    for i in range(len(all_sets)):
        all_sets[i] = list(set(all_sets[i]))
    return all_sets

def find_new_hitting_set(hitting_sets, Omega, checked_hs):
    if len(hitting_sets) == 0:
        return None
    while(len(hitting_sets) > 0 and (hitting_sets[-1] in Omega or hitting_sets[-1] in checked_hs)):
        hitting_sets.pop()
    if len(hitting_sets) > 0:
        return hitting_sets.pop()
    else:
        return None

def powerset(seq):
    """
    Returns all the subsets of this set. This is a generator.
    """
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item


def minimize(dp, complement, mode):
    while True:
        not_changed = True
        for elem in complement:
            comp_copy = deepcopy(complement)
            comp_copy.remove(elem)
            if cons_check(dp, comp_copy, mode):
                complement.remove(elem)
                not_changed = False
        if not_changed:
            break
    return complement


def PDDS(DP, Omega, Lambda, mode, K=1):
    print(f'Running PDDS in {mode} mode')
    checked_hs = []
    for i in range(K):
        print(f'PDDS: Omega for stage {i}: {Omega}')
        hitting_sets = get_hitting_set(Lambda)
        hitting_set = find_new_hitting_set(hitting_sets, Omega, checked_hs)
        checked_hs.append(hitting_set)
        if hitting_set is not None:
            if cons_check(DP, hitting_set, mode):
                toadd = True
                toremove = []
                for subset in Omega:
                    if subset.issuperset(hitting_set):
                        toremove.append(subset)
                    elif subset.issubset(hitting_set):
                        toadd = False
                        break
                if toadd:
                    Omega.append(set(hitting_set))
                for item_to_remove in toremove:
                    Omega.remove(item_to_remove)
            else:
                complement = [item for item in list(DP.healthy_clauses.keys()) if item not in hitting_set]
                if mode == 'conflict':
                    temp_mode = 'diagnosis'
                else:
                    temp_mode = 'conflict'
                min_comp = minimize(DP, complement, temp_mode)
                toadd = True
                to_remove = []
                for subset in Lambda:
                    if subset.issuperset(min_comp):
                        to_remove.append(subset)
                    elif subset.issubset(min_comp):
                        toadd=False
                        break
                if toadd:
                    Lambda.append(set(min_comp))
                for remove_item in to_remove:
                    Lambda.remove(remove_item)
        else:
            return Omega

def cons_check(dp, hitting_set, mode):
    final_clauses = deepcopy(dp.obs_clauses)
    for desc in dp.healthy_clauses.values():
        final_clauses.extend(desc)
    mode_multiplier = {'diagnosis': [-1, 1], 'conflict': [1, -1]}
    for healthy_clause in dp.healthy_clauses.keys():
        if healthy_clause in hitting_set:
            final_clauses.append([mode_multiplier[mode][0]*dp.map_comp_to_int[healthy_clause]])
        else:
            final_clauses.append([mode_multiplier[mode][1]*dp.map_comp_to_int[healthy_clause]])
    # sol = pycosat.solve(final_clauses)
    g = Minisat22()
    for clause in final_clauses:
        g.add_clause(clause)
    return g.solve()

    if mode == 'diagnosis':
        return sol != 'UNSAT'
    else:
        return sol == 'UNSAT'


def stop_timer():
    global timer
    timer = False


def SDE(dp):
    global Lambda, Omega
    """
    Run Switching Diagnostic Engine
    :param dp: diagnostic problem
    :return: The set of conflicts
    """
    Lambda = [set(dp.healthy_clauses.keys())]
    Omega = [set(dp.healthy_clauses.keys())]
    step = 1
    while True :
        print('Running Switching Diagnostic Engine')
        prev_lambda = deepcopy(Lambda)
        prev_omega = deepcopy(Omega)
        PDDS(dp, Lambda, Omega, 'conflict')
        if PDDS(dp, Omega,Lambda, 'diagnosis') == True:
            break
        print(f'step {step}:\nLambda: {Lambda}\nOmega: {Omega}')
        step += 1
        if prev_lambda == Lambda and prev_omega == Omega:
            break
    return Omega


def run_PDDS(dp):
    global Lambda, Omega
    Lambda = [set(dp.healthy_clauses.keys())]
    Omega = [set(dp.healthy_clauses.keys())]
    return PDDS(dp, Omega, Lambda, 'diagnosis', K=1000)


def get_num_obs(path):
    with open(path) as fp:
        obs_items = list(filter(lambda l: len(l) > 0, fp.read().replace('\n', '').split('.')))
        return len(obs_items)


def run_exp(diag_func, seconds):
    sys_files = [f for f in os.listdir(sys_folder) if os.path.isfile(f'{sys_folder}/{f}')]
    obs_files = [f for f in os.listdir(obs_folder) if os.path.isfile(f'{obs_folder}/{f}')]
    for sys_file in sys_files:
        if 'c17' not in sys_file:
            continue
        results = []
        for obs_file in obs_files:
            if 'c17' not in sys_file:
                continue
            if sys_file[:-4] in obs_file:
                num_obs = get_num_obs(f'{obs_folder}/{obs_file}')
                for obs_idx in range(num_obs):
                    res_dict = {}
                    dp = DP(f'{sys_folder}/{sys_file}', f'{obs_folder}/{obs_file}', obs_idx)
                    SDE_run = Diagnosis_on_timer(diag_func, dp)
                    SDE_run.start()
                    time.sleep(seconds)
                    SDE_run.raise_exception()
                    SDE_run.join()
                    res_dict['system'] = sys_file[:-4]
                    res_dict['observation'] = str(obs_idx)
                    res_dict['#diagnoses'] = len(Omega)
                    res_dict['sat_solver'] = sat_solver.__name__
                    res_dict['min_cardinality'] = min([len(x) for x in Omega])
                    res_dict['time_limit'] = seconds
                    res_dict['algorithm'] = diag_func.__name__
                    results.append(res_dict)
                pd.DataFrame(results).to_csv(f'SDE_results_{sys_file[:-4]}_{seconds}_sec_{diag_func.__name__}_{sat_solver.__name__}.csv')


class Diagnosis_on_timer(threading.Thread):
    def __init__(self, diag_func, dp):
        threading.Thread.__init__(self)
        self.diag_func = diag_func
        self.dp = dp

    def run(self):
        # target function of the thread class
        try:
            self.diag_func(self.dp)
        finally:
            print('observation ended')

    def get_id(self):

        # returns id of the respective thread
        if hasattr(self, '_thread_id'):
            return self._thread_id
        for id, thread in threading._active.items():
            if thread is self:
                return id

    def raise_exception(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id,
                                                         ctypes.py_object(SystemExit))
        if res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
            print('Exception raise failure')


if __name__ == "__main__":
    time_limits = [30]
    algorithms = [run_PDDS]
    sat_solvers = [Glucose3, Minisat22, Lingeling, MapleCM, Minicard, Cadical]
    for algorithm in algorithms:
        for time_limit in time_limits:
            for sat_solver in sat_solvers:
                try:
                    SDE_results = run_exp(algorithm, time_limit)
                except Exception as e:
                    pass




