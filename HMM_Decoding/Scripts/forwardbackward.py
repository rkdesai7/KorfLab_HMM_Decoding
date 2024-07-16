import math
import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = "Return states of a genetic sequence and\
                                                it's probability using viterbi decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
arg = parser.parse_args()

# Functions
def log(transition_probs):
    log_transition_probs = {}
    for from_state, to_probs in transition_probs.items():
        log_to_probs = {}
        for to_state, prob in to_probs.items():
            log_prob = math.log(prob)
            log_to_probs[to_state] = log_prob
        log_transition_probs[from_state] = log_to_probs
    return log_transition_probs

def calc_forward(j, k, n, matrix, transitions, emission):
    return matrix[k][n-1]+transitions[states[k]][states[j]]+emission

#Initialization
#read in fasta file and json
sequence = arg.sequence
with open(arg.HMM, 'r') as f:
    data = json.load(f)   
#get values from json in log-scale
states = data["states"]
transitions = log(data["transitions"])
emits = log(data["emissions"])
orders = []
for i in states:
    dict = emits[i]
    order = len(list(dict.keys())[0]) - 1
    orders.append(order)

#Forward Fill
#initialize matrix and emmissions
forward_matrix = []
p = math.log(1/len(states)) #initial state probability for index 0
for i in range(len(states)):
    forward_matrix.append([p])
temp_emit = math.log(.25)
#fill matrix
for i in range(1, len(sequence)):
    for j in range(len(states)):
        sum = 0
        for k in range(len(states)):
            if i < orders[j]+1:
                sum += calc_forward(j, k, i, forward_matrix, transitions, temp_emit)
            else:
                seq = sequence[i-orders[j]: i+1]
                sum += calc_forward(j, k, i, forward_matrix, transitions, emits[states[j]][seq])
        forward_matrix[j].append(sum)

#Backward Fill
#Initialize Empty
backward_matrix = []
for i in range(len(states)):
    temp = [0]*(len(sequence)+1)
    backward_matrix.append(temp)
#Fill matrix
for i in range(len(sequence), -1, -1):
    if i == len(sequence):
        for j in backward_matrix:
            j[i] = 1
    else:
        for j in range(len(states)):
            sum = 0
            for k in range(len(states)):
                if (i+1) < (orders[k]+1) or (i+1)==len(sequence):
                    sum += transitions[states[j]][states[k]]+temp_emit+backward_matrix[k][i+1]
                else:
                  #Turn this into a function
                    seq = sequence[i-(orders[k]-1):i+2]
                    sum += transitions[states[j]][states[k]]+backward_matrix[k][i+1]+emits[states[k]][seq]
            backward_matrix[j][i] = sum

#True probs
#Get rid of temp values from forward and backward matrices
for i in forward_matrix:
    i = i.pop(0)
for i in backward_matrix:
    i = i.pop()

