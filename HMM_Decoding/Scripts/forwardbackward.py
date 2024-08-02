import math
import json
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description = "Return states of a genetic sequence and\
                                                it's probability using viterbi decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
parser.add_argument('--state', type=str, default = "exon", help="Name of state you want a probability graph for")

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
#Preallocate ahead of time
forward_matrix = []
p = math.log(1/len(states)) #initial state probability for index 0
for i in range(len(states)):
    forward_matrix.append([p])
temp_emit = math.log(.25)
#fill matrix
for i in range(1, len(sequence)+1):
    for j in range(len(states)):
        sum = 0
        for k in range(len(states)):
            if i < orders[j]+1:
                sum += calc_forward(j, k, i, forward_matrix, transitions, temp_emit)
            else:
                seq = sequence[i-orders[j]-1: i]
                sum += calc_forward(j, k, i, forward_matrix, transitions, emits[states[j]][seq])
        forward_matrix[j].append(sum)

#Backward Fill
#Initialize Empty
#Preallocate ahead of time
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
#Get rid of beginning and ending values from forward and backward matrices
#Don't need to pop
for i in forward_matrix:
    i = i.pop(0)
for i in backward_matrix:
    i = i.pop()

#Decode
#Calculate new probabiliities
true_probs = []
for i in range(len(states)):
    true_probs.append([])
for i in range(len(forward_matrix[0])):
    denominator = 0
    for k in range(len(states)):
        denominator += forward_matrix[k][i]*backward_matrix[k][i]
    for j in range(len(states)):
        numerator = forward_matrix[j][i]*backward_matrix[j][i]
        prob = numerator/denominator
        true_probs[j].append(prob)

#Could make a graph?
s = states.index(arg.state)
x = np.array(list(range(1, len(sequence)+1)))
y = np.array(true_probs[s])
plt.plot(x, y)
state_name = "Probability that each position is in state " + arg.state
plt.title(state_name)
plt.xlabel("Position")
plt.ylabel("Probability")
print(arg.state)
plt.show()

#Find the maximums
decoded = []
for i in range(len(true_probs[0])):
    state = ""
    prob = 0
    for j in range(len(states)):
        if true_probs[j][i] > prob:
            prob = true_probs[j][i]
            state = states[j]
    decoded.append((sequence[i], state, prob))
for i in decoded:
    print(i)

