import math
import json
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
import gzip
   
parser = argparse.ArgumentParser(description = "Return states of a genetic sequence and\
                                                it's probability using viterbi decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
parser.add_argument('--state', type=str, default = "exon1", help="Name of state you want a probability graph for")
parser.add_argument('--output', type=str, default = "GFF", help="Output format. Option include 'Wiggle', 'Bed', and 'GFF'")
arg = parser.parse_args()


def log(transition_probs):
    log_transition_probs = {}
    for from_state, to_probs in transition_probs.items():
        log_to_probs = {}
        for to_state, prob in to_probs.items():
            if prob > 0:
                log_prob = math.log(prob)
            else:
                log_prob = -99
            log_to_probs[to_state] = log_prob
        log_transition_probs[from_state] = log_to_probs
    return log_transition_probs
def filter_sequence(seq):
    valid = ['A', 'G', 'C', 'T']
    for i in seq:
        if i not in valid:
            print("Sequence has an incorrect format and must only contain nucleotides (A, C, G, or T).")
            sys.exit()
    return True
#def self():
    #read in data
#def main():
    #run and return desired output
#def forward():
    #run forward fill
#def backward():
    #run backward fill
#def true():
    #calc true probs
#def graph():
    #return graph
#def wiggle():
    #return output
#def bed():
    #return output
#def gff():
    #return output


#Initialization
#read in fasta file and json
if filter_sequence(arg.sequence.upper()):
    sequence = arg.sequence.upper()
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
p = 1/len(states) #initial state probability for index 0
for i in range(len(states)):
    forward_matrix.append([p])
temp_emit = math.log(.25)
#fill matrix
for i in range(1, len(sequence)+1):
    for j in range(len(states)):
        sum = 0
        for k in range(len(states)):
            if i < orders[j]+1:
                prev_prob = math.log(forward_matrix[k][i-1])
                prob = prev_prob+transitions[states[k]][states[j]]+temp_emit
                sum += math.exp(prob)
            else:
                seq = sequence[i-orders[j]-1: i]
                try:
                    prev_prob = math.log(forward_matrix[k][i-1])
                except:
                    print(forward_matrix[k][i-1])
                prob = prev_prob+transitions[states[k]][states[j]]+emits[states[j]][seq]
                sum += math.exp(prob)
        forward_matrix[j].append(sum)
for i in forward_matrix:
    i = i.pop(0)

#Backward Fill
#Preallocate ahead of time
backward_matrix = []
for i in range(len(states)):
    temp = [0]*(len(sequence))
    backward_matrix.append(temp)
#Fill matrix
for i in range(len(sequence), -1, -1):
    if i == len(sequence):
        continue
    elif (i+1)==len(sequence):
        for j in backward_matrix:
            j[i] = 1
    else:
        for j in range(len(states)):
            sum = 0
            for k in range(len(states)):
                if (i+1) < (orders[k]+1):
                    prev = math.log(backward_matrix[k][i+1])
                    prob = transitions[states[j]][states[k]]+prev+temp_emit
                    sum += math.exp(prob)
                else:
                    seq = sequence[i-(orders[k]-1):i+2]
                    prev = math.log(backward_matrix[k][i+1])
                    prob = transitions[states[j]][states[k]]+prev+emits[states[k]][seq]
                    sum += math.exp(prob)
            backward_matrix[j][i] = sum

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

#Graph?
s = states.index(arg.state)
x = np.array(list(range(1, len(sequence)+1)))
y = np.array(true_probs[s])
plt.plot(x, y)
state_name = "Probability that each position is in state " + arg.state
plt.title(state_name)
plt.xlabel("Position")
plt.ylabel("Probability")
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
