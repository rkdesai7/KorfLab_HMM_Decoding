import math
import json
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
import gzip
   
parser = argparse.ArgumentParser(description = "Return states of a genetic sequence and it's probability using forward backward decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
parser.add_argument('--state', type=str, default = "exon1", help="Name of state you want to see in the output")
parser.add_argument('--output', type=str, default = "GFF", help="Output format. Options include 'Wiggle', 'Bed', and 'GFF'")
arg = parser.parse_args()


def log_dict(probs_dict):
    log_probs = {}
    for from_state, to_probs in probs_dict.items():
        log_to_probs = {}
        for to_state, prob in to_probs.items():
            if prob > 0:
                log_prob = math.log(prob)
            else:
                log_prob = -99
            log_to_probs[to_state] = log_prob
        log_probs[from_state] = log_to_probs
    return log_probs   
def add_logspace(a, b, thresh = 40):
    a = math.exp(a)
    b = math.exp(b)
    if abs(a - b) > thresh: return max(a, b)
    if a < b: return math.log(1 + math.exp(a - b)) + b
    return math.log(1 + math.exp(b - a)) + a
def readfasta(filename):

	name = None
	seqs = []

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()
def run(states, sequence, orders, emits):
    #Forward Fill
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
                    prev_prob = math.log(forward_matrix[k][i-1])
                    prob = prev_prob+transitions[states[k]][states[j]]+emits[states[j]][seq]
                    sum += math.exp(prob)
            forward_matrix[j].append(sum)
    for i in forward_matrix:
        i = i.pop(0)
    #Backward Fill
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
    return true_probs
def find_probable_state(true_probs, states, sequence):
    decoded = []
    for i in range(len(true_probs[0])):
        state = ""
        prob = 0
        for j in range(len(states)):
            if true_probs[j][i] > prob:
                prob = true_probs[j][i]
                state = states[j]
        decoded.append((sequence[i], state))
    return decoded
def gff(dataset):
    for name, probs in dataset.items():
        ranges = []
        curr_state = None
        start_index = 0
        for i, (pos, state) in enumerate(probs):
            if i == 0:
                curr_state = state
                start_index = i + 1
                continue
            if state != curr_state:
                ranges.append((name, curr_state, start_index, i))
                start_index = i+1
            #Start new range
            curr_state = state
        if curr_state is not None:
            ranges.append((name, curr_state, start_index+1, i+1))
    return pd.DataFrame(ranges, columns = {"Name":[], "State":[], "Start":[], "End":[]})
#def wiggle():
    #return output
#def bed():
    #return output


#read in json
with open(arg.HMM, 'r') as f:
    data = json.load(f)   
#get values from json in log-scale
states = data["states"]
transitions = log_dict(data["transitions"])
emits = log_dict(data["emissions"])
orders = []
for i in states:
    dict = emits[i]
    order = len(list(dict.keys())[0]) - 1
    orders.append(order)
list_of_seq = readfasta(arg.sequence)
dataset = {}
for i in list_of_seq:
    sequence = i[1]
    name = i[0]
    true_probs = run(states, sequence, orders, emits)
    dataset[name] =  find_probable_state(true_probs, states, sequence)
    print(dataset)
#Output
if arg.output == "GFF":
    print(gff(dataset))
