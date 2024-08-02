import math
import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = "Return states of a genetic sequence and\
                                                it's probability using viterbi decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
arg = parser.parse_args()

# Function to convert probabilities to log scale
def log(transition_probs):
    log_transition_probs = {}
    for from_state, to_probs in transition_probs.items():
        log_to_probs = {}
        for to_state, prob in to_probs.items():
            log_prob = math.log(prob)
            log_to_probs[to_state] = log_prob
        log_transition_probs[from_state] = log_to_probs
    return log_transition_probs

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
#initialize matrix and emmissions
matrix = []
p = math.log(1/len(states)) #initial state probability for index 0
for i in range(len(states)):
    matrix.append([(p, -1)])
temp_emit = .25

#fill matrix
n = 1
calc = [0]*len(states) #to find the previous state
for i in range(0, len(sequence)+1):
    for j in range(len(states)):
        for k in range(len(states)):
            if i < (orders[j]+1):
                calc[k] = matrix[k][n-1][0]+transitions[states[k]][states[j]]+temp_emit
            else:
                seq = sequence[i-orders[j]: i+1]
                calc[k] = matrix[k][n-1][0]+transitions[states[k]][states[j]]+emits[states[j]][seq]
        matrix[j].append((max(calc), calc.index(max(calc))))
    n += 1

#Run viterbi
length = len(sequence)
output = [] #to store decoding
#determine starting point state
calc = [0]*len(states)
for i in range(len(states)):
    calc[i] = matrix[i][length][0]
num = calc.index(max(calc))
state = states[num]
prob = matrix[num][length][0]
prev_state = matrix[num][length][1]
output.append((sequence[length-1], state, prob))
#run for rest of bases
for i in range(1, length):
    #extract current state, probability, and base
    state = states[prev_state]
    prob = matrix[prev_state][(i+1)*-1][0]
    base = sequence[(i+1)*-1]
    #add to output
    output.append((base, state, prob))
    #update traceback state
    prev_state = matrix[prev_state][(i+1)*-1][1]

#Display Results
output = output[::-1]
results = pd.DataFrame(output, columns = ["Base", "State", "Log Probability"])
print(results)
