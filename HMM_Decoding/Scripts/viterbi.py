import math
import json

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

#read in fasta file
sequence = 'GTATGATTTTTTAAAATTATAATTGTTTCTTTCGAAAAAAAAAATTTCATTTACAG'

#load json
with open("C:/Users/riyak/Documents/Honors_Thesis/HMMs/HMM_AlternativeSplicing.json", 'r') as f:
    data = json.load(f)

#get values from json in log-scale
o = data["order"][0] #order
states = data["states"]
transitions = log(data["transitions"])
emits = log(data["emissions"])

#initialize matrix and emmissions
matrix = []
p = math.log(1/len(states)) #initial state probability for index 0
for i in range(len(states)):
    matrix.append([(p, -1)])
temp_emit = .25

#fill matrix
n = 1
calc = [0]*len(states) #to find the previous state
for i in range(0, len(sequence)):
    
    if i < o:
        for j in range(len(states)):
            for k in range(len(states)):
                calc[k] = matrix[k][n-1][0]+transitions[states[k]][states[j]]+temp_emit
            matrix[j].append((max(calc), calc.index(max(calc))))
    else:
        seq = sequence[i-o: i+1]
        for j in range(len(states)):
            for k in range(len(states)):
                calc[k] = matrix[k][n-1][0]+transitions[states[k]][states[j]]+emits[states[j]][seq]
            matrix[j].append((max(calc), calc.index(max(calc))))      
    n+=1

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
    
#reverse output so it displays in order
output = output[::-1]
print(output)
