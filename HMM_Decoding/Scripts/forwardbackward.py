import math
import json
import pandas as pd
import argparse
import numpy as np
import sys
import gzip

parser = argparse.ArgumentParser(description="Return states of a genetic sequence and its probability using forward backward decoding")
parser.add_argument('HMM', type=str, help="path to json file that describes Hidden Markov Model")
parser.add_argument('sequence', type=str, help="Gene sequence you are trying to decode")
parser.add_argument('--state', type=str, default="exon1", help="Name of state you want to see in the output")
parser.add_argument('--output', type=str, default="GFF", help="Output format. Options include 'Wiggle', 'Bed', and 'GFF'")
arg = parser.parse_args()

NEG_INF = float('-inf')


def log_dict(probs_dict):
    """Convert a nested probability dict to log-space."""
    log_probs = {}
    for from_state, to_probs in probs_dict.items():
        log_to_probs = {}
        for to_state, prob in to_probs.items():
            log_to_probs[to_state] = math.log(prob) if prob > 0 else NEG_INF
        log_probs[from_state] = log_to_probs
    return log_probs


def logsumexp(log_values):
    """Numerically stable log-sum-exp over a list of log-space values."""
    log_values = [v for v in log_values if v != NEG_INF]
    if not log_values:
        return NEG_INF
    max_val = max(log_values)
    return max_val + math.log(sum(math.exp(v - max_val) for v in log_values))


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
                yield (name, seq)
                name = line[1:]
                seqs = []
            else:
                name = line[1:]
        else:
            seqs.append(line)
    yield (name, ''.join(seqs))
    fp.close()


def run(states, sequence, orders, emits):
    n_states = len(states)
    n_obs = len(sequence)
    LOG_UNIFORM = math.log(0.25)  # fallback emission for low-order positions


    # Initialise: uniform start probability, log(1/n_states)
    log_init = math.log(1.0 / n_states)

    # forward[j] is a list of length n_obs
    forward = [[NEG_INF] * n_obs for _ in range(n_states)]

    for j in range(n_states):
        # Emission at position 0
        if orders[j] == 0:
            seq_key = sequence[0]
            log_emit = emits[states[j]].get(seq_key, NEG_INF)
        else:
            log_emit = LOG_UNIFORM  # not enough context yet
        forward[j][0] = log_init + log_emit

    for i in range(1, n_obs):
        for j in range(n_states):
            # Emission: need orders[j] preceding characters + current
            if i < orders[j]:
                log_emit = LOG_UNIFORM
            else:
                seq_key = sequence[i - orders[j]: i + 1]
                log_emit = emits[states[j]].get(seq_key, NEG_INF)

            # Sum over all previous states k
            log_trans_terms = [
                forward[k][i - 1] + transitions[states[k]][states[j]]
                for k in range(n_states)
            ]
            forward[j][i] = logsumexp(log_trans_terms) + log_emit

    backward = [[NEG_INF] * n_obs for _ in range(n_states)]

    # Base case: log(1) = 0 at the last position
    for j in range(n_states):
        backward[j][n_obs - 1] = 0.0

    for i in range(n_obs - 2, -1, -1):
        for j in range(n_states):
            log_terms = []
            for k in range(n_states):
                # Emission of state k at position i+1
                if (i + 1) < orders[k]:
                    log_emit = LOG_UNIFORM
                else:
                    seq_key = sequence[(i + 1) - orders[k]: i + 2]
                    log_emit = emits[states[k]].get(seq_key, NEG_INF)

                term = (transitions[states[j]][states[k]]
                        + log_emit
                        + backward[k][i + 1])
                log_terms.append(term)
            backward[j][i] = logsumexp(log_terms)
    #decoding
    true_probs = [[0.0] * n_obs for _ in range(n_states)]

    for i in range(n_obs):
        # log denominator = log P(obs) via forward + backward at any position
        log_denom = logsumexp([forward[k][i] + backward[k][i] for k in range(n_states)])

        for j in range(n_states):
            log_num = forward[j][i] + backward[j][i]
            true_probs[j][i] = math.exp(log_num - log_denom) if log_denom != NEG_INF else 0.0

    return true_probs


def find_probable_state(true_probs, states, sequence):
    """At each position pick the state with the highest posterior probability."""
    decoded = []
    for i in range(len(true_probs[0])):
        best_state = max(range(len(states)), key=lambda j: true_probs[j][i])
        decoded.append((sequence[i], states[best_state]))
    return decoded


def gff(dataset):
    rows = []
    for name, probs in dataset.items():
        curr_state = None
        start_index = 0
        for i, (pos, state) in enumerate(probs):
            if i == 0:
                curr_state = state
                start_index = 1  # GFF is 1-based
                continue
            if state != curr_state:
                rows.append((name, curr_state, start_index, i))  # end is i (1-based inclusive)
                curr_state = state
                start_index = i + 1
        # Append final segment
        if curr_state is not None:
            rows.append((name, curr_state, start_index, len(probs)))
    return pd.DataFrame(rows, columns=["Name", "State", "Start", "End"])

with open(arg.HMM, 'r') as f:
    data = json.load(f)

states = data["states"]
transitions = log_dict(data["transitions"])
emits = log_dict(data["emissions"])

# Infer emission order from key length (order-n means key has n+1 chars)
orders = [len(next(iter(emits[s]))) - 1 for s in states]

dataset = {}
for name, sequence in readfasta(arg.sequence):
    true_probs = run(states, sequence, orders, emits)
    dataset[name] = find_probable_state(true_probs, states, sequence)

if arg.output == "GFF":
    print(gff(dataset).to_string(index=False))
