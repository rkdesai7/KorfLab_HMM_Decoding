import argparse
import gzip
import json


def read_fasta(filename):

	label = None
	seq = []

	fp = None
	if    filename == '-':         fp = sys.stdin
	elif filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seq) > 0:
				seq = ''.join(seq)
				yield(label, seq)
				label = line[1:]
				seq = []
			else:
				label = line[1:]
		else:
			seq.append(line)
	yield(label, ''.join(seq))
	fp.close()

def read_mm(file, order):
	counts = {}
	total_len = 0
	total_seq = 0
	for defline, seq in read_fasta(file):
		if 'N' in seq: continue
		total_len += len(seq)
		total_seq += 1
		for i in range(len(seq) - order):
			ctx = seq[i:i+order]
			nt = seq[i+order]
			if ctx not in counts: counts[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
			counts[ctx][nt] += 1
	
	model = {}
	for ctx in counts:
		total = sum(counts[ctx].values())
		for nt in counts[ctx]: model[f'{ctx}{nt}'] = counts[ctx][nt] / total
	
	return model, total_len/total_seq


parser = argparse.ArgumentParser(description='HMM trainer')
parser.add_argument('exon', help='fasta file')
parser.add_argument('intron', help='fasta file')
parser.add_argument('-e', '--emm', type=int, default=3,
	help='Markov model order [%(default)i])')
parser.add_argument('-i', '--imm', type=int, default=3,
	help='Markov model order [%(default)i])')
parser.add_argument('-d', '--don', type=int, default=5,
	help='donor length [%(default)i])')
parser.add_argument('-a', '--acc', type=int, default=6,
	help='acceptor length [%(default)i])')
parser.add_argument('--eet', type=float, default=0.98,
	help='exon to exon transition probability [%(default).3f])')
arg = parser.parse_args()

## Markov model for Exon ##
emm, elen = read_mm(arg.exon, arg.emm)
imm, ilen = read_mm(arg.intron, arg.imm)


# counts
don = [ {'A':0, 'C':0, 'G':0, 'T':0} for i in range(arg.don)]
acc = [ {'A':0, 'C':0, 'G':0, 'T':0} for i in range(arg.acc)]
total = 0
for defline, seq in read_fasta(arg.intron):
	if 'N' in seq: continue               # some introns have Ns
	if not seq.startswith('GT'): continue # canonical only
	if not seq.endswith('AG'): continue   # canonical only
	dseq = seq[:arg.don]
	aseq = seq[-arg.acc:]
	for i, nt in enumerate(dseq): don[i][nt] += 1
	for i, nt in enumerate(aseq): acc[i][nt] += 1
	total += 1

# probs
for i in range(arg.don):
	for nt in don[i]:
		don[i][nt] /= total
for i in range(arg.acc):
	for nt in acc[i]:
		acc[i][nt] /= total

hmm = {}

# states
hmm['states'] = {}
hmm['states']['exon1'] = arg.emm
for i in range(arg.don): hmm['states'][f'don-{i}'] = 0
hmm['states']['intron'] = arg.imm
for i in range(arg.acc): hmm['states'][f'acc-{i}'] = 0
hmm['states']['exon2'] = arg.emm

# transitions
hmm['transitions'] = {}
for s1 in hmm['states']:
	if s1 not in hmm['transitions']: hmm['transitions'][s1] = {}
	for s2 in hmm['states']:
		if s2 not in hmm['transitions'][s1]: hmm['transitions'][s1][s2] = 0
hmm['transitions']['exon1']['exon1'] = arg.eet
hmm['transitions']['exon1']['don-0'] = 1 - arg.eet
for i in range(arg.don -1):
	hmm['transitions'][f'don-{i}'][f'don-{i+1}'] = 1.0
hmm['transitions'][f'don-{arg.don -1}']['intron'] = 1.0
hmm['transitions']['intron']['intron'] = 1 - (1/ilen)
hmm['transitions']['intron']['acc-0'] = 1 / ilen
for i in range(arg.acc -1):
	hmm['transitions'][f'acc-{i}'][f'acc-{i+1}'] = 1
hmm['transitions'][f'acc-{arg.acc -1}']['exon2'] = 1
hmm['transitions']['exon2']['exon2'] = 1

# emissions
hmm['emissions'] = {}
hmm['emissions']['exon1'] = emm
for i in range(arg.don): hmm['emissions'][f'don-{i}'] = don[i]
hmm['emissions']['intron'] = imm
for i in range(arg.acc): hmm['emissions'][f'acc-{i}'] = acc[i]
hmm['emissions']['exon2'] = emm

print(json.dumps(hmm, indent=4))

