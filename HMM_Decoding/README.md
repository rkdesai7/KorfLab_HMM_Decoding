# Overview

This is a project meant to Hidden Markov Models of exon and intron states of genetic sequences. 

### Data

Example Hidden Markov Models in the form of a .json can be found in the Data folder. The term "order" describes the length of sequences used in calculating the emmission probability. For example, an order of 1 would indicate emmission probabilities for sequences of length 2 (AA = .1, AG = .2, ...). 

### Scripts

The scripts folder contains python programs of decoding implementations. Currently, viterbi is implemented.

### In Progress

Currently, I am working on implementing forward backward decoding, as well as a more streamlined FASTA file reading. 
