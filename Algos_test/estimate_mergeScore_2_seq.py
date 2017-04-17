#! /usr/bin/python3

import sys #Arguments gathering, exit script, ...
from Bio import SeqIO #Biopython, use to read fasta file sequences
import time
from memory_profiler import memory_usage #Use to get the memory usages. /!\ You may need to install memory_profiler.
import gc #garbage collector. Use to call it only at good moments

'''
usage: ./estimate_mergeScore_2_seq.py <file.fasta>
This script is used to test algorithms to give an estimate percent score
between 2 sequences. 
It wil compare each pair of a set of sequences (given by a fasta file)
The main goal is to find a suitable algorithm for a GPU
usage. Therefore, a fairly precise algorithm is required (the reverse
complement is not an obligation as we can re-calculate it),
using as little memory as possible (the computation time is less important).
The ideal algorithm would not care if the beginning of the first read and the
ending of the second one are not merged with the other read if a lot of the
reads are merged, example :
 -read1: CCCGATTCA
 -merge:  ||||| ||
 -read2:  CCGATCCAAT
 -score: high percent
'''

def main(seqs):
	# Correct the seqs
	nbSeqs = len(seqs)
	#print("\n==== Sequences ====")
	for s in range(nbSeqs):
		seqs[s] = seqs[s].upper() #To have only upper cases
		#print("[-- SEQ "+str(s)+" --]\n"+seqs[s])
	
	#Run algorithms on each couple
	gc.disable() # Disable automatic garbage collector
	#print("\n==== Algorithms run ====")
	print("algo\tseqID1\tseqID2\tscore (%)\tmem (MB)\ttime (s)")

	for i in range(nbSeqs):
		for j in range(nbSeqs):
			coupleDescript = "\t"+str(i)+"\t"+str(j)+"\t"
			#needleman based
			print("needle"+coupleDescript+analyzeFct(needle,seqs[i],seqs[j]))

			#dot_cut
			print("dot_cut"+coupleDescript+analyzeFct(dot_cut,seqs[i],seqs[j]))

'''
Tool function to run an cmp function and get (and output) its informations
as time, value, etc...
input: function, seq1, seq2
output: string "score(%)\tmem(?)\ttime(s)"
'''
def analyzeFct(fct, seq1, seq2):
	# Execute and get delta time
	t0 = time.time()
	score = fct(seq1, seq2)
	t1 = time.time()

	# Get memory usage during the function exection 
	gc.collect() #Run garbage collector
	mem = memory_usage((fct,(seq1,seq2)),interval=0.1)

	return str(round(score,2))+"\t"+str(round(max(mem),2))+"\t"+str(round(t1-t0,2))

'''
Algorithm 'Needleman&Wunsch'
Based on an needleman principe: each case is the best of
 - score[i-1][j-1] + (seq1[i]==seq2[j]?1:-1) #Match / mismatch
 - score[i-1][j]-1 #Indel
 - score[i][j-1] #Indel
The first nuc of seq2 is either a match of score 1 or a mismatch of score -1: the gap before this on the seq1 doesn't count.
After the end seq1, indels on seq2 doesn't count.
The best score is the best of the last score line and the last column
'''
def needle(seq1, seq2):
	# Get sequences length
	len1 = len(seq1)
	len2 = len(seq2)

	# Prepare scoring matrix[len1][len2]
	matrix = [[0 for j in range(len2)] for i in range(len1)]

	# Fill first column: first nuc of seq2 against all seq1 => gap is allowed, only match/mismatch of score 1/-1
	for i in range(len1): matrix[i][0] = 1 if seq1[i]==seq2[0] else -1

	# Fill first line: it is a gap so the -1 increase
	for j in range(1,len2): matrix[0][j] = -j

	#Fill the remaining of the matrix
	#Keep the best possibilities with the end of seq2 (seq2 is inside seq1)
	bestIN=0
	for i in range(1,len1):
		for j in range(1,len2):
			# Match / mismatch ?
			match = matrix[i-1][j-1] + (1 if seq1[i]==seq2[j] else -1)

			# indel ?
			i1 = matrix[i-1][j]-1
			i2 = matrix[i][j-1]-1

			# keep best
			matrix[i][j] = max(match, i1, i2)

		# Best seq2 inside seq 1 ?
		if(matrix[i][len2-1] > bestIN): bestIN = matrix[i][len2-1]
			
	# Get the list of best concurrents
	best_co = matrix[len1-1][:]
	for i in range(len1):
		best_co.append(matrix[i][-1])
	#Add the se2 is in seq1 possibility
	best_co.append(bestIN)

	# Get score
	return 100*max(best_co)/min(len1,len2)

''''
Algorithm 'Needle3Arrays
Same algorithm but using a score matrix of size (seq1.size(),2)
as we only need the previous line to calculate the current one.
'''
	


'''
Algorithm 'dot_cut'.
The principe is to increase the score when the merge is contiguous (diagonal
comparaison) but an difference make a score=0 :
	IF seq1[i] == seq2[j]: score_actual[j] = score_previous[j-1]+1
	ELSE: score_actual[j] = 0
	...do stuff...
	score_previous = score_actual

Then the final percent score is the sum of higher scores before a difference
compared to a perfect match.
This way you have only to memorize 2 'comparaison' column score and the best
scores.
When a part of seq2 has been 'scored', only check the remaining (on the right)
of seq2 to avoid to count a same parts multiple times and to do "come-backs".
'''

def dot_cut(seq1, seq2):
	#Get sequences length
	len1 = len(seq1)
	len2 = len(seq2)

	# Init the scores arrays
	previous = [0]*len2
	actual = [0]*len2
	previous_best = 0
	best_scores =  list()
	begin_len2 = 0 #From which nuc of len2 the cmp start

	# Do cmp
	for i in range(len1):
		#Cmp seq2 to actual nuc of seq1
		for j in range(begin_len2, len2):
			if seq1[i] == seq2[j]:
				if j>0: actual[j] = previous[j-1]+1
				else: actual[j] = 1
			else: actual[j]=0

		# Get the best. Do we have a drop ?
		best = max(actual)

		if best < previous_best:
			#Best found ! (the previous_best ;) )
			best_scores.append(previous_best)

			#Now, start on the next len2 nuc from the last nuc of the
			#previous_best for further check
			begin_len2 = previous.index(previous_best)+1

			#Clear the buffers
			best = 0
			actual = [0]*len2

		
		previous_best = best
		
		#Move the scores buffers
		previous = actual
		actual = [0]*len2

	# Add the last best score
	best_scores.append(previous_best)
		
	# Results
	return 100*sum(best_scores)/min(len1,len2)


'''
Python will run the following part at start.
It gets the sequences then run the main function.
'''

# Run the main function
if __name__ == "__main__":
	# Get the argument (the fasta file)
	if len(sys.argv) != 2:
		print("Expect only one argument : the fasta file of sequences")
		sys.exit()

	#Read fasta file
	seqs = []
	for record in  SeqIO.parse(open(sys.argv[1]), 'fasta'):
		seqs.append(record.seq)


	# Run the main function
	main(seqs)

