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
			
			#needleman 2 arrays based
			print("need2a"+coupleDescript+analyzeFct(need2arrays,seqs[i],seqs[j]))
			
                        #needleman 2 arrays based
			print("need1a"+coupleDescript+analyzeFct(need1array,seqs[i],seqs[j]))

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
	mem = memory_usage(proc=(fct,(seq1,seq2)),interval=0.1)

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
	for i in range(1,len1):
		for j in range(1,len2):
			# Match / mismatch ?
			match = matrix[i-1][j-1] + (1 if seq1[i]==seq2[j] else -1)

			# indel ?
			i1 = matrix[i-1][j]-1
			i2 = matrix[i][j-1]-1

			# keep best
			matrix[i][j] = max(match, i1, i2)

	# Get the list of best concurrents (end of seq1 and seq2 continue after, or end of seq2 (seq2 is inside seq1))
            #end of seq1
	best_co = matrix[-1][:]
            #end of seq2
	for i in range(len1):
	    best_co.append(matrix[i][-1])
	#Get scores
	best = max(best_co)
	if(best<0):best=0
	return 100*best/min(len1,len2)

'''
Algorithm 'Need2Arrays
Same algorithm (needleman&wunsch) but using a score matrix of size (seq1.size(),2)
as we only need the previous line to calculate the current one.
'''
def need2arrays(seq1, seq2):
	#Get sequences length
	len1 = len(seq1)
	len2 = len(seq2)

	# Prepare scoring arrays
	previous = [0]*len2
	current = [0]*len2

	# bestIN is the best case of seq2 being INSIDE seq1

	# Calculate the first line (seq2 against seq1[0])
	# Note: we can't 'allow' a gap in seq2 'first nuc'
	previous[0] = (1 if seq1[0]==seq2[0] else -1)
	for j in range(1,len2): previous[j] = -j
	bestIN = previous[len2-1]

	# Complete the needle
	for i in range(1, len1):
		# Test first nuc, a gap is allowed So it's just a (mis)match test
		current[0] = (1 if seq1[i]==seq2[0] else -1)

		# Finish the line
		for j in range(1,len2):
			scores=[0,0,0]
			# Match/mismatch ?
			scores[0] = previous[j-1] + (1 if seq1[i]==seq2[j] else -1)
			#Indels ?
			scores[1] = previous[j]-1
			scores[2] = current[j-1]-1
			#keep best
			current[j]=max(scores)
		
		# Best seq2 inside possibility ?
		if(current[j] > bestIN): bestIN=current[j]

		# move buffers
		previous = current
		current = [0]*len2
	
	# Get the best possibilities score which is inside the last line + the bestIN (pseudo last column)
	previous.append(bestIN) #Previous is last because we moved buffers. We add bestIN to it
	best=max(previous)
	if(best<0):best=0
	return 100*max(previous)/min(len1,len2) #Percent and normalize by theorical max best score (full match)

'''
Algorithm 'Need1Array'
Same algorithm (needleman&wunsch) but using a uniq array of size seq1.size()
as we only need the previous/actual line to calculate the current one.
'''
def need1array(seq1, seq2):
    #Get sequences length
    len1 = len(seq1)
    len2 = len(seq2)

    # Prepare scoring array and previous (i-1, j-1) value buffer
    scores = [0]*len1
    previous_align_score = 0
    best = 0

    # Calculate the first line (seq2[0] against seq1)
    for i in range(len1): scores[i] = (1 if seq1[i]==seq2[0] else -1)
    best = scores[-1]

    #Complete the needle
    for j in range(1, len2):
        previous_align_score = scores[0]
        # first nuc of seq1 : Is an indel (as the first is forced to be seq2[0])
        scores[0]=-j
        
        #Do other nucs
        for i in range (1, len1):
            # It can be a match, a mismatch, an indel (of seq1), an indel (of seq2)
            s = max(
                previous_align_score + (1 if seq1[i]==seq2[j] else -1), #match/mismatch
                scores[i-1]-1, #indel from seq1
                scores[i]-1, #indel from seq2
            )
            #get 'old' align
            previous_align_score = scores[i]
            #input new score
            scores[i]=s

        #best ? (end of seq1)
        if scores[-1] > best: best = scores[-1]

    #Check the best is not when seq2 ends in seq1
    bs = max(scores)
    if(bs > best): best = bs
    if(best<0):best=0

    #return best
    return 100*best/min(len1,len2)

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

