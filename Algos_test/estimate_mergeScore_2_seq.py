#! /usr/bin/python3

import sys #Arguements gathering


'''
This script is used to test algorithms to give an estimate percent score
between 2 sequences. The main goal is to find a suitable algorithm for a GPU
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

def main(seq1, seq2):
	#output sequences
	print("Seq1: "+seq1)
	print("Seq2: "+seq2)

	#Run first algorithm (dot_cut)
	print("\ndot_cut: "+str(dot_cut(seq1,seq2))+" %")

'''
Algorithm 'Smith-watermann'
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
	matrix = [[0]*len2]*len1


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
	# Get the sequences
	if len(sys.argv) != 3:
		print("Expect strictly 2 arguments : the 2 sequences")
		sys.exit()

	# Run the main function
	main(str(sys.argv[1]), str(sys.argv[2]))

