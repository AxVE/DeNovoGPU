#! /usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import csv

def main(csvpath):
	#Prepare storage 
	# map [algo_name] -> array2D [seqID1][seqID2] -> list[score,mem,time]
	values = {}
	score_limits = [0.0,100]
	mem_limits = [102400.0,0.0]
	time_limits = [86400.0,0.0]

	# read csv file
	with open(csvpath) as csvfile:
		reader = csv.DictReader(csvfile, delimiter='\t')
		#Headers are 'algo', 'seqID1', 'seqID2', 'score (%)', 'mem (MB)' and 'time (s)'
		for row in reader:
			algo = row['algo']
			id1 = int(row['seqID1'])
			id2 = int(row['seqID2'])
			# New algo ? seq id ? etc ?
			if algo not in values:
				values[algo] = [[[-1,-1,-1]]]
			while len(values[algo]) <= id1:
				values[algo].append([[-1,-1,-1]])
			while len(values[algo][id1]) <= id2:
				values[algo][id1].append([-1,-1,-1])

			# Get data
			score = float(row['score (%)'])
			mem = float(row['mem (MB)'])
			time = float(row['time (s)'])
			values[algo][id1][id2] = [score,mem,time]

			# Get min/max
				# mem
			if mem < mem_limits[0]: mem_limits[0] = mem
			if mem > mem_limits[1]: mem_limits[1] = mem
				# time
			if time < time_limits[0]: time_limits[0] = time
			if time > time_limits[1]: time_limits[1] = time
				
	# Test
	print("Limits:")
	print("\tscores: "+str(score_limits))
	print("\tmem: "+str(mem_limits))
	print("\ttimes: "+str(time_limits))

	# Prepare figure
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_title("Heatmap test")

	# Prepare data
	seqslabels = list(range(len(values[algo])))
	nbSeqs = len(seqslabels)
	data_score = np.empty([nbSeqs, nbSeqs])  #Init the data array
	
	# Example with "needle" algorithm
	algo = "needle"
		# Fill the data array
	for i in range(nbSeqs):
		for j in range(nbSeqs):
			data_score[i][j] = values[algo][i][j][0] # 0 is the id of score

	heatmap = ax.pcolor(data_score, cmap=plt.cm.Blues)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data_score.shape[0])+0.5, minor=False)
	ax.set_yticks(np.arange(data_score.shape[1])+0.5, minor=False)

	# want a more natural, table-like display
	#ax.invert_yaxis()
	#ax.xaxis.tick_top()

	ax.set_xticklabels(seqslabels, minor=False)
	ax.set_yticklabels(seqslabels, minor=False)

	# show plot
	plt.show()


	
###################################

if __name__ == "__main__":
	main("reads_test2.csv")


