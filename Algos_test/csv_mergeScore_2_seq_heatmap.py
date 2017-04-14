#! /usr/bin/python3

import sys
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
	
	# Prepare data
	seqslabels = list(range(len(values[algo])))
	nbSeqs = len(seqslabels)
	data = np.empty([nbSeqs, nbSeqs])  #Init the data array
	fig = plt.figure() # prepare figure
	subplots = ["scores (%)","mem (MB)","time (s)"]
		#Use same limits to normalize heatmaps
	limits=[
		score_limits,
		mem_limits,
		time_limits
	]

	# For each algorithm plot the score, mem and timeA
	algo_list = list(values.keys())
	nbAlgos = len(algo_list)
	s = 0 #subplot id, must be 1 but it will directly be incremented
	for algo in algo_list:
		# For the score=0, mem=1 and time=2
		for e in range(3):	
			s += 1 #Increment subplot id
			# Prepare subplot
			ax = fig.add_subplot(nbAlgos,3,s) # plot on same line
			ax.set_title("Heatmap of "+subplots[e]+" using "+algo)

			# Fill the data array
			for i in range(nbSeqs):
				for j in range(nbSeqs):
					data[i][j] = values[algo][i][j][e] # 0 is the id of score
			
			# Create the heatmap
			heatmap = ax.pcolor(data, cmap=plt.cm.Blues, vmin=limits[e][0], vmax=limits[e][1])

			# put the major ticks at the middle of each cell
			ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
			ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

			# want a more natural, table-like display
			#ax.invert_yaxis()
			#ax.xaxis.tick_top()

			ax.set_xticklabels(seqslabels, minor=False)
			ax.set_yticklabels(seqslabels, minor=False)
			fig.colorbar(heatmap)

	# show plot
	plt.show()


	
###################################

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Expect stricly 1 argument : the csv file")
		sys.exit()

	#run main fct with csv file
	main(sys.argv[1])


