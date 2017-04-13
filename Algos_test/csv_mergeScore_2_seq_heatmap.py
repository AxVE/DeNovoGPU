#! /usr/bin/python3

import matplotlib.pyplot as plt
import csv

def main(csvpath):
	#Prepare storage 
	# map [algo_name] -> array2D [seqID1][seqID2] -> map[score,mem,time]
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


if __name__ == "__main__":
	main("reads_test2.csv")


