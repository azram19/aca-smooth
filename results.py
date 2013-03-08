#!/usr/bin/env python
import os
import re
import subprocess
import logging
import optparse
import json

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

iterations = 20

def call_command(command):
    process = subprocess.Popen(command.split(' '),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.communicate()



reBenchmark = re.compile( 'BENCHMARK: (\d+\.\d+)s' )
rePassed = re.compile( 'Test passed' )

def extract_execution_data(output):
	passed =  True if rePassed.search(output) else False
	time = reBenchmark.findall(output)
	time = time[0]

	data = {
		"time": float(time),
		"passed": passed 
	}

	return data

def extract_statistics(statistics):
	stats = {}

	for k in statistics:
		stats_for_k = statistics[k]
		
		statMax = max((s['time'] for s in stats_for_k)) 
		statMin = min((s['time'] for s in stats_for_k)) 
		statSum = sum((s['time'] for s in stats_for_k))
		statAvg = (statSum - statMin - statMax) / (iterations - 2)

		stats[k] = statAvg

	return stats

def main():
	versions = [
		{
			"comment": "Initial",
			"hash" : "abd8944109bd5d30fb3ae5ad0e8211144eda2b53"
		},
		{	"comment": "GCC compiler flags for the architecture",
			"hash" :  "f84e4c15f48379e29980b8922c7b8d36fb39c5ef"
		},
		{	"comment": "Cache calculations",
			"hash" :  "78dc56b780322cb9170ec3cd036fa6d2b1abe0c2"
		},
		{
			"comment" : "Multithreading",
			"hash" : "4832d252005ea59a25fd63522e5681fda059d161"
		}, 
		{
			"comment" : "Better compiler, more optimazations".
			"hash" : "298074ddda46988df2dd763e8ca580810616bf21"

		}
	]

	statistics = {}

	for k in versions:
		v = k['hash']

		#checkout selected revision
		output, err = call_command( "git checkout %s" % (v,) )
		print output, err

		#build the program
		call_command( "make clean" )
		call_command( "make" )

		#create a revision object
		revision = []

		#execute and take measurements ten times
		for i in xrange(iterations):
			output, _ = call_command( r"./ACA2-2013 small.vtu" ) 
			run_data = extract_execution_data( output )
			revision.append( run_data )	
 		
		statistics[k['comment']] = revision

	#clean up
	call_command( "git checkout master" )
	call_command( "make clean" )
	call_command( "make" )

	print statistics

	print json.dumps( extract_statistics( statistics ) )

if __name__ == "__main__":
    main()
