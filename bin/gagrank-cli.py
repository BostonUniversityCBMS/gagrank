# imports
import argparse
#from   gagrank import gagrank
import gagrank.gagrank as gr
import time
import sys
import os
import re
import numpy as np
from   gagrank import species
from   scipy.stats import rankdata # for getting actual GAG ranking

debug = False # variable for debugging

# main function
def main():
	start_time = time.time()
	
	################################
	# Step 1: check user arguments #
	################################
	
	print "Checking user arguments...",
	
	# initiate parser
	parser = argparse.ArgumentParser(description='Rank GAG structures based on found isotopic clusters.')
	
	# add arguments
	parser.add_argument('-c', required=True, help='GAG class (required)')
	parser.add_argument('-i', required=True, help='Input GAGfinder results file (required)')
	parser.add_argument('-r', required=False, help='Reducing end derivatization (optional)')
	parser.add_argument('-m', type=float, required=True, help='Precursor m/z (required)')
	parser.add_argument('-z', type=int, required=True, help='Precursor charge (required)')
	parser.add_argument('-s', type=int, required=False, help='Number of sulfate losses to consider (optional, default 0)')
	parser.add_argument('-a', required=False, help='Actual sequence, for testing purposes (optional)')
	
	# parse arguments
	args = parser.parse_args()
	
	print "Checking user arguments...",
	
	# get arguments into proper variables
	gClass = args.c
	rFile  = args.i
	fmla   = args.r
	pre_mz = args.m
	pre_z  = args.z
	s_loss = args.s
	actual = args.a
	
	# check to make sure a proper GAG class was added
	if gClass not in ['HS', 'CS', 'KS']:
		print "You must denote a GAG class, either HS, CS, or KS. Try 'python gagfinder.py --help'"
		sys.exit()
	
	# pick a proper class number
	if gClass == 'HS':
		cNum = 3
	elif gClass == 'CS':
		cNum = 1
	else:
		cNum = 4
	
	# check to see if the user wants to consider sulfate loss
	if not s_loss:
		s_loss = 0
	
	# print the system arguments back out to the user
	if debug:
		print "class: %s" % (gClass)
		print "GAGfinder results file: %s" % (rFile)
	
	print "Done!"
	
	# full_path will be the absolute path to the packaged GAGfragDB.db
	full_path = os.path.join(sys._MEIPASS, './GAGfragDB.db/GAGfragDB.db')
	
	# run the guts of GAGrank
	result = gr.rank_gags(rFile, gClass, cNum, fmla, pre_mz, pre_z, s_loss, actual, full_path)
	
	# for debugging
	if actual:
		actual = actual.split(':')
		for sqn in actual:
			print '\tBiRank ranking #' + str(int(len(result['br']) - rankdata(result['br'], 'max')[np.where(result['seq'] == sqn)[0][0]] + 1)) + ' to #' + str(int(len(result['br']) - rankdata(result['br'], 'min')[np.where(result['seq'] == sqn)[0][0]] + 1)) + ' out of ' + str(len(result['seq']))
	
	#########################
	# Step 9: write to file #
	#########################
	
	# write result to file
	gr.write_result_to_file(rFile, result)
	
	print "Finished!"
	print time.time() - start_time

main()
