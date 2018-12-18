#!/usr/bin/python

#################################
# Step 0: imports and functions #
#################################

print "Importing modules...",

import time # for timing functions
start_time = time.time()

import re # for converting chemical formulae to dictionaries and vice versa
import sys # for getting user inputs
import argparse # for getting user inputs
import sqlite3 as sq # for accessing the database
import networkx as nx # for network analysis
import scipy as sp # for getting actual GAG ranking
debug = False # variable for debugging

# import numpy
try:
	import numpy as np
except:
	print "You need to install the numpy module to use this script. Please install and try again."
	sys.exit()

# import scipy
try:
	#import scipy as sp
	from scipy.stats import rankdata
	from scipy.stats import fisher_exact as fe_test
except:
	print "You need to install the scipy module to use this script. Please install and try again."
	sys.exit()

# get local packages
from species import *

### CLASSES AND FUNCTIONS ###

# dict for going between monosaccharide abbreviations
ms_hash = {'D': 'dHexA', 'U': 'HexA', 'X': 'Hex', 'N': 'HexN', 'dHexA': 'D', 'HexA': 'U', 'Hex': 'X', 'HexN': 'N'}

# function for converting dictionary to chemical formula or composition
def dict2fmla (dt, type):
	fs = '' # formula string
	
	# check whether chemical formula or composition
	if type == 'formula':
		symbols = elems
	elif type == 'composition':
		symbols = ['D', 'U', 'X', 'N', 'A', 'S']
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	for sym in symbols:
		if sym in dt:
			if dt[sym] > 0:
				if dt[sym] > 1:
					fs += sym + str(dt[sym])
				else:
					fs += sym
	
	# return
	return fs

# function for converting chemical formula or composition to dictionary
def fmla2dict (fm, type):
	# check whether chemical formula or composition
	if type == 'formula':
		symbols = elems
		dt = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0} #, 'Na':0, 'K':0, 'Li':0, 'Mg':0, 'Ca':0}
	elif type == 'composition':
		symbols = ['D', 'U', 'X', 'N', 'A', 'S']
		dt = {'D':0, 'U':0, 'X':0, 'N':0, 'A':0, 'S':0}
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fm.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in dt: # invalid symbol entered
			if q[0] not in symbols:
				print "Invalid chemical formula entered."
				sys.exit()
			else:
				dt[q[0]] = 0
		
		if q[1] == '': # only one of this atom
			dt[q[0]] += 1
		else:
			dt[q[0]] += int(q[1])
	
	# return
	return dt

# function for shuffling
def choose_iter(elements, length):
	for i in xrange(len(elements)):
		if length == 1:
			yield (elements[i],)
		else:
			for next in choose_iter(elements[i+1:len(elements)], length-1):
				yield (elements[i],) + next

# function for getting a shuffled list
def choose(l, k):
	return list(choose_iter(l, k))

# function for getting information about the precursor
def get_precursor(charge, mz, cursor, deriv_wt, weights, gag_class, adct=None, n_adct=None):
	mass = (mz*abs(charge)) - (charge*weights['H'])
	
	# get proper precursor mass to test
	if not adct:
		n_adct    = 0
		test_mass = mass - deriv_wt
	else:
		test_mass = mass - deriv_wt - (n_adct*weights[adct]) + (n_adct*weights['H'])
	
	# get precursor info from database
	cursor.execute('''SELECT   cpm.id, f.value, p.value, f.monoMass
   	       	 		  FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
	 	       	      WHERE    p.id = cpm.pId 
    	       		  AND      f.id = p.fmId
        	     	  AND      cpm.cId = ?
   	        	 	  ORDER BY ABS(f.monoMass - ?) ASC
	       	     	  LIMIT 1;''', (gag_class, test_mass))
	row = cursor.fetchone()
	
	# return
	return row

# function for getting the reducing and non-reducing end info
def get_ends(pd, gag_class):
	n = pd['D'] + pd['U'] + pd['X'] + pd['N'] # length of GAG
	
	# check if dHexA exists (has to be NR end)
	if pd['D'] > 0:
		nonred = 'D'
		
		if pd['U'] == pd['N']: # we know that the reducing end is HexA because (HexA+dHexA > HexN)
			redend = 'U'
		else: # we know that the reducing end is HexN because (HexA+dHexA == HexN)
			redend = 'N'
	else:
		if gag_class == 4: # KS
			if pd['N'] > pd['X']: # we know that both ends are HexN
				nonred = 'N'
				redend = 'N'
			elif pd['N'] < pd['X']: # we know that both ends are Hex
				nonred = 'X'
				redend = 'X'
			else: # we cannot know which end is which just yet
				nonred = '?'
				redend = '?'
		else: # HS or CS
			if pd['N'] > pd['U']: # we know that both ends are HexN
				nonred = 'N'
				redend = 'N'
			elif pd['N'] < pd['U']: # we know that both ends are HexA
				nonred = 'U'
				redend = 'U'
			else: # we cannot know which end is which just yet
				nonred = '?'
				redend = '?'
	
	# return
	return [nonred, redend, n]

# function for getting GAG chemical name from backbone, ac positions, and so3 positions
def generateGAGstring(bbone, adict=None, sdict=None):
	monos = [] # list to store each mono
	
	# range through length
	for i in range(len(bbone)):
		j  = str(i+1) # dictionary keys are strings
		tf = False # boolean for whether to add position on HexN or not
		
		if bbone[j] == 'D': # current mono is dHexA
			cur = 'dHexA'
		elif bbone[j] == 'U': # current mono is HexA
			cur = 'HexA'
		elif bbone[j] == 'X': # current mono is Hex
			cur = 'Hex'
		else: # current mono is HexN
			cur = 'HexN'
			tf  = True
		
		# there is at least one acetyl group in this GAG and it is on this mono
		if adict and j in adict:
			cur += 'Ac'
		
		# there is at least one sulfate group in this GAG and it is on this mono
		if sdict and j in sdict:
			for p in sorted(sdict[j]):
				# 
				if p == '2' and tf:
					cur += 'S'
				else:
					cur += p + 'S'
		
		monos.append(cur)
	
	gstring = ''
	
	for m in monos:
		gstring += m + '-'
	
	return gstring[0:len(gstring)-1]

# function for making a GAG string a GAG dict
def gstring2gdict(gcl, gstring):
	# split string into monos
	monos = gstring.split('-')
	
	# GAG dictionary
	gdict = {'backbone': {}, 'Ac':{}, 'SO3':{}}
	
	# go through each mono
	for m in range(len(monos)):
		j = str(m+1)
		
		# string manipulation
		if monos[m][:5] == 'dHexA':
			curm = 'dHexA'
			mods = monos[m][5:]
		elif monos[m][:4] == 'HexA' or monos[m][:4] == 'HexN':
			curm = monos[m][:4]
			mods = monos[m][4:]
		else: # Hexose
			curm = 'Hex'
			mods = monos[m][3:]
		
		# add to backbone
		gdict['backbone'][j] = ms_hash[curm]
		
		# temporary dicts
		if curm == 'HexN':
			tadict = {'2': 0}
			tsdict = {}
			
			# add positions to temporary SO3 dict
			for p in modlocs[gcl]['N']['SO3']:
				tsdict[str(p)] = 0
			
			# actually add modifications
			if 'Ac' in mods: # acetyl present
				tadict['2'] = 1
			
			if mods != '':
				start = 0 # to determine where to start
				if mods[0] == 'S':
					tsdict['2'] = 1
					start = 1 # skip the first letter
				
				for x in [mods[i:i+2] for i in range(start, len(mods), 2)]:
					if x == 'Ac':
						tadict['2'] = 1
					else:
						tsdict[x[0]] = 1
				
			# add dicts to gdict
			gdict['Ac'][j]  = tadict
			gdict['SO3'][j] = tsdict
		else:
			tsdict = {}
			
			# add positions to temporary SO3 dict
			for p in modlocs[gcl][ms_hash[curm]]['SO3']:
				tsdict[str(p)] = 0
			
			# actually add modifications
			for x in [mods[i:i+2] for i in range(0, len(mods), 2)]:
				tsdict[x[0]] = 1
			
			# add dicts to gdict
			gdict['SO3'][j] = tsdict
		
	return gdict

# function for encoding BiRank algorithm
def BiRank(G, alpha=0.85, beta=0.85, n_iter=1000, s_0=None, f_0=None):
	# check if G is a bipartite graph
	if not nx.is_bipartite(G):
		#print 'The graph for BiRank must be bipartite'
		return 0
	
	# get the two sides of the bipartite graph
	grp1, grp2 = nx.bipartite.sets(G)
	grp1 = list(grp1)
	grp2 = list(grp2)
	
	# get list of sequences and fragments
	if 'frag' in grp1[0]:
		fragments = grp1
		sequences = grp2
	else:
		sequences = grp1
		fragments = grp2
	
	# add initial sequence info
	if s_0 is None:
		s_0 = []
		for idx in sequences:
			s_0.append(1./len(sequences))
	elif type(s_0) is dict:
		tmp = []
		for idx in sequences:
			tmp.append(s_0[idx])
		
		s_0 = np.array(tmp)
		
	seqs = [s_0]
	
	# add initial fragment info
	if f_0 is None:
		f_0 = []
		for idx in fragments:
			f_0.append(1./len(fragments))
	elif type(f_0) is dict:
		tmp = []
		for idx in fragments:
			tmp.append(f_0[idx])
		
		f_0 = np.array(tmp)
		
	frgs = [f_0]
	
	# get degree info
	a_deg = nx.degree(G, weight='weight')
	s_deg = []
	f_deg = []
	
	deg = a_deg
	
	# populate sequence degrees and fragment degrees
	for s in sequences:
		s_deg.append(a_deg[s] ** -0.5)
	
	for f in fragments:
		f_deg.append(a_deg[f] ** -0.5)
	
	# make degree matrices
	Df = sp.sparse.diags(f_deg)
	Ds = sp.sparse.diags(s_deg)
	
	# get weight matrix
	W = nx.bipartite.biadjacency_matrix(G, sequences, fragments)
	
	# get symmetrically normalized version of weight matrix
	S = Ds * W * Df
	
	# iterate
	go = True
	k  = 1 # for index control
	
	while go:
		# update ranks
		cur_seq = (alpha * S * frgs[k-1]) + ((1 - alpha) * s_0)
		cur_frg = (beta * S.transpose() * seqs[k-1]) + ((1 - beta) * f_0)
		
		if np.array_equal(cur_seq, seqs[len(seqs)-1]) or np.sum(abs(cur_seq - seqs[len(seqs)-1])) < 1.e-10 or k == n_iter:
			go = False
		else:
			# add new rankings to the list of rankings
			seqs.append(cur_seq)
			frgs.append(cur_frg)
			k += 1 # add to iterator
	
	# get back into dictionary format
	s_dict, f_dict = ({}, {})
	
	for idx in range(len(sequences)):
		s_dict[sequences[idx]] = seqs[len(seqs)-1][idx]
	
	for idx in range(len(fragments)):
		f_dict[fragments[idx]] = frgs[len(frgs)-1][idx]
	
	# return
	return s_dict, f_dict

# return a sequence's likelihood score
def scoreSeq(seq, scf):
	lik = 1 # likelihood score
	monos = seq.split('-') # get each individual monosaccharide
	for m in monos: # loop through monosaccharides
		val = 1
		if 'N3' in m or 'N6' in m or m[len(m)-1] == 'N': # free amine
			val -= 0.6
		if '3S' in m and '6S' not in m: # 3S is more rare than 6S
			val -= 0.3
		
		lik *= val
	
	# return scaled likelihood score
	return lik ** scf

# function for writing to file
def write_result_to_file(input_path, scores):
	print "Printing output to file...",
	
	# write to file
	oFile = input_path[:-4] + '_GAGrank_results.tsv'
	f     = open(oFile, 'w')
	f.write("Sequence\tGAGrank score\n")
	
	for q in scores:
		out = str(q[0]) + '\t' + str(q[1]) + '\n'
		f.write(out)
	
	f.close()
	
	print "Done!"

print "Done!"

################################
# Step 1: check user arguments #
################################

print "Checking user arguments...",

# initiate parser
parser = argparse.ArgumentParser(description='Find isotopic clusters in GAG tandem mass spectra.')

# add arguments
parser.add_argument('-c', required=True, help='GAG class (required)')
parser.add_argument('-i', required=True, help='Input GAGfinder results file (required)')
parser.add_argument('-r', required=False, help='Reducing end derivatization (optional)')
parser.add_argument('-m', type=float, required=False, help='Precursor m/z (optional, but must be in mzML file)')
parser.add_argument('-z', type=int, required=False, help='Precursor charge (optional, but must be in mzML file)')
parser.add_argument('-s', type=int, required=False, help='Number of sulfate losses to consider (optional, default 0)')
parser.add_argument('-a', required=False, help='Actual sequence, for testing purposes')

# parse arguments
args = parser.parse_args()

# get arguments into proper variables
gClass = args.c
rFile  = args.i
fmla   = args.r
pre_mz = args.m
pre_z  = args.z
s_loss = args.s
actual = args.a

# scaling factors
p1 = 5.4 # for scaling edge width
p2 = 5.1 # for scaling fragments' prior probabilities
p3 = 0.4 # for scaling sequences' prior probabilities

# get values ready for reducing end derivatization and reagent
df    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
atoms = ['C','H','O','N','S']

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

# parse the reducing end derivatization formula
if fmla:
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fmla.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in atoms: # invalid symbol entered
			print "Invalid chemical formula entered. Please enter only CHONS. Try 'python gagfinder.py --help'"
			sys.exit()
		else:
			if q[1] == '': # only one of this atom
				df[q[0]] += 1
			else:
				df[q[0]] += int(q[1])

# get derivatization weight
wt = {'C':  12.0,
      'H':  1.0078250322,
      'O':  15.994914620,
      'N':  14.003074004,
      'S':  31.972071174,
      'Na': 22.98976928,
      'K':  38.96370649,
      'Li': 7.01600344,
      'Ca': 39.9625909,
      'Mg': 23.98504170}
dw = 0
for q in df:
	dw += df[q] * wt[q]

# check to see if the user wants to consider sulfate loss
if not s_loss:
	s_loss = 0

# print the system arguments back out to the user
if debug:
	print "class: %s" % (gClass)
	print "GAGfinder results file: %s" % (rFile)
	
	formula = ''
	for key in df:
		val = df[key]
		if val > 0:
			formula += key
			if val > 1:
				formula += str(val)
	print "atoms in reducing end derivatization: %s" % (formula)

print "Done!"

###########################################################
# Step 2: Load GAGfinder results and connect to GAGfragDB #
###########################################################

print "Loading GAGfinder results file...",

# load G scores
try:
	res = open(rFile, 'r')
except:
	print "\nIncorrect or missing file. Please try again."
	sys.exit()

res.readline() # skip first line

# dictionaries
gdict    = {} # store G scores for each nonspecific fragment
fragmap  = {} # map nonspecific fragments to specific fragment/charge state pair
mapfrag  = {} # map specific fragment/charge state pair to nonspecific fragments
f_chgs   = {} # keys are unique fragments and values are each charge state
allfrags = [] # list of all fragments (similar to GAGs

# go through lines
line_ct = 0
for line in res.readlines():
	if line_ct >= 68:
		break
	
	cells = line.strip('\n').split('\t') # split string to get values out
	
	if 'M' in cells[3]: # full molecule fragments are not helpful
		continue
	
	# get G score, fragments, and charge
	G   = float(cells[1])/(float(cells[4]) ** p2)
	fgs = cells[3].split('; ')
	chg = int(cells[2])
	
	# add G score to gdict
	gdict['frag'+str(line_ct)] = G
	
	# add fragment to allfrags list
	allfrags.append('frag'+str(line_ct))
	
	# go through each fragment
	for frg in fgs:
		# add to fragment map
		fragmap[(frg, chg)] = 'frag'+str(line_ct)
		
		# add to map fragment
		if 'frag'+str(line_ct) not in mapfrag:
			mapfrag['frag'+str(line_ct)] = [(frg, chg)]
		else:
			mapfrag['frag'+str(line_ct)].append((frg, chg))
		
		if frg in f_chgs: # this fragment hasn't been added yet
			f_chgs[frg].append(chg)
		else: # this fragment has been added
			f_chgs[frg] = [chg]
	
	line_ct += 1

# convert to a numpy array
gsc = []
for key in gdict:
	gsc.append((key, gdict[key]))

dty = np.dtype([('frags', 'object'), ('G', 'float')])
gsc = np.array(gsc, dtype=dty)
gsc = np.sort(gsc, order='G')[::-1]

# connect to GAGfragDB
conn = sq.connect('../lib/GAGfragDB.db')
c    = conn.cursor()

print "Done!"

######################################
# Step 3: Find precursor composition #
######################################

print "Determining precursor composition...",

# calculate precursor mass
pre_mass = (pre_mz*abs(pre_z)) - (pre_z*wt['H'])

test_mass = pre_mass - dw
	
# get precursor info
c.execute('''SELECT   p.value
   	         FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
       	     WHERE    p.id = cpm.pId 
           	 AND      f.id = p.fmId
             AND      cpm.cId = ?
   	         ORDER BY ABS(f.monoMass - ?) ASC
       	     LIMIT 1;''', (cNum, test_mass))
row = c.fetchone()

# place precursor info into variable
pComp = row[0]

print "Done!"

#########################################################
# Step 4: Get reducing end/non-reducing end information #
#########################################################

print "Determining reducing end/non-reducing end information...",

# convert composition into a dictionary
pDict = fmla2dict(pComp, 'composition')
n_pre = pDict['D'] + pDict['U'] + pDict['X'] + pDict['N']

# check if dHexA exists (has to be NR end)
if pDict['D'] > 0:
	NR = 'D'
	
	if pDict['U'] == pDict['N']: # we know that the reducing end is HexA because (HexA+dHexA > HexN)
		RE = 'U'
	else: # we know that the reducing end is HexN because (HexA+dHexA == HexN)
		RE = 'N'
else:
	if cNum == 4: # KS
		if pDict['N'] > pDict['X']:
			NR = 'N'
			RE = 'N'
		elif pDict['N'] < pDict['X']:
			NR = 'X'
			RE = 'X'
		else:
			NR = '?'
			RE = '?'
	else: # HS or CS
		if pDict['N'] > pDict['U']: # we know that both ends are HexN
			NR = 'N'
			RE = 'N'
		elif pDict['N'] < pDict['U']: # we know that both ends are HexA
			NR = 'U'
			RE = 'U'
		else: # we cannot know which end is which just yet
			NR = '?'
			RE = '?'

print "Done!"

#########################################
# Step 5: Set up possible GAG backbones #
#########################################

print "Generating possible backbones and modification positions...",

# variables to keep possibilities
backbones = []
all_Ac    = []
all_SO3   = []

# generate backbones and modification possibile locations
if gClass == 'KS': # working with KS
	poss = ['X','N'] # possible monosaccharides
	
	if NR == '?': # unsure about end monosaccharides
		for i in poss: # we need two backbones
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			odd  = i # odd-numbered monosaccharide
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			# go through each position in backbone
			for j in range(n_pre):
				if (j+1) % 2 == 0: # even
					cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
					
					if even == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-6')
						
				else: # odd
					cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
					
					if odd == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-6')
			
			# add current backbone and modifications
			backbones.append(cur_back)
			all_Ac.append(cur_Ac)
			all_SO3.append(cur_SO3)
	else: # we know end monosaccharides
		cur_back = {} # current backbone
		cur_Ac   = [] # current acetyl modification positions
		cur_SO3  = [] # current sulfate modification positions
		
		odd  = NR
		even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
		
		# go through each position in backbone
		for j in range(n_pre):
			if (j+1) % 2 == 0: # even
				cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
				
				if even == 'N':
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-6')
			else: # odd
				cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
				
				if odd == 'N':
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-6')
		
		# add current backbone and modifications
		backbones.append(cur_back)
		all_Ac.append(cur_Ac)
		all_SO3.append(cur_SO3)
elif gClass == 'CS': # working with CS
	poss = ['U','N'] # possible monosaccharides
	
	if NR == '?': # unsure about end monosaccharides
		for i in poss: # we need two backbones
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			odd  = i # odd-numbered monosaccharide
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			# go through each position in backbone
			for j in range(n_pre):
				if (j+1) % 2 == 0: # even
					cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
					
					if even == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-4')
						cur_SO3.append(str(j+1)+'-6')
					else: # HexA
						cur_SO3.append(str(j+1)+'-2')
				else: # odd
					cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
					
					if odd == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-4')
						cur_SO3.append(str(j+1)+'-6')
					else: # HexA
						cur_SO3.append(str(j+1)+'-2')
			
			# add current backbone and modifications
			backbones.append(cur_back)
			all_Ac.append(cur_Ac)
			all_SO3.append(cur_SO3)
	else: # we know end monosaccharides
		cur_back = {} # current backbone
		cur_Ac   = [] # current acetyl modification positions
		cur_SO3  = [] # current sulfate modification positions
		
		# do position 1 first
		cur_back['1'] = NR
		
		# change odd to U if it's a dHexA
		if NR == 'D':
			odd = 'U'
		else:
			odd = NR
		
		even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
		
		if odd == 'N': # HexN
			cur_Ac.append('1-2')
			cur_SO3.append('1-4')
			cur_SO3.append('1-6')
		else: # HexA
			cur_SO3.append('1-2')
		
		# go through each position in backbone
		for j in range(1, n_pre):
			if (j+1) % 2 == 0: # even
				cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
				
				if even == 'N': # HexN
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-4')
					cur_SO3.append(str(j+1)+'-6')
				else: # HexA
					cur_SO3.append(str(j+1)+'-2')
			else: # odd
				cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
				
				if odd == 'N': # HexN
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-4')
					cur_SO3.append(str(j+1)+'-6')
				else: # HexA
					cur_SO3.append(str(j+1)+'-2')
		
		# add current backbone and modifications
		backbones.append(cur_back)
		all_Ac.append(cur_Ac)
		all_SO3.append(cur_SO3)
else: # working with HS
	poss = ['U','N'] # possible monosaccharides
	
	if NR == '?': # unsure about end monosaccharides
		for i in poss: # we need two backbones
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			odd  = i # odd-numbered monosaccharide
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			# go through each position in backbone
			for j in range(n_pre):
				if (j+1) % 2 == 0: # even
					cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
					
					if even == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-3')
						cur_SO3.append(str(j+1)+'-6')
					else: # HexA
						cur_SO3.append(str(j+1)+'-2')
				else: # odd
					cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
					
					if odd == 'N': # HexN
						cur_Ac.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-2')
						cur_SO3.append(str(j+1)+'-3')
						cur_SO3.append(str(j+1)+'-6')
					else: # HexA
						cur_SO3.append(str(j+1)+'-2')
			
			# add current backbone and modifications
			backbones.append(cur_back)
			all_Ac.append(cur_Ac)
			all_SO3.append(cur_SO3)
	else: # we know end monosaccharides
		cur_back = {} # current backbone
		cur_Ac   = [] # current acetyl modification positions
		cur_SO3  = [] # current sulfate modification positions
		
		# do position 1 first
		cur_back['1'] = NR
		
		# change odd to U if it's a dHexA
		if NR == 'D':
			odd = 'U'
		else:
			odd = NR
		
		even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
		
		if odd == 'N': # HexN
			cur_Ac.append('1-2')
			cur_SO3.append('1-2')
			cur_SO3.append('1-3')
			cur_SO3.append('1-6')
		else: # HexA
			cur_SO3.append('1-2')
		
		# go through each position in backbone
		for j in range(1, n_pre):
			if (j+1) % 2 == 0: # even
				cur_back[str(j+1)] = even # add even-numbered monosaccharide to backbone
				
				if even == 'N': # HexN
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-3')
					cur_SO3.append(str(j+1)+'-6')
				else: # HexA
					cur_SO3.append(str(j+1)+'-2')
			else: # odd
				cur_back[str(j+1)] = odd # add odd-numbered monosaccharide to backbone
				
				if odd == 'N': # HexN
					cur_Ac.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-2')
					cur_SO3.append(str(j+1)+'-3')
					cur_SO3.append(str(j+1)+'-6')
				else: # HexA
					cur_SO3.append(str(j+1)+'-2')
		
		# add current backbone and modifications
		backbones.append(cur_back)
		all_Ac.append(cur_Ac)
		all_SO3.append(cur_SO3)

print "Done!"

################################################
# Step 6: Generate all possible GAG structures #
################################################

print "Generating possible GAG structures...",

# lists for acetyl and sulfate permutations
GAGs = []

# loop through backbones (either one or two)
for idx in range(len(backbones)):
	these_Ac  = choose(all_Ac[idx], pDict['A'])
	these_SO3 = []
	
	if len(these_Ac) == 0:
		qqq = choose(all_SO3[idx], pDict['S'])
		these_SO3.append(qqq)
		
		Ac_dict = None
		
		# check if there is at least one sulfate
		if len(qqq) > 0: # there is at least one sulfate
			# loop through sulfate mods
			for sfp in qqq:
				SO3_dict = {}
				
				# loop through sulfate mods
				for sspot in sfp:
					sunit, sloc = sspot.split('-')
					
					if sunit in SO3_dict:
						SO3_dict[sunit].append(sloc)
					else:
						SO3_dict[sunit] = [sloc]
				
				GAGs.append(generateGAGstring(backbones[idx], Ac_dict, SO3_dict))
		else: # there are no sulfates
			GAGs.append(generateGAGstring(backbones[idx], Ac_dict, None))
	else:
		for a in range(len(these_Ac)):
			qqq = choose(list(set(all_SO3[idx]) - set(these_Ac[a])), pDict['S'])
			these_SO3.append(qqq)
			
			Ac_dict = {}
			
			# loop through Ac mods
			for aspot in these_Ac[a]:
				aunit, aloc = aspot.split('-')
				Ac_dict[aunit] = [aloc]
			
			# check if there is at least one sulfate
			if len(qqq) > 0: # there is at least one sulfate
				# loop through sulfate mods
				for sfp in qqq:
					SO3_dict = {}
					
					for sspot in sfp:
						sunit, sloc = sspot.split('-')
						
						if sunit in SO3_dict:
							SO3_dict[sunit].append(sloc)
						else:
							SO3_dict[sunit] = [sloc]
				
					GAGs.append(generateGAGstring(backbones[idx], Ac_dict, SO3_dict))
			else: # there are no sulfates
				GAGs.append(generateGAGstring(backbones[idx], Ac_dict, None))

print "Done!"

#################################
# Step 7: Build bipartite graph #
#################################

print "Building bipartite graph...",

# initialize graph
gr = nx.Graph()

# edge weights
glyco = 1.
intrl = 0.2
xring = 1.

# loop through all structures
for g in GAGs:
	# set up variables
	cur_dict = gstring2gdict(gClass, g) # convert GAG string to GAG dictionary
	monos    = sorted(cur_dict['backbone'].keys(), key=lambda x: float(x)) # list of monosaccharides
	backbone = cur_dict['backbone']
	actl     = cur_dict['Ac']
	slft     = cur_dict['SO3']
	
	# go through structure
	nn = 1 # variable for fragment length
	while nn < len(monos):
		# loop through fragments of length nn
		for i in range(len(monos)-nn+1):
			j = monos[i:i+nn] # monosaccharides in this fragment
			
			# get the complete monosaccharides
			fml = {'D':0, 'U':0, 'X':0, 'N':0, 'A':0, 'S':0}
			for k in j:
				fml[backbone[k]] += 1 # add current monosaccharide to formula
				
				# count the sulfates
				for p in slft[k]:
					fml['S'] += slft[k][p]
				
				# count the acetyls
				if k in actl:
					for p in actl[k]:
						fml['A'] += actl[k][p]
			
			## glycosidic fragments first
			# loop through possible sulfate losses
			for s in range(s_loss+1):
				fml1 = dict(fml) # copy of original
				fml1['S'] = max(0, fml1['S']-s) # remove sulfate (if possible)
				
				# copy dictionary
				sf = dict2fmla(fml1, 'composition')
				if i == (len(monos)-nn) and fmla is not None:
					sf += '+RE' # add RE to let user know it's reducing end fragment
				
				if sf in f_chgs.keys(): # make sure this fragment has been found in the spectrum
					for cg in f_chgs[sf]:
						if i == 0 or i == (len(monos)-nn): # terminal fragment
							if gr.has_edge(g, fragmap[(sf, cg)]): # edge exists
								gr[g][fragmap[(sf, cg)]]['weight'] = max(gr[g][fragmap[(sf,cg)]]['weight'], glyco ** p1)
							else: # edge does not exist
								gr.add_edge(g, fragmap[(sf, cg)], weight=glyco ** p1)
						else: # internal fragment
							if gr.has_edge(g, fragmap[(sf, cg)]): # edge exists
								gr[g][fragmap[(sf, cg)]]['weight'] = max(gr[g][fragmap[(sf,cg)]]['weight'], intrl ** p1)
							else: # edge does not exist
								gr.add_edge(g, fragmap[(sf, cg)], weight=intrl ** p1)
				
				# loop through neutral losses
				for k in range(2): # water loss
					for l in range(3): # H loss
						sf1 = sf # avoid repeating
						
						# water loss
						if k == 1:
							sf1 += '-H2O'
						
						# H loss(es)
						if l == 1:
							sf1 += '-H'
						elif l == 2:
							sf1 += '-2H'
						
						if sf1 in f_chgs.keys(): # make sure this fragment has been found in the spectrum
							for cg in f_chgs[sf1]:
								if i == 0 or i == (len(monos)-nn): # terminal fragment
									if gr.has_edge(g, fragmap[(sf1, cg)]): # edge exists
										gr[g][fragmap[(sf1, cg)]]['weight'] = max(gr[g][fragmap[(sf1,cg)]]['weight'], glyco ** p1)
									else: # edge does not exist
										gr.add_edge(g, fragmap[(sf1, cg)], weight=glyco ** p1)
								else: # internal fragment
									if gr.has_edge(g, fragmap[(sf1, cg)]): # edge exists
										gr[g][fragmap[(sf1, cg)]]['weight'] = max(gr[g][fragmap[(sf1,cg)]]['weight'], intrl ** p1)
									else: # edge does not exist
										gr.add_edge(g, fragmap[(sf1, cg)], weight=intrl ** p1)
				
				# check if we're done with sulfate losses
				if fml1['S'] == 0:
					break
			
			## crossring fragments second
			if i == 0: # non-reducing end
				addto = monos[nn] # which index to add to
				msch  = backbone[addto] # which monosaccharide to add to
				
				# loop through cross-ring cleavages
				for x in xmodlocs[gClass][msch]['NR']:
					fml1 = dict(fml) # copy of original
					
					# count the sulfates
					for y in xmodlocs[gClass][msch]['NR'][x]['SO3']:
						fml1['S'] += slft[addto][str(y)]
					
					# count the acetyls
					if len(xmodlocs[gClass][msch]['NR'][x]['Ac']) > 0: # acetyl possible
						fml1['A'] += actl[addto]['2']
					
					# sulfate losses
					for s in range(s_loss+1):
						fml2 = dict(fml1) # copy of altered copy
						fml2['S'] = max(0, fml2['S'] - s) # remove sulfate (if possible)
						
						# generate string
						sf = dict2fmla(fml2, 'composition') + '+' + msch + 'NR' + x
						
						if sf in f_chgs.keys(): # make sure this fragment has been found in the spectrum
							for cg in f_chgs[sf]:
								if gr.has_edge(g, fragmap[(sf, cg)]): # edge exists
									gr[g][fragmap[(sf, cg)]]['weight'] = max(gr[g][fragmap[(sf,cg)]]['weight'], xring ** p1)
								else: # edge does not exist
									gr.add_edge(g, fragmap[(sf, cg)], weight=xring ** p1)
						
						# loop through neutral losses
						for l in range(1): # H loss
							sf1 = sf # avoid repeating
							
							# H loss(es)
							if l == 1:
								sf1 += '-H'
							elif l == 2:
								sf1 += '-2H'
							
							if sf1 in f_chgs.keys(): # make sure this fragment has been found in the spectrum
								for cg in f_chgs[sf1]:
									if gr.has_edge(g, fragmap[(sf1, cg)]): # edge exists
										gr[g][fragmap[(sf1, cg)]]['weight'] = max(gr[g][fragmap[(sf1,cg)]]['weight'], xring ** p1)
									else: # edge does not exist
										gr.add_edge(g, fragmap[(sf1, cg)], weight=xring ** p1)
				
						# check if we're done with sulfate losses
						if fml2['S'] == 0:
							break
			elif i == (len(monos)-nn): # reducing end
				addto = monos[i-1] # which index to add to
				msch  = backbone[addto] # which monosaccharide to add to
				
				# loop through cross-ring cleavages
				for x in xmodlocs[gClass][msch]['RE']:
					fml1 = dict(fml) # copy of original
					
					# count the sulfates
					for y in xmodlocs[gClass][msch]['RE'][x]['SO3']:
						fml1['S'] += slft[addto][str(y)]
					
					# count the acetyls
					if len(xmodlocs[gClass][msch]['RE'][x]['Ac']) > 0: # acetyl possible
						fml1['A'] += actl[addto]['2']
					
					# sulfate losses
					for s in range(s_loss+1):
						fml2 = dict(fml1) # copy of altered copy
						fml2['S'] = max(0, fml2['S'] - s) # remove sulfate (if possible)
						
						# generate string
						sf = dict2fmla(fml2, 'composition') + '+' + msch + 'RE' + x
						if fmla is not None:
							sf += '+RE'
						
						if sf in f_chgs.keys(): # make sure this fragment has been found in the spectrum
							for cg in f_chgs[sf]:
								if gr.has_edge(g, fragmap[(sf, cg)]): # edge exists
									gr[g][fragmap[(sf, cg)]]['weight'] = max(gr[g][fragmap[(sf,cg)]]['weight'], xring ** p1)
								else: # edge does not exist
									gr.add_edge(g, fragmap[(sf, cg)], weight=xring ** p1)
						
						# loop through neutral losses
						for l in range(1): # H loss
							sf1 = sf # avoid repeating
							
							# H loss(es)
							if l == 1:
								sf1 += '-H'
							elif l == 2:
								sf1 += '-2H'
							
							if sf1 in f_chgs.keys(): # make sure this fragment has been found in the spectrum
								for cg in f_chgs[sf1]:
									if gr.has_edge(g, fragmap[(sf1, cg)]): # edge exists
										gr[g][fragmap[(sf1, cg)]]['weight'] = max(gr[g][fragmap[(sf1,cg)]]['weight'], xring ** p1)
									else: # edge does not exist
										gr.add_edge(g, fragmap[(sf1, cg)], weight=xring ** p1)
				
						# check if we're done with sulfate losses
						if fml2['S'] == 0:
							break
		
		# increment fragment length
		nn += 1

print "Done!"

########################################################
# Step 8: Get the different rankings for each sequence #
########################################################

print "Calculating enrichment score for each GAG structure..."

# sequences that are in the network
poss_seqs = []
for seq in GAGs:
	if seq in gr.nodes():
		poss_seqs.append(seq)

### BiRank
# get sum of g-scores
g_sum = np.sum(gsc['G'])

# get prior for fragments
prior_frg = {}
for frag in gsc['frags']:
	prior_frg[frag] = gdict[frag] / g_sum

# get prior for sequences
prior_seq = {}
seq_sum   = 0.
for key in poss_seqs:
	this = scoreSeq(key, p3)
	
	prior_seq[key] = this
	seq_sum       += this

for key in prior_seq:
	prior_seq[key] /= seq_sum

# run BiRank
br = BiRank(gr, alpha=.98, beta=.94, f_0=prior_frg, s_0=prior_seq)[0]

# convert to a numpy array for downstream analysis
bsc = []
for key in br:
	bsc.append((key, br[key]))

dty = np.dtype([('seq', 'object'), ('br', 'float')])
bsc = np.array(bsc, dtype=dty)
bsc = np.sort(bsc, order='br')[::-1]

# for debugging
if actual:
	actual = actual.split(':')
	for sqn in actual:
		print '\tBiRank ranking #' + str(int(len(bsc['br']) - rankdata(bsc['br'], 'max')[np.where(bsc['seq'] == sqn)[0][0]] + 1)) + ' to #' + str(int(len(bsc['br']) - rankdata(bsc['br'], 'min')[np.where(bsc['seq'] == sqn)[0][0]] + 1)) + ' out of ' + str(len(bsc['seq']))

print "Done!"

##########################
# Step 9: Output to file #
##########################

# output to file
write_result_to_file(rFile, bsc)

print '\nFinished!'
print time.time() - start_time