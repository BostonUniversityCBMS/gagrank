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
import numpy as np # for numeric arrays
import scipy as sp # for handling degree matrices
from   scipy.stats import rankdata # for getting actual GAG ranking

print "Done!"

## arrays for info about GAG parts
# initialize dictionaries for weights and formulae
wt = {}
fm = {}

# monoisotopic element weights
wt['monoH'] = 1.00782504
wt['monoC'] = 12.0
wt['monoO'] = 15.99491461956
wt['monoN'] = 14.0030740048
wt['monoS'] = 31.97207100

# monoisotopic compound weights
wt['monoHex']   = 6*wt['monoC'] + 12*wt['monoH'] + 6*wt['monoO']
wt['monoHexA']  = 6*wt['monoC'] + 10*wt['monoH'] + 7*wt['monoO']
wt['monoHexN']  = 6*wt['monoC'] + 13*wt['monoH'] + wt['monoN'] + 5*wt['monoO']
wt['monodHexA'] = 6*wt['monoC'] + 8*wt['monoH'] + 6*wt['monoO']
wt['monoH2O']   = 2*wt['monoH'] + wt['monoO']
wt['monoSO3']   = wt['monoS'] + 3*wt['monoO']
wt['monoAc']    = 2*wt['monoC'] + 2*wt['monoH'] + wt['monoO']

# formulae
fm['Hex']   = {'C':6, 'H':12, 'O':6, 'N':0, 'S':0}
fm['HexA']  = {'C':6, 'H':10, 'O':7, 'N':0, 'S':0}
fm['HexN']  = {'C':6, 'H':13, 'O':5, 'N':1, 'S':0}
fm['dHexA'] = {'C':6, 'H':8, 'O':6, 'N':0, 'S':0}
fm['H2O']   = {'C':0, 'H':2, 'O':1, 'N':0, 'S':0}
fm['SO3']   = {'C':0, 'H':0, 'O':3, 'N':0, 'S':1}
fm['Ac']    = {'C':2, 'H':2, 'O':1, 'N':0, 'S':0}

### modification locations
# initialize dictionary
modlocs       = {}
modlocs['HS'] = {'D':{}, 'U':{}, 'N':{}}
modlocs['CS'] = {'D':{}, 'U':{}, 'N':{}}
modlocs['KS'] = {'X':{}, 'N':{}}

# Acetyl
modlocs['HS']['N'] = {'Ac': [2]}
modlocs['CS']['N'] = {'Ac': [2]}
modlocs['KS']['N'] = {'Ac': [2]}

# Sulfate
modlocs['HS']['D']['SO3'] = [2]
modlocs['HS']['U']['SO3'] = [2]
modlocs['HS']['N']['SO3'] = [2,3,6]
modlocs['CS']['D']['SO3'] = [2]
modlocs['CS']['U']['SO3'] = [2]
modlocs['CS']['N']['SO3'] = [4,6]
modlocs['KS']['X']['SO3'] = [6]
modlocs['KS']['N']['SO3'] = [6]

### cross-ring fragments ###
# initialize dictionaries for weights and formulae
xwt  = {}
xfm  = {}
xmod = {}

### cross-ring formulae ###
## HexA
xfm['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['HexA']['NR'] = {}
xfm['HexA']['RE'] = {}

# add fragments
xfm['HexA']['NR']['0,2'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['NR']['1,3'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['NR']['1,5'] = {'C':5, 'H':8, 'O':5, 'N':0, 'S':0}
xfm['HexA']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['NR']['3,5'] = {'C':3, 'H':4, 'O':3, 'N':0, 'S':0}
xfm['HexA']['NR']['0,3'] = {'C':3, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['HexA']['NR']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexA']['NR']['2,5'] = {'C':4, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['RE']['1,3'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['HexA']['RE']['2,4'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['RE']['3,5'] = {'C':3, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexA']['RE']['1,4'] = {'C':3, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}

## HexA
xfm['dHexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['dHexA']['RE'] = {}

# add fragments
xfm['dHexA']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,3'] = {'C':4, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['dHexA']['RE']['2,4'] = {'C':4, 'H':5, 'O':5, 'N':0, 'S':0}
xfm['dHexA']['RE']['3,5'] = {'C':3, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,4'] = {'C':3, 'H':3, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}

## HexN
xfm['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['HexN']['NR'] = {}
xfm['HexN']['RE'] = {}

# add fragments
xfm['HexN']['NR']['0,2'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['HexN']['NR']['1,3'] = {'C':2, 'H':5, 'O':1, 'N':1, 'S':0}
xfm['HexN']['NR']['1,5'] = {'C':5, 'H':11, 'O':3, 'N':1, 'S':0}
xfm['HexN']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexN']['NR']['3,5'] = {'C':3, 'H':6, 'O':2, 'N':0, 'S':0}
xfm['HexN']['NR']['0,3'] = {'C':3, 'H':7, 'O':2, 'N':1, 'S':0}
xfm['HexN']['NR']['1,4'] = {'C':3, 'H':7, 'O':2, 'N':1, 'S':0}
xfm['HexN']['NR']['2,5'] = {'C':4, 'H':5, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['0,2'] = {'C':2, 'H':5, 'O':1, 'N':1, 'S':0}
xfm['HexN']['RE']['1,3'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['HexN']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['HexN']['RE']['2,4'] = {'C':4, 'H':9, 'O':3, 'N':1, 'S':0}
xfm['HexN']['RE']['3,5'] = {'C':3, 'H':7, 'O':3, 'N':1, 'S':0}
xfm['HexN']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['2,5'] = {'C':2, 'H':5, 'O':2, 'N':1, 'S':0}

## Hex
xfm['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['Hex']['NR'] = {}
xfm['Hex']['RE'] = {}

# add fragments
xfm['Hex']['NR']['0,2'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['NR']['1,3'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['NR']['1,5'] = {'C':5, 'H':10, 'O':4, 'N':0, 'S':0}
xfm['Hex']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['NR']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['Hex']['NR']['2,5'] = {'C':4, 'H':8, 'O':3, 'N':0, 'S':0}
xfm['Hex']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['RE']['1,3'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['Hex']['RE']['2,4'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['RE']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['Hex']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}



### cross-ring weights ###
## HexA
xwt['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['HexA']['NR'] = {}
xwt['HexA']['RE'] = {}

# add weights
xwt['HexA']['NR']['0,2'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['NR']['1,3'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['NR']['1,5'] = 5*wt['monoC'] + 8*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['NR']['3,5'] = 3*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['NR']['0,3'] = 3*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['NR']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['NR']['2,5'] = 4*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['RE']['1,3'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['RE']['1,5'] = 1*wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['RE']['2,4'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['RE']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['RE']['1,4'] = 3*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']

## dHexA
xwt['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xwt['dHexA']['RE']  = {}

# add weights
xwt['dHexA']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['dHexA']['RE']['1,3'] = 4*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['1,5'] = wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['dHexA']['RE']['2,4'] = 4*wt['monoC'] + 5*wt['monoH'] + 5*wt['monoO']
xwt['dHexA']['RE']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['dHexA']['RE']['1,4'] = 3*wt['monoC'] + 3*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']

## HexN
xwt['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['HexN']['NR'] = {}
xwt['HexN']['RE']  = {}

# add weights
xwt['HexN']['NR']['0,2'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['HexN']['NR']['1,3'] = 2*wt['monoC'] + 5*wt['monoH'] + wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['1,5'] = 5*wt['monoC'] + 11*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['NR']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['NR']['0,3'] = 3*wt['monoC'] + 7*wt['monoH'] + 2*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['1,4'] = 3*wt['monoC'] + 7*wt['monoH'] + 2*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['2,5'] = 4*wt['monoC'] + 5*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['0,2'] = 2*wt['monoC'] + 5*wt['monoH'] + 1*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['1,3'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['HexN']['RE']['1,5'] = 1*wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['RE']['2,4'] = 4*wt['monoC'] + 9*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['3,5'] = 3*wt['monoC'] + 7*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['2,5'] = 2*wt['monoC'] + 5*wt['monoH'] + 2*wt['monoO'] + wt['monoN']

## Hex
xwt['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['Hex']['NR'] = {}
xwt['Hex']['RE'] = {}

# add fragments
xwt['Hex']['NR']['0,2'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['NR']['1,3'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['NR']['1,5'] = 5*wt['monoC'] + 10*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['NR']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['NR']['2,5'] = 4*wt['monoC'] + 8*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['RE']['1,3'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['RE']['1,5'] = wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['RE']['2,4'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['RE']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']



### cross-ring sulfation/acetylation/adduct possibilities, and COOH places
## HS/Heparin
xmod['HS'] = {}

# HexA
xmod['HS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexA']['NR'] = {}
xmod['HS']['HexA']['RE'] = {}

# add modification possibilities
xmod['HS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# dHexA
xmod['HS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['HS']['dHexA']['RE']  = {}

# add weights
xmod['HS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# HexN
xmod['HS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexN']['NR'] = {}
xmod['HS']['HexN']['RE'] = {}

# add modification possibilities
xmod['HS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['1,5'] = {'SO3':3, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['1,4'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['0,2'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['2,4'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['3,5'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['0,3'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['2,5'] = {'SO3':1, 'Ac':1, 'COOH':0}


## CS/DS
xmod['CS'] = {}

# HexA
xmod['CS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexA']['NR'] = {}
xmod['CS']['HexA']['RE'] = {}

# add modification possibilities
xmod['CS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# dHexA
xmod['CS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['CS']['dHexA']['RE']  = {}

# add weights
xmod['CS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# HexN
xmod['CS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexN']['NR'] = {}
xmod['CS']['HexN']['RE'] = {}

# add modification possibilities
xmod['CS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,3'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,5'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,3'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['2,5'] = {'SO3':0, 'Ac':1, 'COOH':0}


## HS
xmod['KS'] = {}

# Hex
xmod['KS']['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['Hex']['NR'] = {}
xmod['KS']['Hex']['RE'] = {}

# add modification possibilities
xmod['KS']['Hex']['NR']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,3'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':0}

# HexN
xmod['KS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['HexN']['NR'] = {}
xmod['KS']['HexN']['RE'] = {}

# add modification possibilities
xmod['KS']['HexN']['NR']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['1,5'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['RE']['3,5'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':0, 'Ac':1, 'COOH':0}



### cross-ring modification locations
# initialize dictionaries
xmodlocs = {}
xmodlocs['HS'] = {'D':{'NR':{}, 'RE':{}}, 'U':{'NR':{}, 'RE':{}}, 'N':{'NR':{}, 'RE':{}}}
xmodlocs['CS'] = {'D':{'NR':{}, 'RE':{}}, 'U':{'NR':{}, 'RE':{}}, 'N':{'NR':{}, 'RE':{}}}
xmodlocs['KS'] = {'X':{'NR':{}, 'RE':{}}, 'N':{'NR':{}, 'RE':{}}}

## HS
# dHexA
xmodlocs['HS']['D']['NR']['0,2'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['NR']['1,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['NR']['3,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['NR']['0,3'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['NR']['1,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['NR']['2,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['RE']['0,2'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['RE']['2,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['RE']['3,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['RE']['0,3'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['D']['RE']['1,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['D']['RE']['2,5'] = {'SO3':[2], 'Ac':[]}

# HexA
xmodlocs['HS']['U']['NR']['0,2'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['NR']['1,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['NR']['3,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['NR']['0,3'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['NR']['1,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['NR']['2,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['RE']['0,2'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['RE']['2,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['RE']['3,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['RE']['0,3'] = {'SO3':[2], 'Ac':[]}
xmodlocs['HS']['U']['RE']['1,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['U']['RE']['2,5'] = {'SO3':[2], 'Ac':[]}

# HexN
xmodlocs['HS']['N']['NR']['0,2'] = {'SO3':[3,6], 'Ac':[]}
xmodlocs['HS']['N']['NR']['1,5'] = {'SO3':[2,3,6], 'Ac':[2]}
xmodlocs['HS']['N']['NR']['2,4'] = {'SO3':[3], 'Ac':[]}
xmodlocs['HS']['N']['NR']['3,5'] = {'SO3':[6], 'Ac':[]}
xmodlocs['HS']['N']['NR']['0,3'] = {'SO3':[6], 'Ac':[]}
xmodlocs['HS']['N']['NR']['1,4'] = {'SO3':[2,3], 'Ac':[2]}
xmodlocs['HS']['N']['NR']['2,5'] = {'SO3':[3,6], 'Ac':[]}
xmodlocs['HS']['N']['RE']['0,2'] = {'SO3':[2], 'Ac':[2]}
xmodlocs['HS']['N']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['HS']['N']['RE']['2,4'] = {'SO3':[2,6], 'Ac':[2]}
xmodlocs['HS']['N']['RE']['3,5'] = {'SO3':[2,3], 'Ac':[2]}
xmodlocs['HS']['N']['RE']['0,3'] = {'SO3':[2,3], 'Ac':[2]}
xmodlocs['HS']['N']['RE']['1,4'] = {'SO3':[6], 'Ac':[]}
xmodlocs['HS']['N']['RE']['2,5'] = {'SO3':[2], 'Ac':[2]}

## CS/DS
# dHexA
xmodlocs['CS']['D']['NR']['0,2'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['NR']['1,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['NR']['3,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['NR']['0,3'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['NR']['1,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['NR']['2,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['RE']['0,2'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['RE']['2,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['RE']['3,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['RE']['0,3'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['D']['RE']['1,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['D']['RE']['2,5'] = {'SO3':[2], 'Ac':[]}

# HexA
xmodlocs['CS']['U']['NR']['0,2'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['NR']['1,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['NR']['3,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['NR']['0,3'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['NR']['1,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['NR']['2,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['RE']['0,2'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['RE']['2,4'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['RE']['3,5'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['RE']['0,3'] = {'SO3':[2], 'Ac':[]}
xmodlocs['CS']['U']['RE']['1,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['U']['RE']['2,5'] = {'SO3':[2], 'Ac':[]}

# HexN
xmodlocs['CS']['N']['NR']['0,2'] = {'SO3':[4,6], 'Ac':[]}
xmodlocs['CS']['N']['NR']['1,3'] = {'SO3':[], 'Ac':[2]}
xmodlocs['CS']['N']['NR']['1,5'] = {'SO3':[4,6], 'Ac':[2]}
xmodlocs['CS']['N']['NR']['2,4'] = {'SO3':[4], 'Ac':[]}
xmodlocs['CS']['N']['NR']['1,4'] = {'SO3':[4], 'Ac':[2]}
xmodlocs['CS']['N']['NR']['2,5'] = {'SO3':[4,6], 'Ac':[]}
xmodlocs['CS']['N']['RE']['0,2'] = {'SO3':[], 'Ac':[2]}
xmodlocs['CS']['N']['RE']['1,3'] = {'SO3':[4,6], 'Ac':[]}
xmodlocs['CS']['N']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['CS']['N']['RE']['2,4'] = {'SO3':[6], 'Ac':[2]}
xmodlocs['CS']['N']['RE']['1,4'] = {'SO3':[6], 'Ac':[]}
xmodlocs['CS']['N']['RE']['2,5'] = {'SO3':[], 'Ac':[2]}

## KS
# HexA
xmodlocs['KS']['X']['NR']['0,2'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['NR']['1,3'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['X']['NR']['1,5'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['X']['NR']['1,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['X']['NR']['2,5'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['RE']['0,2'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['X']['RE']['1,3'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['X']['RE']['2,4'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['RE']['1,4'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['X']['RE']['2,5'] = {'SO3':[], 'Ac':[]}

# HexN
xmodlocs['KS']['N']['NR']['0,2'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['N']['NR']['1,5'] = {'SO3':[6], 'Ac':[2]}
xmodlocs['KS']['N']['NR']['2,4'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['N']['NR']['3,5'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['N']['NR']['0,3'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['N']['NR']['1,4'] = {'SO3':[], 'Ac':[2]}
xmodlocs['KS']['N']['NR']['2,5'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['N']['RE']['0,2'] = {'SO3':[], 'Ac':[2]}
xmodlocs['KS']['N']['RE']['1,5'] = {'SO3':[], 'Ac':[]}
xmodlocs['KS']['N']['RE']['2,4'] = {'SO3':[6], 'Ac':[2]}
xmodlocs['KS']['N']['RE']['3,5'] = {'SO3':[], 'Ac':[2]}
xmodlocs['KS']['N']['RE']['0,3'] = {'SO3':[], 'Ac':[2]}
xmodlocs['KS']['N']['RE']['1,4'] = {'SO3':[6], 'Ac':[]}
xmodlocs['KS']['N']['RE']['2,5'] = {'SO3':[], 'Ac':[2]}

# scaling factors
p1 = 5.4 # for scaling edge width
p2 = 5.1 # for scaling fragments' prior probabilities
p3 = 0.4 # for scaling sequences' prior probabilities

# variable for debugging
debug = False

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

# function for reading gagfinder output
def read_gf_results(res_file):		
	# load G scores
	try:
		res = open(res_file, 'r')
	except:
		print "\nIncorrect or missing file. Please try again."
		sys.exit()
	
	res.readline() # skip first line
	
	# dictionaries
	gd = {} # store G scores for each nonspecific fragment
	fm = {} # map nonspecific fragments to specific fragment/charge state pair
	mf = {} # map specific fragment/charge state pair to nonspecific fragments
	fc = {} # keys are unique fragments and values are each charge state
	af = [] # list of all fragments (similar to GAGs
	
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
		
		# add G score to dictionary
		gd['frag'+str(line_ct)] = G
		
		# add fragment to allfrags list
		af.append('frag'+str(line_ct))
		
		# go through each fragment
		for frg in fgs:
			# add to fragment map
			fm[(frg, chg)] = 'frag'+str(line_ct)
			
			# add to map fragment
			if 'frag'+str(line_ct) not in mf:
				mf['frag'+str(line_ct)] = [(frg, chg)]
			else:
				mf['frag'+str(line_ct)].append((frg, chg))
			
			if frg in fc: # this fragment hasn't been added yet
				fc[frg].append(chg)
			else: # this fragment has been added
				fc[frg] = [chg]
		
		line_ct += 1
	
	# convert to a numpy array
	gs = []
	for key in gd:
		gs.append((key, gd[key]))
	
	dty = np.dtype([('frags', 'object'), ('G', 'float')])
	gs  = np.array(gs, dtype=dty)
	gs  = np.sort(gs, order='G')[::-1]
	
	# return dicts etc.
	return [gd, fm, mf, fc, af, gs]

# function for getting precursor composition
def get_pre(cl, pm, pz, rw, crsr):
	# calculate precursor mass
	pre_mass = (pm*abs(pz)) - (pz*wt['monoH'])
	test_mass = pre_mass - rw
	
	# get precursor info
	crsr.execute('''SELECT   p.value
   		            FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
       		        WHERE    p.id = cpm.pId 
           		    AND      f.id = p.fmId
	                AND      cpm.cId = ?
   		            ORDER BY ABS(f.monoMass - ?) ASC
       		        LIMIT 1;''', (cl, test_mass))
	row = crsr.fetchone()
	
	# return
	return row[0]

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
			if pd['N'] > pd['X']:
				nonred = 'N'
				redend = 'N'
			elif pd['N'] < pd['X']:
				nonred = 'X'
				redend = 'X'
			else:
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

# function for getting backbone and modification info
def get_bb(cl, nonred, lng):
	# variables to keep possibilities
	bb = []
	ac = []
	so = []
	
	# generate backbones and modification possibile locations
	if cl == 'KS': # working with KS
		poss = ['X','N'] # possible monosaccharides
		
		if nonred == '?': # unsure about end monosaccharides
			for i in poss: # we need two backbones
				cur_back = {} # current backbone
				cur_Ac   = [] # current acetyl modification positions
				cur_SO3  = [] # current sulfate modification positions
				
				odd  = i # odd-numbered monosaccharide
				even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
				
				# go through each position in backbone
				for j in range(lng):
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
				bb.append(cur_back)
				ac.append(cur_Ac)
				so.append(cur_SO3)
		else: # we know end monosaccharides
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			odd  = nonred
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			# go through each position in backbone
			for j in range(lng):
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
			bb.append(cur_back)
			ac.append(cur_Ac)
			so.append(cur_SO3)
	elif cl == 'CS': # working with CS
		poss = ['U','N'] # possible monosaccharides
		
		if nonred == '?': # unsure about end monosaccharides
			for i in poss: # we need two backbones
				cur_back = {} # current backbone
				cur_Ac   = [] # current acetyl modification positions
				cur_SO3  = [] # current sulfate modification positions
				
				odd  = i # odd-numbered monosaccharide
				even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
				
				# go through each position in backbone
				for j in range(lng):
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
				bb.append(cur_back)
				ac.append(cur_Ac)
				so.append(cur_SO3)
		else: # we know end monosaccharides
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			# do position 1 first
			cur_back['1'] = nonred
			
			# change odd to U if it's a dHexA
			if nonred == 'D':
				odd = 'U'
			else:
				odd = nonred
			
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			if odd == 'N': # HexN
				cur_Ac.append('1-2')
				cur_SO3.append('1-4')
				cur_SO3.append('1-6')
			else: # HexA
				cur_SO3.append('1-2')
			
			# go through each position in backbone
			for j in range(1, lng):
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
			bb.append(cur_back)
			ac.append(cur_Ac)
			so.append(cur_SO3)
	else: # working with HS
		poss = ['U','N'] # possible monosaccharides
		
		if nonred == '?': # unsure about end monosaccharides
			for i in poss: # we need two backbones
				cur_back = {} # current backbone
				cur_Ac   = [] # current acetyl modification positions
				cur_SO3  = [] # current sulfate modification positions
				
				odd  = i # odd-numbered monosaccharide
				even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
				
				# go through each position in backbone
				for j in range(lng):
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
				bb.append(cur_back)
				ac.append(cur_Ac)
				so.append(cur_SO3)
		else: # we know end monosaccharides
			cur_back = {} # current backbone
			cur_Ac   = [] # current acetyl modification positions
			cur_SO3  = [] # current sulfate modification positions
			
			# do position 1 first
			cur_back['1'] = nonred
			
			# change odd to U if it's a dHexA
			if nonred == 'D':
				odd = 'U'
			else:
				odd = nonred
			
			even = list(set(poss) - set(odd))[0] # even-numbered monosaccharide
			
			if odd == 'N': # HexN
				cur_Ac.append('1-2')
				cur_SO3.append('1-2')
				cur_SO3.append('1-3')
				cur_SO3.append('1-6')
			else: # HexA
				cur_SO3.append('1-2')
			
			# go through each position in backbone
			for j in range(1, lng):
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
			bb.append(cur_back)
			ac.append(cur_Ac)
			so.append(cur_SO3)
	
	# return
	return [bb, ac, so]

# function for generating GAG strings
def get_GAGs(bb, ac, so, pd):
	# lists for acetyl and sulfate permutations
	g = []
	
	# loop through backbones (either one or two)
	for idx in range(len(bb)):
		these_Ac  = choose(ac[idx], pd['A'])
		these_SO3 = []
		
		if len(these_Ac) == 0:
			qqq = choose(so[idx], pd['S'])
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
					
					g.append(generateGAGstring(bb[idx], Ac_dict, SO3_dict))
			else: # there are no sulfates
				g.append(generateGAGstring(bb[idx], Ac_dict, None))
		else:
			for a in range(len(these_Ac)):
				qqq = choose(list(set(so[idx]) - set(these_Ac[a])), pd['S'])
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
					
						g.append(generateGAGstring(bb[idx], Ac_dict, SO3_dict))
				else: # there are no sulfates
					g.append(generateGAGstring(bb[idx], Ac_dict, None))
	
	# return
	return g

# function for building bipartite graph
def get_graph(gg, cl, loss, ref, fc, fm):
	# initialize graph
	gph = nx.Graph()
	
	# edge weights
	glyco = 1.
	intrl = 0.2
	xring = 1.
	
	# loop through all structures
	for g in gg:
		# set up variables
		cur_dict = gstring2gdict(cl, g) # convert GAG string to GAG dictionary
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
				for s in range(loss+1):
					fml1 = dict(fml) # copy of original
					fml1['S'] = max(0, fml1['S']-s) # remove sulfate (if possible)
					
					# copy dictionary
					sf = dict2fmla(fml1, 'composition')
					if i == (len(monos)-nn) and ref is not None:
						sf += '+RE' # add RE to let user know it's reducing end fragment
					
					if sf in fc.keys(): # make sure this fragment has been found in the spectrum
						for cg in fc[sf]:
							if i == 0 or i == (len(monos)-nn): # terminal fragment
								if gph.has_edge(g, fm[(sf, cg)]): # edge exists
									gph[g][fm[(sf, cg)]]['weight'] = max(gph[g][fm[(sf,cg)]]['weight'], glyco ** p1)
								else: # edge does not exist
									gph.add_edge(g, fm[(sf, cg)], weight=glyco ** p1)
							else: # internal fragment
								if gph.has_edge(g, fm[(sf, cg)]): # edge exists
									gph[g][fm[(sf, cg)]]['weight'] = max(gph[g][fm[(sf,cg)]]['weight'], intrl ** p1)
								else: # edge does not exist
									gph.add_edge(g, fm[(sf, cg)], weight=intrl ** p1)
					
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
							
							if sf1 in fc.keys(): # make sure this fragment has been found in the spectrum
								for cg in fc[sf1]:
									if i == 0 or i == (len(monos)-nn): # terminal fragment
										if gph.has_edge(g, fm[(sf1, cg)]): # edge exists
											gph[g][fm[(sf1, cg)]]['weight'] = max(gph[g][fm[(sf1,cg)]]['weight'], glyco ** p1)
										else: # edge does not exist
											gph.add_edge(g, fm[(sf1, cg)], weight=glyco ** p1)
									else: # internal fragment
										if gph.has_edge(g, fm[(sf1, cg)]): # edge exists
											gph[g][fm[(sf1, cg)]]['weight'] = max(gph[g][fm[(sf1,cg)]]['weight'], intrl ** p1)
										else: # edge does not exist
											gph.add_edge(g, fm[(sf1, cg)], weight=intrl ** p1)
					
					# check if we're done with sulfate losses
					if fml1['S'] == 0:
						break
				
				## crossring fragments second
				if i == 0: # non-reducing end
					addto = monos[nn] # which index to add to
					msch  = backbone[addto] # which monosaccharide to add to
					
					# loop through cross-ring cleavages
					for x in xmodlocs[cl][msch]['NR']:
						fml1 = dict(fml) # copy of original
						
						# count the sulfates
						for y in xmodlocs[cl][msch]['NR'][x]['SO3']:
							fml1['S'] += slft[addto][str(y)]
						
						# count the acetyls
						if len(xmodlocs[cl][msch]['NR'][x]['Ac']) > 0: # acetyl possible
							fml1['A'] += actl[addto]['2']
						
						# sulfate losses
						for s in range(loss+1):
							fml2 = dict(fml1) # copy of altered copy
							fml2['S'] = max(0, fml2['S'] - s) # remove sulfate (if possible)
							
							# generate string
							sf = dict2fmla(fml2, 'composition') + '+' + msch + 'NR' + x
							
							if sf in fc.keys(): # make sure this fragment has been found in the spectrum
								for cg in fc[sf]:
									if gph.has_edge(g, fm[(sf, cg)]): # edge exists
										gph[g][fm[(sf, cg)]]['weight'] = max(gph[g][fm[(sf,cg)]]['weight'], xring ** p1)
									else: # edge does not exist
										gph.add_edge(g, fm[(sf, cg)], weight=xring ** p1)
							
							# loop through neutral losses
							for l in range(1): # H loss
								sf1 = sf # avoid repeating
								
								# H loss(es)
								if l == 1:
									sf1 += '-H'
								elif l == 2:
									sf1 += '-2H'
								
								if sf1 in fc.keys(): # make sure this fragment has been found in the spectrum
									for cg in fc[sf1]:
										if gph.has_edge(g, fm[(sf1, cg)]): # edge exists
											gph[g][fm[(sf1, cg)]]['weight'] = max(gph[g][fm[(sf1,cg)]]['weight'], xring ** p1)
										else: # edge does not exist
											gph.add_edge(g, fm[(sf1, cg)], weight=xring ** p1)
					
							# check if we're done with sulfate losses
							if fml2['S'] == 0:
								break
				elif i == (len(monos)-nn): # reducing end
					addto = monos[i-1] # which index to add to
					msch  = backbone[addto] # which monosaccharide to add to
					
					# loop through cross-ring cleavages
					for x in xmodlocs[cl][msch]['RE']:
						fml1 = dict(fml) # copy of original
						
						# count the sulfates
						for y in xmodlocs[cl][msch]['RE'][x]['SO3']:
							fml1['S'] += slft[addto][str(y)]
						
						# count the acetyls
						if len(xmodlocs[cl][msch]['RE'][x]['Ac']) > 0: # acetyl possible
							fml1['A'] += actl[addto]['2']
						
						# sulfate losses
						for s in range(loss+1):
							fml2 = dict(fml1) # copy of altered copy
							fml2['S'] = max(0, fml2['S'] - s) # remove sulfate (if possible)
							
							# generate string
							sf = dict2fmla(fml2, 'composition') + '+' + msch + 'RE' + x
							if ref is not None:
								sf += '+RE'
							
							if sf in fc.keys(): # make sure this fragment has been found in the spectrum
								for cg in fc[sf]:
									if gph.has_edge(g, fm[(sf, cg)]): # edge exists
										gph[g][fm[(sf, cg)]]['weight'] = max(gph[g][fm[(sf,cg)]]['weight'], xring ** p1)
									else: # edge does not exist
										gph.add_edge(g, fm[(sf, cg)], weight=xring ** p1)
							
							# loop through neutral losses
							for l in range(1): # H loss
								sf1 = sf # avoid repeating
								
								# H loss(es)
								if l == 1:
									sf1 += '-H'
								elif l == 2:
									sf1 += '-2H'
								
								if sf1 in fc.keys(): # make sure this fragment has been found in the spectrum
									for cg in fc[sf1]:
										if gph.has_edge(g, fm[(sf1, cg)]): # edge exists
											gph[g][fm[(sf1, cg)]]['weight'] = max(gph[g][fm[(sf1,cg)]]['weight'], xring ** p1)
										else: # edge does not exist
											gph.add_edge(g, fm[(sf1, cg)], weight=xring ** p1)
					
							# check if we're done with sulfate losses
							if fml2['S'] == 0:
								break
			
			# increment fragment length
			nn += 1
	
	return gph

# function for ranking the nodes in the graph
def rank_nodes(gg, gph, g_scores, gd):
	# sequences that are in the network
	poss_seqs = []
	for seq in gg:
		if seq in gph.nodes():
			poss_seqs.append(seq)
	
	### BiRank
	# get sum of g-scores
	g_sum = np.sum(g_scores['G'])
	
	# get prior for fragments
	prior_frg = {}
	for frag in g_scores['frags']:
		prior_frg[frag] = gd[frag] / g_sum
	
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
	br = BiRank(gph, alpha=.98, beta=.94, f_0=prior_frg, s_0=prior_seq)[0]
	
	# convert to a numpy array for downstream analysis
	b_score = []
	for key in br:
		b_score.append((key, br[key]))
	
	dty     = np.dtype([('seq', 'object'), ('br', 'float')])
	b_score = np.array(b_score, dtype=dty)
	b_score = np.sort(b_score, order='br')[::-1]
	
	return b_score

# function for running the guts of GAGrank
def rank_gags(gf, gc, cn, rf, mz, z, sl, a=None, db_path='../lib/GAGfragDB.db'):	
	# get values ready for reducing end derivatization and reagent
	df    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
	atoms = ['C','H','O','N','S']
	
	# parse the reducing end derivatization formula
	if rf:
		parts = re.findall(r'([A-Z][a-z]*)(\d*)', rf.upper()) # split formula by symbol
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
	
	# print the reducing end derivatization back out to the user
	if debug:
		formula = ''
		for key in df:
			val = df[key]
			if val > 0:
				formula += key
				if val > 1:
					formula += str(val)
		print "atoms in reducing end derivatization: %s" % (formula)
	
	###########################################################
	# Step 2: Load GAGfinder results and connect to GAGfragDB #
	###########################################################
	
	print "Loading GAGfinder results file...",
	
	gdict, fragmap, mapfrag, f_chgs, allfrags, gsc = read_gf_results(gf)
	
	# connect to GAGfragDB
	conn = sq.connect(db_path)
	c    = conn.cursor()
	
	print "Done!"
	
	######################################
	# Step 3: Find precursor composition #
	######################################
	
	print "Determining precursor composition...",
	
	pComp = get_pre(cn, mz, z, dw, c)
	
	print "Done!"
	
	#########################################################
	# Step 4: Get reducing end/non-reducing end information #
	#########################################################
	
	print "Determining reducing end/non-reducing end information...",
	
	# convert composition into a dictionary
	pDict         = fmla2dict(pComp, 'composition')
	NR, RE, n_pre = get_ends(pDict, gc)
	
	print "Done!"
	
	#########################################
	# Step 5: Set up possible GAG backbones #
	#########################################
	
	print "Generating possible backbones and modification positions...",
	
	backbones, all_Ac, all_SO3 = get_bb(gc, NR, n_pre)
	
	print "Done!"
	
	################################################
	# Step 6: Generate all possible GAG structures #
	################################################
	
	print "Generating possible GAG structures...",
	
	GAGs = get_GAGs(backbones, all_Ac, all_SO3, pDict)
	
	print "Done!"
	
	#################################
	# Step 7: Build bipartite graph #
	#################################
	
	print "Building bipartite graph...",
	
	gr = get_graph(GAGs, gc, sl, rf, f_chgs, fragmap)
	
	print "Done!"
	
	########################################################
	# Step 8: Get the different rankings for each sequence #
	########################################################
	
	print "Calculating enrichment score for each GAG structure...",
	
	bsc = rank_nodes(GAGs, gr, gsc, gdict)
	
	print "Done!"
	
	# return
	return bsc

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

# main function
def main():
	################################
	# Step 1: check user arguments #
	################################
	
	# initiate parser
	parser = argparse.ArgumentParser(description='Find isotopic clusters in GAG tandem mass spectra.')
	
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
	
	# run the guts of GAGrank
	result = rank_gags(rFile, gClass, cNum, fmla, pre_mz, pre_z, s_loss, actual)
	
	# for debugging
	if actual:
		actual = actual.split(':')
		for sqn in actual:
			print '\tBiRank ranking #' + str(int(len(result['br']) - rankdata(result['br'], 'max')[np.where(result['seq'] == sqn)[0][0]] + 1)) + ' to #' + str(int(len(result['br']) - rankdata(result['br'], 'min')[np.where(result['seq'] == sqn)[0][0]] + 1)) + ' out of ' + str(len(result['seq']))
	
	#########################
	# Step 9: write to file #
	#########################
	
	# write result to file
	write_result_to_file(rFile, result)
	
	print "Finished!"
	print time.time() - start_time

# run main
if __name__ == '__main__':
	main()