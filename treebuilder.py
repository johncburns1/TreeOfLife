from __future__ import division
import os
import sys
import copy
import numpy as np
import nwalign as nw
import time
import matplotlib
matplotlib.use('Agg')
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import scipy
import pylab
from scipy.spatial.distance import pdist
import swalign as sw
import math as mat

if len(sys.argv) < 3:
    print 'Usage: python treebuilder.py [sequences.fa] [keyword1] | [keyword2] > [output.fa]'

MAX_SEQUENCES = 50 # the number of sequences to run on!
name = ''
seq = []
allsequences = {}

def check(name, sequence):
    found = 0
    for keyword in sys.argv[2:]:
        
        if (keyword[0] == '-') != (keyword.replace('-','').lower() in name.lower()):
            found += 1

    if found == len(sys.argv[2:]):
        return True
    return False

with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == '>':
            if check(name, seq):
            	print name + ''.join(seq).strip()
                allsequences[name] = copy.copy(seq)
            seq = []
            name = line
        else:
            seq.append(line)

    if check(name, seq):
        print name + ''.join(seq).strip()
        allsequences[name] = copy.copy(seq)

sys.stderr.write("Total number of sequences matching all keywords:%d" % len(allsequences))

for name in allsequences.keys():
    allsequences[name] = ''.join(map(str.strip,allsequences[name]))

#exit(1)
#filter the sequences down to 1 per species:
unique_species = []
unique_sequences = []
minimum_sequence_length = 1500 #minimum length of 1500, because we are doing global and not local alignment...

for name in allsequences.keys():
    if ';' in name:
        speciesname = name.rpartition(';')[2]
        if len(allsequences[name]) > minimum_sequence_length and speciesname not in unique_species: #minimum length of 1500, because we are doing global and not local alignment...
            unique_species.append(speciesname)
            unique_sequences.append(name)
unique_sequences = sorted(unique_sequences)

print 'Number of unique species: %d'%len(unique_species)



#exit(1)


def scorefunc(a,b,match = 2, mismatch = -1, gap = -3):
	myscore = 0
	    
	if len(a) != len(b):
		print 'Sequences to score must have identical lengths!'
		    	
	for i in xrange(len(a)):
		if a[i] == '-' or b[i] == '-':
			myscore += gap
		if a[i] == b[i]:
			myscore += match
		else:
			if(a[i] == 'C'):
				if(b[i] == 'T' or b[i] == 'U'):
					mismatch = 0
			elif(a[i] == 'T' or a[i] == 'U'):
				if(b[i] == 'C'):
					mismatch = 0
			elif(a[i] == 'A'):
				if(b[i] == 'G'):
					mismatch = 0
			elif(a[i] == 'G'):
				if(b[i] == 'A'):
					mismatch = 0
			mismatch += mismatch
			mismatch = 1
			
	return myscore

t0 = time.clock()


unique_sequences = unique_sequences[0:MAX_SEQUENCES] # only do the first 50 for speed..

similarity_matrix = np.zeros((len(unique_sequences),len(unique_sequences)))

dist_matrix = np.zeros((len(unique_sequences),len(unique_sequences)))

scoring = sw.ScoringMatrix('scoring_matrix.txt')
sw = sw.LocalAlignment(scoring)

match = 2
n = 0

for x,seq1 in enumerate(unique_sequences):
    for y, seq2 in enumerate(unique_sequences):
        
        alignment = nw.global_align(allsequences[seq1],allsequences[seq2])
        
        score = float(nw.score_alignment(alignment[0],alignment[1], gap_open=-5, gap_extend=-2, matrix='scoring_matrix.txt'))
        
        n = float(len(alignment[0])* match)
        if abs(score) > n:
        	score = 0 
        	
        similarity_matrix[x,y] = int(score)
        dist_matrix[x,y] = float(score/n)
       
        	
        
    print '%d/%d Calculated %d alignments in %f seconds'%(x,len(unique_sequences),len(unique_sequences),time.clock()-t0)
    
    t0 = time.clock()

#print 'similarity matrix %s with %s elements in each' %(len(similarity_matrix), len(similarity_matrix[0]))
#print similarity_matrix 

#actually create a nice .tsv similarity file
outfile = open('similarity.tsv', 'w')
outfile.write(' \t'+'\t'.join(unique_sequences)+'\n')
for x,name in enumerate(unique_sequences):
    outfile.write(name +'\t')
    for y in range(len(unique_sequences)):
        outfile.write('\t' + str(similarity_matrix[x,y]))
    outfile.write('\n')

#Draw a large dendrogram:

shortnames = [n.rpartition(';')[2] for n in unique_sequences]

# Compute and plot first dendrogram.
fig = plt.figure(figsize=(16,16))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
Y = sch.linkage(similarity_matrix, method='weighted', metric='euclidean')
Z1 = sch.dendrogram(Y, orientation='right')#, truncate_mode = 'level', p = 50)
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
Y = sch.linkage(similarity_matrix, method='weighted', metric='euclidean')
Z2 = sch.dendrogram(Y)#, truncate_mode = 'level', p = 50)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
similarity_matrix = similarity_matrix[idx1,:]
similarity_matrix = similarity_matrix[:,idx2]
im = axmatrix.matshow(similarity_matrix, aspect='auto', origin='lower', cmap=plt.cm.Spectral)
# axmatrix.set_xticks([])
# axmatrix.set_yticks([])

#moar
axmatrix.set_xticks(range(len(shortnames)))
axmatrix.set_xticklabels([shortnames[i] for i in idx1], minor=False, size =5)
axmatrix.xaxis.set_label_position('bottom')
axmatrix.xaxis.tick_bottom()

plt.xticks(rotation=-90, fontsize=5)

axmatrix.set_yticks(range(len(shortnames)))
axmatrix.set_yticklabels([shortnames[i] for i in idx2], minor=False, size = 5)
axmatrix.yaxis.set_label_position('right')
axmatrix.yaxis.tick_right()

axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

# Plot colorbar.
plt.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('Dendrogram_%s.png'%('_'.join(sys.argv[2:])), dpi = 300)


