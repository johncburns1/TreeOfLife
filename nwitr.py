from __future__ import print_function, unicode_literals
import sys

score_match = -1
score_miss = 1
score_gap = 2


def traverse(t, str1, str2):
    for i in range(1, len(str1)):
    	for j in range(1, len(str2)):
	    c = t[j][i]
	    u = c == (t[j - 1][i] + score_gap)
	    l = c == (t[j][i - 1] + score_gap)
	    b = str1[i - 1] == str2[j - 1]
	    ul = c == (t[j - 1][i - 1] + score_match if b else score_miss)
		
	    if ul:
        	t[i][j] = t[i - 1][j - 1]
	    if l:
	        t[i][j] = t[i - 1][j]
	    if u:
	        t[i][j] = t[i][j-1]	    
        
def diff(str1, str2):
    
    # Initialization
    l1 = len(str1) + 1
    l2 = len(str2) + 1
    t = [[0 for _ in range(l1)] for _ in range(l2)]
    for i in range(l1):
        t[0][i] = score_gap * i
    for i in range(l2):
        t[i][0] = score_gap * i
    
    # Calc
    for y in range(1, l2):
        for x in range(1, l1):
            b = str1[x - 1] == str2[y - 1]
            t[y][x] = min(
                t[y][x - 1] + score_gap,
                t[y - 1][x] + score_gap,
                t[y - 1][x - 1] + score_match if b else score_miss)
    score = t[l2 - 1][l1 - 1]
    
    # Traceback
    results = []
    traverse(t, str1, str2)
    # Result
    #from pprint import pprint
    #pprint(t)
    print("Score: {0}".format(score))
        
    return score


