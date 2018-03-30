from __future__ import print_function, unicode_literals
import sys

score_match = -2
score_miss = 1
score_gap = 2


def _traceback(t, r, str1, str2, x, y, s1='', s2=''):
    if x > 0 or y > 0:
        c = t[y][x]
        u = c == (t[y - 1][x] + score_gap)
        l = c == (t[y][x - 1] + score_gap)
        b = str1[x - 1] == str2[y - 1]
        ul = c == (t[y - 1][x - 1] + score_match if b else score_miss)
        if ul:
        	_traceback(t, r, str1, str2,
                       x - 1, y - 1, str1[x - 1] + s1, str2[y - 1] + s2)
        if l:
        	_traceback(t, r, str1, str2,
                       x - 1, y, str1[x - 1] + s1, '-' + s2)
        if u:
        	_traceback(t, r, str1, str2,
                       x, y - 1, '-' + s1, str2[y - 1] + s2)
    else:
	r.append((s1, s2))
        
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
    _traceback(t, results, str1, str2, x, y)
    # Result
    #from pprint import pprint
    #pprint(t)
    print("Score: {0}".format(score))
    for i, result in enumerate(results):
        print("Result {0}:\n{1}\n{2}".format(i + 1, result[0], result[1]))
    
    return result[0] 																																																																																																																									

