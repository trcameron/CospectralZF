# Minimum Rank
# Author: Thomas R. Cameron
# Date: 11/23/2021
from sympy import *
from nauty_geng_reader import graph6
from networkx import Graph, draw, path_graph
from itertools import combinations, tee
import traceback

###############################################
###             all minors                  ###
###############################################
def all_minors(mat,n,k):
    rows = list(combinations(range(n),k))
    cols = rows.copy()
    mlist = [None]*len(rows)*len(cols)
    k = 0
    for i in rows:
        for j in cols:
            mlist[k] = mat[i,j].det()
            k += 1
    return mlist
    
###############################################
###             skew_min_rank               ###
###############################################
def skew_min_rank(g):
    # variable definitions
    n = g.order()
    m = g.size()
    vr = [var('e%d'%i) for i in range(m)]
    vr.append(var('s'))
    # matrix definition
    mat = zeros(n,n)
    k = 0
    for e in g.edges:
        mat[e[0],e[1]] = vr[k]
        mat[e[1],e[0]] = -vr[k]
        k += 1
    # find rank
    for k in range(1,n+1):
        l = all_minors(mat,n,k)
        l.append(prod(vr)-1)
        groeb = groebner(l)
        if(1 not in groeb):
            return k-1
    return n
###############################################
###             main                        ###
###############################################
def main():
    try:
        line = 'F?qb?'
        a = graph6(bytes(line.rstrip(),'utf-8'))
        g = Graph(a)
        
        r = skew_min_rank(g)
        print(r)
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main()