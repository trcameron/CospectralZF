# Zero Forcing
# Author: Thomas R. Cameron
# Date: 11/23/2021
from os import popen
from sys import argv
from networkx import Graph, draw, draw_shell, is_isomorphic, petersen_graph
from collections import deque
from nauty_geng_reader import graph6
from matplotlib import pyplot as plt
from time import time
from numpy import array, inf
from numpy.linalg import matrix_rank
import traceback
        
###############################################
###             has_path                    ###
###############################################
def has_path(g,source,target):
    queue = deque()
    queue.append(source)
    visited = {u : False for u in g.nodes}      # dict comprehension
    visited[source] = True
    while(queue):
        u = queue.popleft()
        for v in g.neighbors(u):
            if(v==target):
                return True
            if(not visited[v]):
                visited[v] = True
                queue.append(v)
    return False
###############################################
###             psd_forcing_rule            ###
###############################################
def psd_forcing_rule(g,colored,u):
    f = set()
    if(not colored[u]):
        return f
    gw = g.subgraph([w for w in g.nodes if not colored[w]])
    s = set([w for w in g.neighbors(u) if not colored[w]])
    for w in s:
        add = True
        for v in (s - set([w])):
            if(has_path(gw,w,v)):
                add = False
                break
        if(add):
            f.add(w)
    return f
###############################################
###             std_forcing_rule            ###
###############################################
def std_forcing_rule(g,colored,u):
    if(not colored[u]):
        return set()
    f = set([w for w in g.neighbors(u) if not colored[w]])
    if(len(f)==1):
        return f
    else:
        return set()
###############################################
###             skw_forcing_rule            ###
###############################################
def skw_forcing_rule(g,colored,u):
    f = set([w for w in g.neighbors(u) if not colored[w]])
    if(len(f)==1):
        return f
    else:
        return set()
###############################################
###             closure                     ###
###############################################
def closure(g,s,f):
    # initialize colored list and stack
    colored = {u : False for u in g.nodes}  # dict comprehension
    stack = deque()
    for u in s:
        colored[u] = True                   # only the vertices in s are colored
    for u in g.nodes:
        forcings = f(g,colored,u)
        if(len(forcings)>0):                # vertices that can do forcing
            stack.append((u,forcings))
    # while there are still vertices that can do forcing
    while(stack):
        pair = stack.pop()                  # (forcing vertex,forcings)
        for v in pair[1]:                   # for each vertex that is forcerd
            colored[v] = True               # color new vertex and update neighbor forcings
            for w in g.neighbors(v):
                forcings = f(g,colored,w)
                if(len(forcings)>0):
                    stack.append((w,forcings))
            forcings = f(g,colored,v)       # update v forcings
            if(len(forcings)>0):
                stack.append((v,forcings))
                      
    # return
    return set([u for u in g.nodes if colored[u]])
###############################################
###     zero-forcing number (brute force)   ###
###############################################
def zf_num_bf(g,f):
    zf_sets = [set()]   # list of possible zero-forcing sets
    # iterate over all possible cardinalities
    for i in range(1,g.order()+1):
        # iterate over all possible zero-forcing sets
        temp = []
        for s in zf_sets:
            # iterate over all nodes not in s:
            for v in (g.nodes - s):
                t = s.union(set([v]))
                if(t not in temp):
                    temp.append(t)
                    if(len(closure(g,temp[-1],f))==g.order()):
                        return temp[-1]
        # update zf_sets
        zf_sets = temp
###############################################
###     zero-forcing number (wavefront)     ###
###############################################
def zf_num_wf(g,f):
    cl_pairs = [(set(),0)] # list of closure pairs (s,r), where s is the closure of some subset of V, 
    # and r is the cardinality of a subset of V whose closure is s. 
    # iterate over all possible cardinalities
    for i in range(1,g.order()+1):
        # iterate over all closure pairs
        for (s,r) in cl_pairs:
            # iterate over all nodes
            for v in g.nodes:
                v_set = set([v])
                v_nbhd = set(g.neighbors(v))
                s_new = closure(g,s.union(v_nbhd,v_set),f)
                r_new = r + int(v not in s) + max(len(v_nbhd-s)-1,0)
                if(r_new<=i and sum([(s_new,j) in cl_pairs for j in range(1,i+1)])==0):
                    cl_pairs.append((s_new,r_new))
                    if(len(s_new)==g.order()):
                        return r_new
###############################################
###             zero-forcing polynomial     ###
###############################################
def zf_poly(g,f):
    # zero-forcing polynomials
    poly = [0 for i in range(g.order()+1)]
    zf_sets = [set(g.nodes)]    # list of zero-forcing sets
    # build zero-forcing sets of each order
    for k in range(g.order(),-1,-1):
        poly[k] = len(zf_sets)
        temp = []
        for s in zf_sets:
            for u in s:
                diff = s - set([u])
                if(len(closure(g,diff,f))==g.order() and diff not in temp):
                    temp.append(diff)
        zf_sets = temp
    # return
    return poly
###############################################
###         time step                       ###
###############################################
def time_step(g,s,f):
    # initialize colored list and stack
    colored = {u : False for u in g.nodes}  # dict comprehension
    stack = deque()
    for u in s:
        colored[u] = True                   # only the vertices in s are colored
    for u in g.nodes:
        forcings = f(g,colored,u)
        if(len(forcings)>0):                # vertices that can do forcing
            stack.append((u,forcings))
    # while there are still vertices that can do forcing
    while(stack):
        pair = stack.pop()                  # (forcing vertex,forcings)
        for v in pair[1]:                   # for each vertex that is forcerd
            colored[v] = True               # color v
    # return
    return set([u for u in g.nodes if colored[u]])
###############################################
###         propogation time helper         ###
###############################################
def prop_time_helper(g,s,f):
    # initialize list of sets constaining colored vertices at each time step
    pt = [s]
    pt.append(time_step(g,pt[-1],f))
    # wihle new vertices are being forced
    while(len(pt[-1])!=len(pt[-2])):
        pt.append(time_step(g,pt[-1],f))
    # if the final coloring is the whole vertex set
    if(len(pt[-1])==g.order()):
        return len(pt)-2
    else:
        return inf
###############################################
###         propogation time                ###
###############################################
def prop_time(g,f):
    # initialize list of zero-forcing sets
    zf_sets = [set(g.nodes)]
    # build zero-forcing sets of each order
    for k in range(g.order(),-1,-1):
        temp = []
        for s in zf_sets:
            for u in s:
                diff = s - set([u])
                if(len(closure(g,diff,f))==g.order() and diff not in temp):
                    temp.append(diff)
        # if there are zero-forcing sets of this order
        if(len(temp)!=0):
            zf_sets = temp
        else:
            break
    # find propogation time for all minimum zero-forcing sets
    pt = []
    for s in zf_sets:
        pt.append(prop_time_helper(g,s,f))
    # return minimum of pt
    return min(pt)
###############################################
###         throttling number helper        ###
###############################################
def thro_num_helper(g,s,f):
    return len(s) + prop_time_helper(g,s,f)
###############################################
###         throttling number               ###
###############################################
def thro_num(g,f):
    # initialize list of zero-forcing sets
    zf_sets = [set(g.nodes)]
    all_zf_sets = [zf_sets[0]]
    # build zero-forcing sets of each order
    for k in range(g.order(),-1,-1):
        temp = []
        for s in zf_sets:
            for u in s:
                diff = s - set([u])
                if(len(closure(g,diff,f))==g.order() and diff not in temp):
                    temp.append(diff)
        # if there are zero-forcing sets of this order
        if(len(temp)!=0):
            zf_sets = temp
            all_zf_sets.extend(temp)
        else:
            break
    # find throttling number for all zero-forcing sets
    th = []
    for s in all_zf_sets:
        th.append(thro_num_helper(g,s,f))
    # return minimum of th
    return min(th)
###############################################
###             main                        ###
###############################################
def main(arg):
    try:
        # graph 1
        line = 'ICR`uiwY_'
        a1 = graph6(bytes(line.rstrip(),'utf-8'))
        g1 = Graph(a1)
        zf1 = len(zf_num_bf(g1,std_forcing_rule))
        # graph 2
        g2 = g1.copy()
        g2.add_node(10)
        for i in range(1,6):
            g2.add_edge(i,10)
        zf2 = len(zf_num_bf(g2,std_forcing_rule))
        # graph 3
        line = 'ICdedhkY_'
        a3 = graph6(bytes(line.rstrip(),'utf-8'))
        g3 = Graph(a3)
        zf3 = len(zf_num_bf(g3,std_forcing_rule))
        # graph 4
        g4 = g3.copy()
        g4.add_node(10)
        for i in range(2,7):
            g4.add_edge(i,10)
        zf4 = len(zf_num_bf(g4,std_forcing_rule))
        # draw graphs
        fig,ax = plt.subplots(nrows=2,ncols=2)
        draw(g1,with_labels=True,ax=ax[0,0])
        ax[0,0].title.set_text('std_zf_number = %d'%zf1)
        draw(g2,with_labels=True,ax=ax[0,1])
        ax[0,1].title.set_text('std_zf_number = %d'%zf2)
        draw(g3,with_labels=True,ax=ax[1,0])
        ax[1,0].title.set_text('std_zf_number = %d'%zf3)
        draw(g4,with_labels=True,ax=ax[1,1])
        ax[1,1].title.set_text('std_zf_number = %d'%zf4)
        plt.show()
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main(argv[1:])