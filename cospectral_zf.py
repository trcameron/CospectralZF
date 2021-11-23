# Cospectral Zero Forcing
# Author: Thomas R. Cameron
# Date: 11/23/2021
from os import popen
from sys import argv
from nauty_geng_reader import graph6
from zero_forcing import zf_num_bf, std_forcing_rule, psd_forcing_rule, skw_forcing_rule
from numpy.linalg import eigvalsh, norm
from networkx import Graph, draw, is_regular
from matplotlib import pyplot as plt
from min_rank import skew_min_rank
import traceback

TOL = 2**(-52)

###############################################
###             main                        ###
###############################################
def main(arg):
    try:
        # graph order
        n = int(arg[0])
        # open nauty geng
        f = popen('/usr/local/Cellar/nauty/27r1/bin/geng %d'%n)
        # storage for graphs, eigenvalues, zero-forcing numbers, and throttling numbers
        graphs = []
        eigvals = []
        lines = []
        # read lines
        line = f.readline()
        while line!= '':
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            
            graphs.append(g)
            eigvals.append(eigvalsh(a))
            lines.append(line)
            
            line = f.readline()
        # find co-spectral digraphs
        count_zf = 1; count_th = 1
        for i in range(len(eigvals)):
            for j in range(i+1,len(eigvals)):
                if(norm(eigvals[i]-eigvals[j])<n*TOL*max(1,norm(eigvals[i]))):
                    # compare skew-zero-forcing numbers
                    skw_zf1 = zf_num_bf(graphs[i],skw_forcing_rule)
                    skw_zf2 = zf_num_bf(graphs[j],skw_forcing_rule)
                    if(len(skw_zf1)!=len(skw_zf2)):
                        # compute max skew nullity
                        skw_mn1 = graphs[i].order() - skew_min_rank(graphs[i])
                        skw_mn2 = graphs[j].order() - skew_min_rank(graphs[j])
                        # check hypothesis
                        if(skw_mn1==len(skw_zf1) and skw_mn2==len(skw_zf2)):
                            print("cs with distinct skew_zf=min_rank")
                            if(is_regular(graphs[i]) or is_regular(graphs[j])):
                                print(lines[i].rstrip()+"\t"+lines[j].rstrip())
                        
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main(argv[1:])