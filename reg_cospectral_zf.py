# Regular Cospectral Zero Forcing
# Author: Thomas R. Cameron
# Date: 11/23/2021
from os import popen
from sys import argv
from nauty_geng_reader import graph6
from zero_forcing import zf_num_wf, thro_num, std_forcing_rule, psd_forcing_rule, skw_forcing_rule
from numpy.linalg import eigvalsh, norm
from numpy import array
from networkx import Graph, draw, draw_circular, draw_shell, is_regular
from matplotlib import pyplot as plt
import traceback

TOL = 2**(-52)

###############################################
###             main                        ###
###############################################
def main(arg):
    a1 = array([[0,-1,0,0,0,-1,0,-1,0,0],[-1,0,-1,0,0,0,-1,-1,0,0],[0,-1,0,-1,0,0,-1,0,-1,0],
        [0,0,-1,0,-1,0,0,-1,-1,0],[0,0,0,-1,0,-1,0,0,-1,-1],[-1,0,0,0,-1,0,-1,0,0,-1],
        [0,-1,-1,0,0,-1,0,0,0,-1],[-1,-1,0,-1,0,0,0,0,-1,0],[0,0,-1,-1,-1,0,0,-1,0,0],
        [-1,0,0,0,-1,-1,-1,0,0,0]])
    g1 = Graph(a1)
    zf_num1 = zf_num_wf(g1,std_forcing_rule)
    psd_zf_num1 = zf_num_wf(g1,psd_forcing_rule)
    skw_zf_num1 = zf_num_wf(g1,skw_forcing_rule)
    print(zf_num1)
    print(psd_zf_num1)
    print(skw_zf_num1)
    draw_shell(g1)
    plt.show()
    quit()
    try:
        # graph order
        n = int(arg[0])
        # forcing rule
        forcing_rule = arg[1]
        # open nauty geng
        f = popen('/usr/local/Cellar/nauty/2.7r3/bin/geng %d'%n)
        # storage for graphs, eigenvalues, zero-forcing numbers, and throttling numbers
        graphs = []
        eigvals = []
        lines = []
        # read lines
        line = f.readline()
        while line!= '':
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            # if graph is regular, store graph and eigenvalues
            if(is_regular(g)):
                graphs.append(g)
                eigvals.append(eigvalsh(a))
                lines.append(line)
            
            line = f.readline()
        print("# regular graphs = %d"%len(eigvals))
        # find regular co-spectral digraphs
        count = 0
        for i in range(len(eigvals)):
            for j in range(i+1,len(eigvals)):
                if(norm(eigvals[i]-eigvals[j])<n*TOL*max(1,norm(eigvals[i]))):
                    count += 1
                    # compare zero-forcing numbers
                    zf_numi = zf_num_wf(graphs[i],std_forcing_rule)
                    zf_numj = zf_num_wf(graphs[j],std_forcing_rule)
                    psd_zf_numi = zf_num_wf(graphs[i],psd_forcing_rule)
                    psd_zf_numj = zf_num_wf(graphs[j],psd_forcing_rule)
                    if(zf_numi!=zf_numj and psd_zf_numi != psd_zf_numj):
                        fig,ax = plt.subplots(nrows=1,ncols=2)
                        draw_circular(graphs[i],with_labels=True,ax=ax[0])
                        ax[0].title.set_text('zf = %d, zf+ = %d'%(zf_numi,psd_zf_numi))
                        draw_circular(graphs[j],with_labels=True,ax=ax[1])
                        ax[1].title.set_text('zf = %d, zf+ = %d'%(zf_numj,psd_zf_numj))
                        plt.show()
        print("# regular cospectral graphs = %d"%count)
                        
    except Exception as e:
        print(traceback.format_exc())
if __name__ == '__main__':
    main(argv[1:])