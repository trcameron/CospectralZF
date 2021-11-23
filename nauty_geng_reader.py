# Nauty GENG Reader
# Author: Thomas R. Cameron
# Date: 11/23/2021
from sys import stdin
from numpy import array, zeros
from networkx import Graph, draw
from matplotlib import pyplot as plt

###############################################
###             data_to_n                   ###
###############################################
def data_to_n(data):
    """Read initial one-, four- or eight-unit value from digraph6 integer sequence. Return (value, rest of seq.)"""
    if data[0] <= 62:
        return data[0], data[1:]
    elif data[1] <= 62:
        return (data[1] << 12) + (data[2] << 6) + data[3], data[4:]
    else:
        return (data[2] << 30) + (data[3] << 24) + (data[4] << 18) + (data[5] << 12) + (data[6] << 6) + data[7], data[8:]
###############################################
###             bits                        ###
###############################################
def bits(data):
    """Returns sequence of individual bits from 6-bit-per-value list of data values."""
    for d in data:
        for i in [5, 4, 3, 2, 1, 0]:
            yield (d >> i) & 1
###############################################
###             parseData                   ###
###############################################          
def parseData(n,data):
    """Returns stream of pairs b[i], x[i] for sparse6 format."""
    k = 1
    while 1 << k < n:
        k += 1
    
    chunks = iter(data)
    d = None  # partial data word
    dLen = 0  # how many unparsed bits are left in d

    while 1:
        if dLen < 1:
            try:
                d = next(chunks)
            except StopIteration:
                return
            dLen = 6
        dLen -= 1
        b = (d >> dLen) & 1  # grab top remaining bit

        x = d & ((1 << dLen) - 1)  # partially built up value of x
        xLen = dLen  # how many bits included so far in x
        while xLen < k:  # now grab full chunks until we have enough
            try:
                d = next(chunks)
            except StopIteration:
                return
            dLen = 6
            x = (x << 6) + d
            xLen += 6
        x = x >> (xLen - k)  # shift back the extra bits
        dLen = xLen - k
        yield b, x
###############################################
###             sparse6                     ###
###############################################
def sparse6(string):
    # check for sparse6 data type
    if(string[0]!=58):
        raise Exception("sparse6 characters must start with :")
    else:
        string = string[1:]
    # subtract 63 from each bit, check for bits that are over 126
    data = [c - 63 for c in string]
    if any(c > 63 for c in data):
        raise Exception("sprase6 characters must be in range(63, 127)")
    # extract size of graph (n) and remaining bits which hold edge information
    n, data = data_to_n(data)
    # build adjacency matrix
    v = 0 
    a = zeros((n,n),dtype=float)
    
    for b, x in parseData(n,data):
        if b == 1:
            v += 1
        # padding with ones can cause overlarge number here
        if x >= n or v >= n:
            break
        elif x > v:
            v = x
        else:
            a[x,v] = 1; a[v,x] = 1
    return a
###############################################
###             graph6                      ###
###############################################
def graph6(bytes_in):
    # subtract 63 from each bit, check for bits that are over 126
    data = [c - 63 for c in bytes_in]
    if any(c > 63 for c in data):
        raise Exception("graph6 characters must be in range(63, 127)")
    # extract size of graph (n) and remaining bits which hold edge information
    n, data = data_to_n(data)
    # check if we have the correct number of bits
    nd = (n*(n-1)//2 + 5) // 6
    if len(data) != nd:
        raise Exception(f"Expected {nd*6} bits and got {len(data) * 6} in graph6")
    # build adjacency matrix    
    a = zeros((n,n),dtype=float)
    for (i, j), b in zip([(i, j) for j in range(1,n) for i in range(j)], bits(data)):
        if b:
            a[i,j] = 1; a[j,i] = 1
    # return
    return a
###############################################
###             main                        ###
###############################################
def main():
    try:
        # read input stream
        for line in stdin:
            a = graph6(bytes(line.rstrip(),'utf-8'))
            g = Graph(a)
            draw(g,ax=plt.subplot(111))
            plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()