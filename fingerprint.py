#!/usr/bin/python

""" --- Select good multi-allelic InDels for fingerprinting --- 
        Bjorn Pieper. MPIPZ. July 2014                          """

import MySQLdb
import csv
import numpy
import re

with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/ma_indel.csv') as ID:
    indels = list(csv.reader(ID, delimiter=","))

db  = MySQLdb.connect(host="localhost", user="robot", passwd="giveusthetechnology", db="chpoly")
dbc = db.cursor()

for loci in range (0,1):
    dbc.execute("SELECT length, count(*) from variants where chromosome=%s and position=%s group by length;", (indels[loci][0], indels[loci][1]))
    locus = dbc.fetchall()
    
    lengths = []
    counts  = []
    for look in locus:
        lengths.append(look[0])
        counts.append(look[1])

    mat_lengths = numpy.zeros([len(lengths),len(lengths)], dtype=int)
    mat_counts  = numpy.zeros([len(lengths),len(lengths)], dtype=int)

    col = 0
    for y in range(0, len(lengths)):
        for x in range(col, len(mat_lengths)):
            if (x == y): 
                mat_lengths[x,y] = lengths[x]
                mat_counts[x,y] = counts[x]
            else:
                mat_lengths[x,y] = lengths[x] - lengths [y]
                mat_counts[x,y] = counts[x] + counts[y]
        col = col + 1
    
    for a,b in numpy.nditer([mat_lengths, mat_counts]):
        print "%d : %d" % (a,b)

print "\n"
print mat_lengths
print "\n"
print mat_counts


with open("/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/test_out.txt", "w") as writer:
    writer.write(u'#'+'\t'.join(str(out) for out in mat_lengths) + '\n')    

"""
for test in range (0,10):
    print indels[test]
    print indels[test][0]
    
    for x,y in numpy.nditer([test3, test_lengths]):
        print "%d:%d" % (x, y)
    
    for look in locus:
        test3.append(look[0])
        lengths = numpy.matrix(locus[0])
        counts  = numpy.matrix(locus[1])
        alles   = numpy.matrix(locus)
    print test3
"""
