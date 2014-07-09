#!/usr/bin/python

""" --- Process scored indels to select the good ones --- 
    Bjorn Pieper. MPIPZ. July 2014    """

import MySQLdb
import csv
import re
import os

score_cutoff = 0.4
begin = 1; end=7

db  = MySQLdb.connect(host="localhost", user="robot", passwd="giveusthetechnology", db="chpoly")
dbc = db.cursor()

with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/fp_v2_out.txt') as ID:
  indels = list(csv.reader(ID, delimiter=","))

for x in range(begin,end):
    res = sorted([indel for indel in indels if (int(indel[3]) == x and float(indel[4]) > score_cutoff)], key= lambda indel: indel[4],reverse=True)
    FILEw = os.path.join("/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/", "results_group_" + str(x) + "_alleles")
    with open(FILEw, "w") as allres: 
        for y in res:
            dbc.execute("select chromosome as chr, position, length as length, count(*) as n, group_concat(concat(Allele) separator ' ') as accessions from variants where chromosome=%s and position=%s group by length;", (y[0],y[1]))  
            snt = dbc.fetchall()
            
            allres.write(','.join(y)); allres.write("\n")
            for z in snt:
                allres.write(''.join( [str(zz) + "," for zz in z[0:len(z)-1]]) + str(z[len(z)-1]) + "\n")
            allres.write("\n")


# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
