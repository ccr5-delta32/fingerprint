#!/usr/bin/python

""" --- Select good multi-allelic InDels for fingerprinting --- 
    Bjorn Pieper. MPIPZ. July 2014    """

import MySQLdb
import csv
import numpy as np
import re

with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/ma_indel.csv') as ID:
  indels = list(csv.reader(ID, delimiter=","))

db  = MySQLdb.connect(host="localhost", user="robot", passwd="giveusthetechnology", db="chpoly")
dbc = db.cursor()
"""
for loci in range (0,1):
  dbc.execute("SELECT length, count(*) from variants where chromosome=%s and position=%s group by length;", (indels[loci][0], indels[loci][1]))
  locus = dbc.fetchall()

  lengths = []
  counts  = []
  for look in locus:
    lengths.append(look[0])
  counts.append(look[1])

  mat_lengths = np.zeros([len(lengths),len(lengths)], dtype=int)
  mat_counts  = np.zeros([len(lengths),len(lengths)], dtype=int)

  col = 0
  for y in range(0, len(lengths)):
    for x in range(col, len(lengths)):
      if (x == y): 
        mat_lengths[x,y] = lengths[x]
        mat_counts[x,y] = counts[x]
      else:
        mat_lengths[x,y] = lengths[x] - lengths [y]
        mat_counts[x,y] = counts[x] + counts[y]
    col = col + 1

for a,b in np.nditer([mat_lengths, mat_counts]):
  print "%d : %d" % (a,b)

print "\n"
print mat_lengths
print "\n"
print mat_counts
print "\n"
print np.amax(mat_lengths)
"""
"""
with open("/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/test_out.txt", "w") as writer:
    writer.write(u'#'+'\t'.join(str(out) for out in mat_lengths) + '\n')    


for test in range (0,10):
    print indels[test]
    print indels[test][0]

    for x,y in np.nditer([test3, test_lengths]):
    print "%d:%d" % (x, y)

    for look in locus:
    test3.append(look[0])
    lengths = np.matrix(locus[0])
    counts  = np.matrix(locus[1])
    alles    = np.matrix(locus)
    print test3
"""

total = float(106)    # 106 accessions other than Oxford [107 in total]
threshold = 20  # threshold length of allelic difference to be useful
#scores = np.zeros((len(indels),2))
scores = np.zeros((3,4))
sstring = []

for loci in range(0,3):
  dbc.execute("SELECT length, count(*) from variants where chromosome=%s and position=%s group by length;", (indels[loci][0], indels[loci][1]))
  locus = np.array(dbc.fetchall())
  string = str(int(total) - sum(locus[:,1]))
  locus = locus[abs(locus[:,0]) > threshold]
  scores[loci,1]=len(locus)+1
  scores[loci,2]=1
  scores[loci,3]=sum(locus[:,1])
  if (len(locus[locus[:,0]<0]) < len(locus)):
    locus[:,0] = locus[:,0] + abs(min(locus[:,0]))
  print string 
  test = 1
  at = 0
  current = np.min(locus[loci,0])
  for x in range(0, len(locus)):
    if (locus[x,0] > (test * threshold)) :
      string = string + "+" + str(sum(locus[at:x-1,1]))
      scores[loci,2] = scores[loci,2] + 1
      at=x
      test = test + 1
  if (at < x):
    string = string + "+" + str(sum(locus[at:x,1])) 
    scores[loci,2] = scores[loci,2] + 1
  sstring.append( string )
  scores[loci,0] = scores[loci,3] / (0.5 * total ) 

  print str(scores[loci]) + " > " + string

"""
if (abs(locus[x][0]) >= (test * threshold)):
      scores[loci, 0] = scores[loci, 0] + ( counts[x] / (0.5 * total) )
      temp = x
    if (abs(locus[x][0]) >= lengths[temp] + (test * threshold)):
      temp1 = locus[x][1]
      test=test+1
      temp=x
      for y in range(x+1, len(lengths)-1):
        if (locus[y][1]>temp1 and abs(locus[y][0]) < (test * threshold)):
          temp=x
      scores[loci, 0] = scores[loci, 0] + ( counts[temp] / (0.5 * total) )
    
    
    if (abs(l engths[x]) >= (test + 1) * threshold):
      test = test + 1
    if (abs(lengths[x]) > test * threshold):
      scores[loci, 0] = scores[loci, 0] + ( counts[x] / total )
    for y in range(x+1,len(lengths)):
      print "y: " + str(scores[loci,0]) 
      if (abs(lengths[x] - lengths[y]) > threshold):
        if (abs(lengths[x]) >= threshold):
          scores[loci, 0] = scores[loci, 0] + ( (counts[x] + counts[y]) / total) 
        else:
          scores[loci, 0] = scores[loci, 0] + ( (counts[y]) / total) 
"""

"""
      elif (scores[locus][1] == 0):
        continue
      else:
        scores[locus][1] = scores[locus]
"""











# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
