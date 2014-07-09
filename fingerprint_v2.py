#!/usr/bin/python

""" --- Select good multi-allelic InDels for finger printing --- 
    Bjorn Pieper. MPIPZ. July 2014    """

import MySQLdb
import csv
import numpy as np
import re
from copy import deepcopy

with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/ma_indel_filtrd.csv') as ID:
  indels = list(csv.reader(ID, delimiter=","))

## debug subset
#indels = [indel for indel in indels if indel[1] == '261863' or indel[1] == '137284']
#indels = indels[0:305]
#db  = MySQLdb.connect(host="localhost", user="robot", passwd="giveusthetechnology", db="chpoly")
#dbc = db.cursor()

n = len(indels) 
m = max([int(x[4]) for x in indels])
total = float(107)    # 107 accessions with Oxford included 
threshold = 30  # threshold length of allelic difference to be useful
scores = np.zeros((n,4))
result1 = np.zeros((n,m))
result2 = np.zeros((n,m))
result3 = np.zeros((n,m))

kromo = []
pos = []
orinal = []
groups = []
egroup = []
score = []

clamp = lambda n, minn, maxn: max(min(maxn, n), minn)

def final(x):
  if (x == len(locus)-1):
    ox.append(cumm)
    result1[loci,x] = group
    result2[loci,x] = loc[x,1]
    result3[loci,x] = loc[x,0]

for loci in range(0,n):
  dbc.execute("SELECT length, count(*) from variants where chromosome=%s and position=%s group by length;", (indels[loci][0], indels[loci][1]))
  locus = np.array(dbc.fetchall())
  kromo.append(indels[loci][0])
  pos.append(indels[loci][1])
  orinal.append(len(locus))
  ox = [total - (sum(locus[:,1]) - sum(locus[abs(locus[:,0]) < threshold,1]))]
  if (max(abs(locus[:,0])) < threshold):
    egroup.append(total)
    score.append(0)
    groups.append(0)
    continue
  locus = locus[abs(locus[:,0]) >= threshold]
  loc = deepcopy(locus)
  scores[loci,1]=len(locus)+1
  scores[loci,2]=1
  scores[loci,3]=sum(locus[:,1])
  if (len(locus[locus[:,0]<0]) != 0 ):
    locus[:,0] = locus[:,0] + abs(min(locus[:,0]))

  group = 1
  cumm = 0
  if (min(locus[:,0]) >= threshold):
    thr = min(locus[:,0]) + threshold
  else:
    thr = threshold
  
  for x in range(0, len(locus)):
    if (locus[x,0] < (thr)):
      result1[loci,x] = group
      result2[loci,x] = loc[x,1]
      result3[loci,x] = loc[x,0]
      cumm = cumm + loc[x,1]
      final(x)
    elif ((locus[x,0] >= (thr)) & (group == result1[loci,x-1])):
      ox.append(cumm)
      cumm = loc[x,1] 
      group = group + 1
      thr = locus[x,0] + threshold
      result1[loci,x] = group
      result2[loci,x] = loc[x,1]
      result3[loci,x] = loc[x,0]
      final(x)
    else:
      result1[loci,x] = group
      result2[loci,x] = loc[x,1]
      result3[loci,x] = loc[x,0]
      cumm = cumm + loc[x,1]
      final(x)
  if (group > 1 ):
    score.append( sum( [clamp( (float(x)/(total)), 0, ((total/(group+1))/total)) for x in ox[0:len(ox)]] ) )
    egroup.append( ''.join([str(int(x)) + ';' for x in ox[0:len(ox)-1]]+[str(ox[len(ox)-1])]) )
  else:
    score.append( (min((total - ox[0]), ox[0])) / (0.5*total) ) 
    egroup.append(str(int(ox[0])) + ';' + str(int(total - ox[0])))
  groups.append(group)

with open("/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/fp_v2_out.txt", "w") as writer:
  for z in range(0,len(score)):
    writer.write(kromo[z] + ",")
    writer.write(str(pos[z]) + ",")
    writer.write(str(orinal[z]) + ",")
    writer.write(str(groups[z]) + ",")
    writer.write(str(score[z]) + ",")
    writer.write("'" + str(egroup[z]) + "',")
    writer.write( ''.join( [str(int(x)) + ',' for x in result1[z][0:len(result1[z])]] ) )
    writer.write( ''.join( [str(int(x)) + ',' for x in result2[z][0:len(result2[z])]] ) )
    writer.write( ''.join( [str(int(x)) + ',' for x in result3[z][0:len(result3[z])-1]] + [str(int(result3[z][len(result3[z])-1]))]) )
    writer.write("\n")

# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
