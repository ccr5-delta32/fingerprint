#!/usr/bin/python

""" --- InDel filtering to prepare the list of indels for 'fingerprint.py' --- 
    Bjorn Pieper. MPIPZ. July 2014    """

import MySQLdb
import csv

# Filter InDels for not having ambiguous calls anywhere [done once]
with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/ma_indel_filtrd.csv', 'a') as filtrd:
  for loci in range(0, len(indels)):
    dbc.execute("SELECT chromosome, position, length, count(*), count(distinct(length)) from variants where chromosome=%s and position=%s and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%N%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%S%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%Y%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%M%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%K%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%R%%' and (select group_concat(concat(Oxford, consensus)) from variants where chromosome=%s and position=%s) not like '%%W%%'  group by position;", (indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1],indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1], indels[loci][0], indels[loci][1]))
    locus = dbc.fetchall()
    print "locus " + str(loci) + str(locus) 
    if locus != ():
      stuff = ''
      for x in range(0, len(locus[0])-1):
        stuff = stuff + str(locus[0][x]) + ','  
      stuff = stuff + str(locus[0][len(locus[0])-1])
      filtrd.write(stuff + "\n")

