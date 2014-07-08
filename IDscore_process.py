#!/usr/bin/python

""" --- Process scored indels to select the good ones --- 
    Bjorn Pieper. MPIPZ. July 2014    """

import MySQLdb
import csv
import numpy as np
import re
from copy import deepcopy

with open('/home/bpuser/ServHome/MPIPZ/genetic_resources/Molec_Markrs+Primers/fingerprint/test_out.txt') as ID:
  indels = list(csv.reader(ID, delimiter=","))

srt6 = sorted([indel for indel in indels if int(indel[3]) == 6], key= lambda indel: indel[4],reverse=True)

for x in range(0, min(4, len(srt6))):
    print srt6[x]
    print "\n"


