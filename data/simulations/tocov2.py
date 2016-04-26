#!/usr/bin/env python

import csv
import sys


###Generate list of zeros
listofzeros = [0] * 48600
oddlistofzeros = [0] * 48600
evenlistofzeros = [0] * 48600

#print listofzeros

print sys.argv[1]
with open(sys.argv[1], 'rb') as f:
    reader = csv.reader(f)
    inputdata = list(reader)

for row in inputdata[1:]:
    if row[8] != "Stopped":
        if row[1]=="J02459_F":
            if ((7000 <=int(row[2])<=12000) or (27000 <= int(row[2])<=32000)):
                #print row, "start", int(row[2]) , "end",(int(row[2])+int(row[7]))
                for i in range(int(row[2]),(int(row[2])+int(row[7]))):
                    listofzeros[i] +=1
                    if int(row[0]) % 2 == 0:
                        evenlistofzeros[i] +=1
                    else:
                        oddlistofzeros[i]+=1
        else:
            if  ((10503 <=int(row[2])<=15503) or (30503 <= int(row[2])<=35503)):
                seqlen=48502
                #print row, "start", seqlen-(int(row[2])+int(row[7])) , "end", seqlen-int(row[2])
                for i in range(seqlen-(int(row[2])+int(row[7])),seqlen-int(row[2])):
                    listofzeros[i] +=1
                    if int(row[0]) % 2 == 0:
                        evenlistofzeros[i] +=1
                    else:
                        oddlistofzeros[i]+=1
        #print "done"

print "NEXT"

for ix, c in enumerate(listofzeros):
    print "simref",ix,c

print "NEXT"

for ix, c in enumerate(oddlistofzeros):
    print "simref",ix,c

print "NEXT"
for ix, c in enumerate(evenlistofzeros):
    print "simref",ix,c

print "NEXT"
