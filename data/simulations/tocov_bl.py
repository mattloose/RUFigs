#!/usr/bin/env python

import csv
import sys


###Generate list of zeros
listofzeros = [0] * 49000
oddlistofzeros = [0] * 49000
evenlistofzeros = [0] * 49000

#print listofzeros

print sys.argv[1]
with open(sys.argv[1], 'rb') as f:
    reader = csv.reader(f)
    inputdata = list(reader)

f = open('ALL2_bl.txt', 'w')
fe = open('EVEN2_bl.txt', 'w')
fo = open('ODD2_bl.txt','w')

for row in inputdata[1:]:
    if row[8] != "Stopped":
        if row[1]=="J02459":
            print row, "start", int(row[2]) , "end",(int(row[2])+int(row[7]))
            for i in range(int(row[2]),(int(row[2])+int(row[7]))):
                listofzeros[i] +=1
                if int(row[0]) % 2 == 0:
                    evenlistofzeros[i] +=1
                else:
                    oddlistofzeros[i]+=1
        else:
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
    stringout = "simref"+"\t"+str(ix)+"\t"+str(c)+"\n"
    #if ix%50==0:
    f.write(stringout)

print "NEXT"

for ix, c in enumerate(oddlistofzeros):
    print "simref",ix,c
    stringout = "simref"+"\t"+str(ix)+"\t"+str(c)+"\n"
    #if ix%50==0:
    fo.write(stringout)

print "NEXT"
for ix, c in enumerate(evenlistofzeros):
    print "simref",ix,c
    stringout = "simref"+"\t"+str(ix)+"\t"+str(c)+"\n"
    #if ix%50==0:
    fe.write(stringout)

print "NEXT"
