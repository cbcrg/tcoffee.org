#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import sys
import re
import random

pstate=0


emit={}
emit['L']=[0]*10
emit['F']=[0]*10

temit={}
temit['F']=0
temit['L']=0

transitions={}
ttrans={}
for s1 in ("L","F"):
    transitions[s1]={}
    ttrans[s1]=0
    for s2 in ("L", "F"):
        transitions[s1][s2]=0


with open(sys.argv[1]) as f:
    for line in f:
        state=line[1]
        emit[state][int(line[0])]+=1
        temit[state]+=1
        if pstate:
            transitions[pstate][state]+=1
            ttrans[state]+=1
        pstate=state

print ("Transitions")        
for s1 in ("L","F"):
        for s2 in ("L", "F"):
            sys.stdout.write ("%s::%s::%.3f\n"%(s1,s2, transitions[s1][s2]/float(ttrans[s1])))

print ("Emission")
for s1 in ("L","F"):
    sys.stdout.write ("%c Dice:\n"%(s1))
    for s2 in range (1,7):
        sys.stdout.write ("%s::%d::%.3f\n"%(s1,s2,emit[s1][s2]/float(temit[s1])))
        


                              
        
        
