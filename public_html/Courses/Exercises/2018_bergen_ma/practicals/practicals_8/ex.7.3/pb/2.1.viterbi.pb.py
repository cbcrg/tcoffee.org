#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import sys
import re
import math

LOG_ZERO=-999999999
SMALL=0.00000001

def log_multiply (v1, v2):
   if v1==LOG_ZERO or v2 ==LOG_ZERO:
      return LOG_ZERO
   return v1+v2

def viterbi (Data,Model):
  maxk=LOG_ZERO
  ptr_k=""
  L=len(Data)
  States=Model.keys()
  

  V   ={}#This will contain your Viterbi MAtrix - the best scores
  PTR ={}#This will contain the traceback
  for i in range (0, L+1):
      V[i]={}
      PTR[i]={}
        
  for k in (States):
      V[0][k]=0
            
  for i in range (1,L+1):
      symbol=Data[i-1]
     
      for l in (States):
        max_k=LOG_ZERO
        ptr_k=""
        for k in (States):
            v=#missing_1
            
            if v>max_k or max_k==LOG_ZERO:
                #use log multiply rather than multiplying probalilities
                #missing_2
                #missing_3
                
        V[i][l]=#missing_4
        PTR[i][l]=#missing_5
        
  max_k=LOG_ZERO
  ptr_k=""

  for k in (States):
    vv=V[L][k]
    if vv>max_k or max_k==LOG_ZERO:
        max_k=vv
        ptr_k=k
  vd={}#this will contain the list of decoded states
  for i in range (L, 0,-1):
      vd[i-1]=ptr_k
      ptr_k=PTR#missing_6

  return vd

def readmodel (name):
   Model={}
   with open(name) as f:
    for line in f:
        mobj=re.match ( r'(\S*)::(\S*)::(\S*)', line)
        if mobj:
            s1=mobj.group(1)
            s2=mobj.group(2)
            p=float (mobj.group(3))
            if (p<SMALL): p=LOG_ZERO
            else: p=math.log(p)

            if s1 not in Model:
                Model[s1]={}

            Model[s1][s2]=p

    return Model

def readdata (name):
   Data={}
   Ref={}
   i=0
   with open(name) as f:
    for line in f:
        Data[i]=line[0]
        Ref [i]=line[1]
   
        i+=1
   
   return (Data,Ref)

if __name__ == '__main__':
    (Data,Ref)=readdata   (sys.argv[1])
    Model=readmodel (sys.argv[2])
    ViterbiD=viterbi (Data, Model)
    
    L=len(Data)
    for i in range (0,L):
       #we will print the original state as provided in data along with the prediction
       sys.stdout.write("%c%c%c\n"%(Data[i],Ref[i],ViterbiD[i]))