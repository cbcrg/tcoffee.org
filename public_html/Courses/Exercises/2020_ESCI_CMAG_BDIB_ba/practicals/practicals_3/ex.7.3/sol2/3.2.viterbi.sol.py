#!/usr/bin/env python

import sys
import re
import math
import random

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
  Ns=len (States)
  V   ={}
  PTR ={}
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
            v=log_multiply(V[i-1][k],Model[k][l]);
            
            if v>max_k or max_k==LOG_ZERO:
                max_k=v
                ptr_k=k
        V[i][l]=log_multiply(Model[l][symbol],max_k)
        PTR[i][l]=ptr_k
        
  max_k=LOG_ZERO
  ptr_k=""

  for k in (States):
    vv=V[L][k]
    if vv>max_k or max_k==LOG_ZERO:
        max_k=vv
        ptr_k=k
  vd={}
  for i in range (L, 0,-1):
      vd[i-1]=ptr_k
      ptr_k=PTR[i][ptr_k]

  return (vd,max_k)

def viterbi_trainning (Data, Model, States, Emits):

  (v,score)=viterbi(Data,Model)
  L=len(Data)

  A={}
  E={}
  for k in (States):
     A[k]={}
     E[k]={}
     for l in (States):
        A[k][l]=0
     for l in (Emits):  
        E[k][l]=0
        
  for i in range (1,L):
    A[v[i-1]][v[i]]+=1
    E[v[i]][Data[i]]+=1

  for k in (States):
    for l in (States):A[k][l]+=1
    for l in (Emits) :E[k][l]+=1

  for k in (States):
     num=float(0)
     for l in (States): num+=A[k][l]
     for l in (States): Model[k][l]=math.log(A[k][l]/num)

  for k in (States):
     num=float (0)
     for l in (Emits): num+=E[k][l]
     for l in (Emits): Model[k][l]=math.log(E[k][l]/num)

  (v,score)=viterbi(Data,Model)
  return (Model,score)   

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

    States=Model.keys()
    
    
    Emits=[]
    EmitsH={}

    for k in (States):
       k1=Model[k].keys()
       for l in (k1):
          if l not in Model and l not in EmitsH:
             Emits.append(l)
             EmitsH[l]=1
    
    return (Model, States, Emits)

def model2random (Model, States, Emits):
   for k in (States):
     n=len(States)

     tr=0
     for l in (States):
        r=random.random()
        tr+=r
        Model[k][l]=r
     for l in (States):
        Model[k][l]/=tr
        
     n=len(Emits)
     tr=0
     for l in (Emits):
        r=random.random()
        tr+=r
        Model[k][l]=r  
     for l in (Emits):
        Model[k][l]/=tr   
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
    (Model, States, Emits)=readmodel (sys.argv[2])

    
    bModel={}
    for l in (States):
       bModel[l]={}

    bscore=LOG_ZERO
    for j in range (0,100):
       pscore=0
       Model=model2random(Model,States,Emits)
       for i in range (0, 100):
          (Model,score)=viterbi_trainning (Data, Model,States,Emits)
          #sys.stderr.write("iteration %d Score=%d\n"%(i,score))
          if (score==pscore):break
          pscore=score
       sys.stderr.write("iteration %d Score=%d\n"%(j,pscore))
       if (pscore>bscore):
          bscore=pscore
          for l in (States):
             for k in (States):
                bModel[l][k]=Model[l][k]
             for k in (Emits):
                bModel[l][k]=Model[l][k]
    
    sys.stderr.write("BScore=%d\n"%(bscore))
    for l in (States):
       for k in (States):
          sys.stdout.write("%c::%c::%.3f\n"%(l,k,math.exp(bModel[l][k])))
       for k in (Emits):
          sys.stdout.write("%c::%c::%.3f\n"%(l,k,math.exp(bModel[l][k])))
                   
