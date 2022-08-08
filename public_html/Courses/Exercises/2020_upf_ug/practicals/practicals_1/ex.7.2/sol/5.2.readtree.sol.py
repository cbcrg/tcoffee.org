#!/usr/bin/env python

import sys
import re
class Node:
    def __init__(self):
        self.name=""
        self.distance=float (-1)
        self.bootstrap=float (-1)
        self.left=0
        self.right=0
        self.parent=0

def declare_new_tree_node(nodes,nn):
    nodes[nn]=Node()
    return (nn, nn+1)



def scan_name_and_dist (i, l):
   name=""
   number=""
   if l[i]==';':
       return ("",-1,i)
     
   while l[i]!=':' and i<len(l) and l[i]!=')' and l[i]!=';' and l[i]!=',':
       name+=l[i]
       i+=1
     
   if l[i]!=':':
       distance=float(0)
       return (name,distance, i)
   else:
       i+=1
     
   while  str.isdigit(l[i]) or l[i]=='e' or l[i]=='-' or l[i]=='.':
       number+=l[i]
       i+=1

   number=float(number)
   return (name, number,i)

def newick2nodes (line):
   
   nodes={}
   nodes[0]=-1
   nn=1 # root starts at 1
   (N,nn)=declare_new_tree_node(nodes,nn)
   T=R=N
   
   c=pi=i=0
   while (line[i])!=';':
      c=line[i]
      i+=1
      if c=='(':
         (N,nn)=declare_new_tree_node(nodes,nn)
         nodes[N].parent=T

         if nodes[T].right==0:
            nodes[T].right=N
         elif nodes[T].left==0:
            nodes[T].left=N
         else:
            nodes[N].right=nodes[T].right
            nodes[nodes[T].right].parent=N
            
            nodes[N].left=nodes[T].left
            nodes[nodes[T].left].parent=N
            
            nodes[T].right=N

            (N,nn)=declare_new_tree_node(nodes,nn)
            
            nodes[T].left=N
            nodes[N].parent=T

         T=N
         lastc=0
        
      elif c==')':
        T=nodes[T].parent
        (nodes[T].name,nodes[T].distance,i)=scan_name_and_dist (i,line)
        if nodes[T].name and nodes[T].name[0]:
            nodes[T].bootstrap=float(nodes[T].name)
            nodes[T].name=""
        lastc=0;
        
      elif c==',':
        T=nodes[T].parent;
        lastc+=1
      else:
        (N,nn)=declare_new_tree_node(nodes,nn)
        nodes[N].parent=T

        if nodes[T].right==0:
            nodes[T].right=N
        elif nodes[T].left==0:
            nodes[T].left=N    
        else:
            nodes[N].right=nodes[T].right
            nodes[nodes[T].right].parent=N

            nodes[N].left=nodes[T].left
            nodes[nodes[T].left].parent=N

            nodes[T].right=N
            
            
            (N,nn)=declare_new_tree_node(nodes,nn)
            nodes[T].left=N
            nodes[N].parent=T
	    	  
        T=N
        i=i-1
        
        (nodes[T].name,nodes[T].distance,i)=scan_name_and_dist (i,line);
        lastc=0
        
   T=nodes[T].parent
   
   if nodes[T].right==0 and nodes[T].left!=0:
      T=nodes[T].left
   elif nodes[T].right!=0 and nodes[T].left==0:
      T=nodes[T].right

   nodes[T].parent=-1
   return (nodes,nn)

#The main code starts here


tree=""

with open (sys.argv[1]) as f:
    for line in f:
        tree+=line

tree=re.sub (r'[ \n\t\r]',"",tree)

#get the tree
nn=0
nodes={}
(nodes, nn)=newick2nodes(tree)


#Print list of sister species
#check that a given node has a left and a right child
#check that these children are terminals
#if they are print them

for n in nodes:
    if n>0:
        l=nodes[n].left
        r=nodes[n].right
        if l and r and nodes[l].name!="" and nodes[r].name!="":
            print ("sisters",nodes[l].name,nodes[r].name)


