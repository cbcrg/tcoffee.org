#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import sys
import re
import random
outcomes = '123456'

    
class Model:
    def __init__(self, pFL, pLF,p6):
        self.emit = list(outcomes)
        self.fairProb = [1.0/6] * 6
        
        self.loadProb = # missing1
        self.tranProb = { 'F':pFL,'L':pLF }
                                
    def setState(self,which=None):
        if which:  self.state = which
        else:      self.state = random.choice('FL')



    def rollDice(self):
        if self.state == 'F':  L = self.fairProb
        else:                  L = self.loadProb
        f = random.random()
        S = 0
        for i in range(len(L)):
            #missing2
            if f < S:  break
        return self.emit[i]

    def transit(self):
        def switch():
            if self.state == 'F':  self.state = 'L'
            else:                  self.state = 'F'

        p = # missing3 
        f = # missing4
        if f < p:  switch()

    def sequence(self,N=100):
        rolls = list()
        states = list()
        for i in range(N):
            sys.stdout.write("%s%s\n"%#missing5
            # missing6
        return 
                



if __name__ == '__main__':
    
    m = Model(pFL=float(sys.argv[1]),pLF=float(sys.argv[2]), p6=float(sys.argv[3]))
    m.setState('F')
    m.sequence(int(sys.argv[4]))
       
