# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 17:00:23 2019

@author: gonca
"""

import sympy as sym
import math
import matplotlib.pyplot as plt
import numpy as np
import random

pi = math.pi
points=100
maxint=20

x, t = sym.symbols("x t") 


'''Zeros de funções Reais'''

class graph:   
    def plot(f,x1,x2,points):
        xaxis = np.linspace(x1,x2,points) # 100 linearly spaced numbers
        yaxis=[]
        for elem in xaxis:
            y=f.subs(x,elem)
            yaxis.append(y)
        plt.axis(option='on')
        plt.grid(True, which='both')
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.plot(xaxis,yaxis)

class zero:
    def init(self):
        self.maxiterations = 20
    def bisec(fun,a,b,e):
        #if the algorithm doesn=t converge you may consider
        #there is a method problem
        k=0
        x1=a
        x2=b
        fa=fun.subs(x,a)
        fb=fun.subs(x,b)
        if fa*fb>0:
            print('Something is Wrong')
            return
        else:
            while True:
                xk=(a+b)/2
                if fb*fa<0:
                    b=xk
                if fb*fa>0:
                    a=xk
                if fb==0:
                    root=b
                    break
                if fa==0:
                    root=a
                    break
                fa=fun.subs(x,a)
                fb=fun.subs(x,b)
                if abs(b-a)<=e:
                    break
                k=k+1
                root = b
                if k>maxint:
                    break
              
        graph.plot(fun,x1,x2,points)
        print (k)
        
        return float(root)
    
    def falsepos(fun,a,b,e):
        k=0
        x1=a
        x2=b
        fa=fun.subs(x,a)
        fb=fun.subs(x,b)
        while True:
            xk=(a*abs(fb+b)*fa)/(fb+fa)
            if abs(fb)>abs(fa):
                b=xk
            if abs(fa)>abs(fb):
                a=xk
            if fb==0:
                root=b
                break
            if fa==0:
                root=a
                break
            fxk=fun.subs(x,xk)
            if abs(fxk)<e:
                root=xk
                break
            fa=fun.subs(x,a)
            fb=fun.subs(x,b)
            if abs(b-a)<=e:
                break
            k=k+1
            root = b
            if k>maxint:
                break
          
        graph.plot(fun,x1,x2,points)
        print(k)
        
        return float(root)
    
    def mpf(fun,ifun,x0,e):
        k=0
        xk=x0
        xkprev=40
        while True:
            g=ifun.subs(x,xk)
            xkprev=xk
            xk=g
            if abs(xk-xkprev)<e:
                break
            k=k+1
            if k>maxint:
                break
        root = xk
        
        y=x
        graph.plot(fun,(-2*x0),(2*x0),points)
        graph.plot(g,(-2*x0),(2*x0),points)
        graph.plot(y,(-2*x0),(2*x0),points)
        print(k)
        
        return float(root)
    
    def nr(fun,x0,e):
        k=0
        dfun=sym.diff(fun,x)
        xk=x0
        while True:
            fxk=fun.subs(x,xk)
            dfxk=dfun.subs(x,xk)
            g=xk-fxk/dfxk
            xkprev=xk
            xk=g
            if abs(xk-xkprev)<e:
                break
            k=k+1
            if k>maxint:
                break
        root=xk
        
        if x0==0:
            x1=-2
            x2=2
        else:
            x1=x0
            x2=x0
            
        graph.plot(fun,(2*x1),(2*x2),points)
        print(k)
        
        return float(root)
    
    def sec(fun,a,b,e1,e2=None):
        if e2==None:
            e2=e1
        xk=b
        xkprev=a
        k=0
        while True:
            fxk=fun.subs(x,xk)
            #if fxk<e1:
                #break
            fxkprev=fun.subs(x,xkprev)
            g=((xkprev*fxk)-(xk*fxkprev))/(fxk-fxkprev)
            xkprev=xk
            xk=g
            if abs(xk-xkprev)<e1:
                break
            k=k+1
            if k>maxint:
                break
        
        root=xk
    
        graph.plot(fun,a,b,points)
        print(k)
        
        return float(root)
                
            
    def compare(fun,a,b,e):
        x0 = random.randint(a,b)
        roots = []
        
        try:
            y=bisec(fun,a,b,e)
            roots.append(y)        
        except ValueError:
            print ('Not good results on bisec method')
    
        try:
            y=falsepos(fun,a,b,e)
            roots.append(y)        
        except ValueError:
            print ('Not good results on false position method')
    
        try:
            y=mpf(fun,x0,e)
            roots.append(y)        
        except ValueError:
            print ('Not good results on mpf method')
    
        try:
            y=nr(fun,x0,e)
            roots.append(y)        
        except ValueError:
            print ('Not good results on Newtons method')
    
        try:
            y=sec(fun,a,b,e)
            roots.append(y)        
        except ValueError:
            print ('Not good results on secante method')
            
        graph.plot(fun,3*a,3*b,points)
            
        return roots



'''Linear Systems'''
class linear:
    def init(self):
        self.maxiterations = 10
    
    def fullPivot(matrix,b):
         p = np.max(abs(matrix))
         indices = np.where(matrix==p)
         i,j=indices
        
    
    def gaussSeidel(A,x0,b,e):
        y=x0.copy() #get a copy of the guess
        k=0 #number of iterations
        i=0 #number of the line
        j=0 #number of the collumn
        t=0 #number of the
        erro=0 #initiates the dif value
        n=np.ma.size(x0) #gets the number of lines of the matrix
        summ=0
        erArray=np.zeros(n)
        
        while True: #infinite loop while not satisfied the stop criterion
            while i<n: #while i<number of lines
                while j<n: #while j<number of columns
                    if j!=i: #not getting [i][j]
                        summ=summ+y[j]*A[i][j] #equation
                    j=j+1 #increases column
                value = (b[i]-summ)/A[i][i]
                y[i] = value #gets the correspondig y or x
                j=0 #initiates again number of columns
                summ=0 #initiates the summ again
                i=i+1 #go to the next line
            k=k+1 #increases the interactions
            
            while t<n: #stop criterion
                er=abs(y[t]-x0[t]) #|xk+1 - xk|
                erArray[t] = er
                t=t+1
                
            erro=np.max(erArray)
            if erro<e:
                break #gets out of the infinite loop
            else:
                t=0
                i=0
                x0=y.copy()#the next x0 will be the last y
            
        print ('The method took ' +str(k)+ ' iterations to get a result')
        print ('x vector is')
        print (y)
        
        return (y,k) 

'''             
            while t<n: #stop criterion
                er=abs(y[t]-x0[t]) #|xk+1 - xk|
                if er>erro:#gets the maximum of the errors
                    erro=er
                t=t+1
            if erro<e:
                break #gets out of the infinite loop
            else:
                t=0
                i=0
                x0=y.copy()#the next x0 will be the last y
'''

        
        
        
    