# -*- coding: utf-8 -*-
"""
Created on Sun May 12 20:22:23 2019

@author: Gonçalo
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:19:34 2019
Gauss Seidel Method for 3*3 Matrix
@author: gonca
"""

import numpy as np


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
    
    return  

#the entry arrays must be float, at least the guess
A=np.array([[6,-3,2],[2,5,-1],[3,1,5]])
b=np.array([4,1,1])
x0=np.zeros(3)
e=10**-7
'''
Another Example
A=np.array([[6,3,1],[4,9,-3],[1,-1,3]])
b=np.array([10,16,14])
x0=np.array([-1.8,5.5,7.3])
e=0.1
'''
gaussSeidel(A,x0,b,e)