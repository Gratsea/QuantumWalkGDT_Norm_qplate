# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 12:15:59 2018

@author: kgratsea
"""

# -*- coding: utf-8 -*-
"""
Quantum walk with GDT. 3 free parameters for each coin operator

function to be minimzed norm

@author: kgratsea
"""

import numpy as np
import math
from scipy import optimize
global final
import random

metrhths=0
metrhths_prNORM=0
previous_NORM = 1

def tensor(vectorA,vectorB) :
    m = np.size(vectorA,0)
    n = np.size(vectorB,0)
    tens=np.zeros((m,n))
    for i in range(m) :
        for j in range(n) :
            tens[i][j] = vectorA[i]*vectorB[j]
    return (tens);

def func(z) :  
    global metrhths 
    global metrhths_prNORM
    global previous_NORM 
    metrhths += 1

    n=51 #number of steps
    k=2*n+1 #number of sites at the final state
    
    initial = np.zeros((2*k,1),dtype=complex)
    #localised on one site
    initial[2*n][0]= 1.
    initial[2*n+1][0]= 1.5
    initial/= np.linalg.norm(initial)

   
    #definition of matrixS  
    #(m,up)--> (m+1,up) (m,down)--> (m-1,down)

    qplate = np.zeros((2*k,2*k),dtype=complex)    
    i=1
    while (i+2) < 2*k :
        qplate[i][i+1] = 1.0
        i += 2
       
    j=1
    while (j+2) < 2*k :
        qplate[j+1][j] = 1.0
        j += 2
        
    matrixS = qplate    
        
    listSt = []
    listc = []
    listC = []

    listSt.append (initial)
    
    #Define coin operators paper Chandrashekar PRA,032326 (2018) 
    
    l = 0 # for corresponding the correct coin parameters at each step n
    for j in range (0,n,+1) : 
        #print ("n",j)
        c=np.zeros((2,2),dtype=complex)
        theta=z[0+l]
        ksi=z[1+l]
        zeta=z[2+l]
        
        c[0][0]= (math.cos(ksi*math.pi) + math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)
        c[0][1]= (math.cos(zeta*math.pi) + math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2) 
        c[1][0]= (math.cos(zeta*math.pi) - math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2)         
        c[1][1]= - (math.cos(ksi*math.pi) - math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)  
        
        listc.append(c)
        matrixC = np.zeros((2*k,2*k),dtype=complex)
        #print (c)
        
        for i in range (0,2*k,2):
            matrixC[0+i][0+i] = c[0][0]
            matrixC[1+i][1+i] = c[1][1]
            matrixC[0+i][1+i] = c[0][1]          
            matrixC[1+i][0+i] = c[1][0]   
         
        listC.append (matrixC)    
        
        #print (initial)
        m1 = np.dot(matrixC,initial)
        m2 = np.dot(matrixS,m1)   #next state
        #print (m1)
        #print (m2)
        listSt.append (m2)
        initial = m2/np.linalg.norm(m2)
        l += 3 # moving to the next coin parameters
        
    Phi=initial    
    #create state of oute sites and find max. ent. there
    Phi_outer=np.zeros((2,2),dtype=complex)
    Phi_outer[0][0] = Phi[0][0]
    Phi_outer[1][0] = Phi[1][0]
    Phi_outer[0][1] = Phi[2*k-2][0]
    Phi_outer[1][1] = Phi[2*k-1][0]

    
    psiA, l, psiB = np.linalg.svd(Phi_outer,full_matrices=1) #decomposition of initial matrix
    #print ("l",l)
    
    NORM=0.0
    p=1.0 # p has to be larger or equal than 1 for the algorithm to work
    
    m = np.size(Phi_outer,0)  #number of rows of initial matrix
    n = np.size(Phi_outer,1)  #number of columns of initial matrix
    sum=np.zeros((m,n))
    
    for i in range(2) :
         tens= tensor(psiA[i],psiB[i])
         sum += math.pow(l[i],p-1)*tens
         NORM = NORM + math.pow(l[i],p) 

    NORM = math.pow(NORM,1./p)
    #print (NORM)
    #print (previous_NORM,NORM)
    #print (previous_NORM,NORM)

    with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_only_NORM.txt', 'a+') as f:
        print (metrhths,",",NORM,file=f)
    f.close()


    if (abs(NORM - previous_NORM) > 0.000001) and (NORM - 1) > 0.01 :
        with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound_only_NORM.txt', 'a+') as f:
            print (metrhths_prNORM,",",NORM,file=f)
            metrhths_prNORM += 1

        f.close()
        if (-NORM+math.sqrt(2)<0.0000001) :
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("initial")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (initial,file=f)
            f.close()
            
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("l,NORM")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (l,NORM,file=f)
            f.close()
        
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("z")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (z,file=f)
            f.close()
            
            
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("listc")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (listc,file=f)
            f.close()
            
            
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("Phi")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (Phi,file=f)
            f.close()
            
            f = open("QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt","a+")
            f.write("Phi_outer")
            f.close()
            with open('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound.txt', 'a+') as f:
                print (Phi_outer,file=f)
            f.close()    
    
    previous_NORM  = NORM

    
    return (-NORM+ math.sqrt(2))

    
f=51
    
my_randoms=[]
for i in range (3*f):
    my_randoms.append(random.randrange(1,10,1))

print (my_randoms)    
    
initial_coin_parameters=my_randoms

        
minimizer_kwargs = {"method": "BFGS"}


ret = optimize.basinhopping(func,initial_coin_parameters, minimizer_kwargs=minimizer_kwargs,niter=40, T=1.0, disp = True)
  
l=0
listc=[]
for j in range (0,f,+1) : 
        #print ("j",j)
        c=np.zeros((2,2),dtype=complex)
        theta=ret.x[0+l]
        ksi=ret.x[1+l]
        zeta=ret.x[2+l]
        
        c[0][0]= (math.cos(ksi*math.pi) + math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)
        c[0][1]= (math.cos(zeta*math.pi) + math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2) 
        c[1][0]= (math.cos(zeta*math.pi) - math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2)         
        c[1][1]= - (math.cos(ksi*math.pi) - math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)  
        
        listc.append(c)
        l+=3
 
    
'''      
import matplotlib.pyplot as plt
import numpy as np


x, y = np.loadtxt('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound_only_NORM.txt', delimiter=',', unpack=True)
plt.plot(x,y,'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Schmidt norm after n=5 steps')
plt.legend()
plt.savefig('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound_only_NORM.png')

a, b = np.loadtxt('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound_only_NORM.txt', delimiter=',', unpack=True)
plt.plot(a, b,'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Schmidt norm after n=5 steps ( points ) ')
plt.legend()
plt.savefig('QW_NORM_Phiouter_qplate-Coin_n20_T1_steps51_bound_only_NORM.png')
'''
