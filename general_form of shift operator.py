import numpy as np
import math




#general Shift operator- q plate

k=3     
delta= math.pi
matrixS = np.zeros((2*k,2*k),dtype=complex)
i=1
while (i+2) < 2*k :
    matrixS[i][i+2] = 1.0
    i += 2
   
j=0
while (j+2) < 2*k :
    matrixS[j+2][j] = 1.0
    j += 2
    
k=3     
#the one I use
qplate = np.zeros((2*k,2*k),dtype=complex)
i=1
while (i+2) < 2*k :
    qplate[i+2][i] = 1.0
    i += 2
   
j=0
while (j+2) < 2*k :
    qplate[j][j+2] = 1.0
    j += 2

k=3     
inv_qplate = np.zeros((2*k,2*k),dtype=complex)
i=1
while (i+2) < 2*k :
    inv_qplate[i][i+2] = 1.0
    i += 2
   
j=0
while (j+2) < 2*k :
    inv_qplate[j+2][j] = 1.0
    j += 2

initial = np.zeros((2*k,1),dtype=complex)

initial[2][0]= 1.
initial[3][0]=1.5
initial/= np.linalg.norm(initial)



'''
for i in range (0,2*k,1) :
    matrixS[0+i][0+i] = 0.0
    if (i+3)< 2*k :
        matrixS[i][i+3] = 
    if (i-3) > 0 :
        matrixS[i][i-3] = - 1j*math.sin(delta/2)'''
    
invS1 = np.zeros((2*k,2*k),dtype=complex)
matrixS1 = np.zeros((2*k,2*k),dtype=complex)
for i in range (0,2*k,2) :
    invS1[0+i][0+i] =1.
    matrixS1[0+i][0+i] =  1.
    if (i+3)< 2*k :
        invS1[1+i][3+i] = 1. #S-1
        matrixS1[3+i][1+i] = 1.