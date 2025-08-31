# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:18:30 2019

@author: ahnux
"""

import numpy as np
#from scipy.integrate import odeint
from scipy import linspace
from scipy.integrate import solve_ivp
"""
定义常微分方程，给出各方向导数,即速度
"""
def reaction(Protein):

    A, B, AmBn, ac_B, C, ac_BnCh, ac_C, ac_acB = Protein

    R0 = k1 * A**m * B**n
    R1 = k2 * AmBn
    R2 = k3 * ac_B**n * C**h
    R3 = k4 * ac_BnCh
  
    return np.array([R0, R1, R2, R3])

def chainode(t, Protein):

    A, B, AmBn, ac_B, C, ac_BnCh, ac_C, ac_acB = Protein
    
    V = reaction(Protein)
    
    dA = -m*V[0]
    dB = -n*V[0]
    dAmBn = V[0]-V[1]
    dac_B = n*V[1]-n*V[2]
    dC = -h*V[2]
    dac_BnCh = V[2]-V[3]
    dac_C = h*V[3]
    dac_acB = n*V[3]
    
    
    return np.array([dA, dB, dAmBn, dac_B, dC, dac_BnCh, dac_C, dac_acB])

if __name__ =='__main__':
    
    dt = linspace(0,36000,3600000)
    
    k1=1e-5
    k2=1.0
    k3=1e-5
    k4=1.0
    
    output = open('1pattern of ratio.dat','w+')
    
    for i in np.arange(0,41,1):
        
        m=1+0.1*i
        
        for j in np.arange(0,41,1):
            
            n=1+0.1*j
            
            for k in np.arange(0,41,1):
                
                h=1+0.1*k
                
                print(h)
                protein_init = [10.0, 10.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0]
                
                P = solve_ivp(chainode, (0, 3600), protein_init, method='BDF')#, t_eval=dt
                time = P.t
                protein = P.y
            
                sum_b = str(protein[3,-1] + protein[7,-1])
                sum_c = str(protein[6,-1])
                sum_bc = str(protein[3,-1] + protein[7,-1] + protein[6,-1])
                
                output.write(str(m))
                output.write('\t')
                output.write(str(n))
                output.write('\t')
                output.write(str(h))
                output.write('\t')
                output.write(sum_b)
                output.write('\t')
                output.write(sum_c)
                output.write('\t')
                output.write(sum_bc)
                output.write('\n')
    output.close()
