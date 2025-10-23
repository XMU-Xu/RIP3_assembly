# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:18:30 2019

@author: ahnux
"""
import numpy as np
from scipy import linspace
from scipy.integrate import solve_ivp

"""
定义常微分方程，给出各方向导数,即速度
"""
#k1=1.0e-5*500 #TNF+TNFR1 = TNF_TNFR1;
k_1=1.0e-3
#k2=1.0e-3
#k3=1.5e-6*500 #TNFR1_active + TRADD = TNFR1_active_TRADD;
k_3=1.0e-3
#k4=1.0e-3
#k5=1.0e-8*500 #TRADD_active + FADD = TRADD_active_FADD;
k_5=1.0e-3
#k6=1.0e-3 #TRADD_active_FADD = TRADD_active_1 + FADD_bind_T;
#k7=1.0e-6*500 #TNFR1_active + RIP1 = TNFR1_active_RIP1;
k_7=1.0e-3
#k8=1.0e-3
#k9=1.0e-6*500 #RIP1_active + FADD = RIP1_active_FADD;
k_9=1.0e-3
#k10=1.0e-3
#k11=1.0e-15*500*500*500 #RIP1_active + 3*RIP3 = RIP1_active_RIP3_3;
k_11=1.0e-3
#k12=1.0e-3
#k13=1.0e-15*500*500*500 #RIP1_active_1 + 3*RIP3 = RIP1_active_1_RIP3_3;
k_13=1.0e-3
#k14=1.0e-3
#k15=1.0e-6*500 #RIP1_active_2 + FADD = RIP1_active_2_FADD;
k_15=1.0e-3
#k16=1.0e-3
#k17=1.0e-6*500 #FADD_bind_T + 2*ProCasp8 = FADD_bind_T_ProCasp8_2;
k_17=1.0e-3
#k18=1.0e-5
#k19=1.0e-6*500     #FADD_bind_R + ProCasp8 = FADD_bind_R_ProCasp8;
k_19=1.0e-3
#k20=1.0e-5
#k21=1.0e-12*500*500  #    3*RIP3_pho + 2*MLKL = RIP3_pho_3_MLKL_2;
k_21=1.0e-3
#k22=1.0e-3
k23=0 #;// 1.0e-2   FADD_bind_R_ProCasp8 + RIP3_pho = FADD_bind_R_ProCasp8_RIP3_pho;
k_23=1.0e-2
k24=0
k25=0
k26=0
k27=0
k28=0
k29=0
k30=0
k31=0
k32=0
k33=0
k34=0
k35=0
k36=0
k37=0
k38=0
k39=0
k40=0
k41=0
k42=0
k43=0
k44=0
k45=0
k46=0
k47=0
k48=0
k49=0
k50=0
k51=0
k52=0
k53=0
k54=0
k55=0
k56=0
k57=0
k58=0
k59=0
k60=0
k61=0
k62=0

def LHSample(D,bounds,N):
    result = np.empty([N, D])
    temp = np.empty([N])
    d = 1.0 / N

    for i in range(D):

        for j in range(N):
            temp[j] = np.random.uniform(low=j * d, high=(j + 1) * d, size = 1)[0]

        np.random.shuffle(temp)

        for j in range(N):
            result[j, i] = temp[j]

    #对样本数据进行拉伸
    b = np.array(bounds)
    lower_bounds = b[:,0]
    upper_bounds = b[:,1]
    if np.any(lower_bounds > upper_bounds):
        print('范围出错')
        return None

    #   sample * (upper_bound - lower_bound) + lower_bound
    np.add(np.multiply(result, (upper_bounds - lower_bounds), out=result),lower_bounds,out=result)
    return np.array(result)

def reaction(Protein):

    TNF, TNFR1, TRADD, FADD, ProCasp8, RIP1, RIP3, MLKL, TNF_TNFR1, TNF_1, TNFR1_active, TNFR1_active_TRADD, \
    TNFR1_active_1, TRADD_active_1, TRADD_active, TRADD_active_FADD,  FADD_bind_T, TNFR1_active_RIP1, \
    TNFR1_active_2, RIP1_active, RIP1_active_FADD, RIP1_active_1, FADD_bind_R, RIP1_active_RIP3_3, RIP1_active_2, \
    RIP3_pho, RIP1_active_1_RIP3_3, RIP1_active_3, RIP1_active_2_FADD,RIP1_active_4, FADD_bind_T_ProCasp8_2, \
    FADD_bind_T_1, Casp8, FADD_bind_R_ProCasp8, FADD_bind_R_1, RIP3_pho_3_MLKL_2, RIP3_pho_1, MLKL_pho, FADD_bind_R_ProCasp8_RIP3_pho = Protein

    R0 = k1 * TNF * TNFR1 - k_1 * TNF_TNFR1
    R1 = k2 * TNF_TNFR1
    R2 = k3 * TNFR1_active * TRADD - k_3 * TNFR1_active_TRADD
    R3 = k4 * TNFR1_active_TRADD
    R4 = k5 * TRADD_active * FADD - k_5 * TRADD_active_FADD
    R5 = k6 * TRADD_active_FADD
    R6 = k7 * TNFR1_active * RIP1 - k_7 * TNFR1_active_RIP1
    R7 = k8 * TNFR1_active_RIP1
    R8 = k9 * FADD * RIP1_active - k_9 * RIP1_active_FADD
    R9 = k10 * RIP1_active_FADD
    R10 = k11 * RIP1_active * RIP3**nRIP3 - k_11 * RIP1_active_RIP3_3
    R11 = k12 * RIP1_active_RIP3_3
    R12 = k13 * RIP1_active_1 * RIP3**nRIP3 - k_13 * RIP1_active_1_RIP3_3
    R13 = k14 * RIP1_active_1_RIP3_3
    R14 = k15 * FADD * RIP1_active_2 - k_15 * RIP1_active_2_FADD
    R15 = k16 * RIP1_active_2_FADD
    R16 = k17 * FADD_bind_T * ProCasp8**nPC8 - k_17 * FADD_bind_T_ProCasp8_2
    R17 = k18 * FADD_bind_T_ProCasp8_2
    R18 = k19 * FADD_bind_R * ProCasp8 - k_19 * FADD_bind_R_ProCasp8
    R19 = k20 * FADD_bind_R_ProCasp8
    R20 = k21 * RIP3_pho**nRIP3 * MLKL**nMLKL - k_21 * RIP3_pho_3_MLKL_2
    R21 = k22 * RIP3_pho_3_MLKL_2
    R22 = k23 * RIP3_pho * FADD_bind_R_ProCasp8 - k_23 * FADD_bind_R_ProCasp8_RIP3_pho
    R23 = k24 * TNF
    R24 = k25 * TNFR1
    R25 = k26 * TNF_TNFR1
    R26 = k27 * TNF_1
    R27 = k28 * TNFR1_active
    R28 = k29 * TRADD
    R29 = k30 * TNFR1_active_TRADD
    R30 = k31 * TNFR1_active_1
    R31 = k32 * TRADD_active
    R32 = k33 * FADD
    R33 = k34 * TRADD_active_FADD
    R34 = k35 * TRADD_active_1
    R35 = k36 * FADD_bind_T
    R36 = k37 * RIP1
    R37 = k38 * TNFR1_active_RIP1
    R38 = k39 * TNFR1_active_2
    R39 = k40 * RIP1_active
    R40 = k41 * RIP1_active_FADD
    R41 = k42 * RIP1_active_1
    R42 = k43 * FADD_bind_R
    R43 = k44 * RIP3
    R44 = k45 * RIP1_active_RIP3_3
    R45 = k46 * RIP1_active_2
    R46 = k47 * RIP3_pho
    R47 = k48 * RIP1_active_1_RIP3_3
    R48 = k49 * RIP1_active_3
    R49 = k50 * RIP1_active_2_FADD
    R50 = k51 * RIP1_active_4
    R51 = k52 * ProCasp8
    R52 = k53 * FADD_bind_T_ProCasp8_2
    R53 = k54 * FADD_bind_T_1
    R54 = k55 * Casp8
    R55 = k56 * FADD_bind_R_ProCasp8
    R56 = k57 * FADD_bind_R_1
    R57 = k58 * MLKL
    R58 = k59 * RIP3_pho_3_MLKL_2
    R59 = k60 * RIP3_pho_1
    R60 = k61 * MLKL_pho
    R61 = k62 * FADD_bind_R_ProCasp8_RIP3_pho

    return np.array([R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, \
                     R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41, R42, R43, R44, R45, R46, R47, \
                     R48, R49, R50, R51, R52, R53, R54, R55, R56, R57, R58, R59, R60, R61])

def chainode(t, Protein):
    """
    p：位置矢量
    sets：其他参数
    """
    TNF, TNFR1, TRADD, FADD, ProCasp8, RIP1, RIP3, MLKL, TNF_TNFR1, TNF_1, TNFR1_active, TNFR1_active_TRADD, \
    TNFR1_active_1, TRADD_active_1, TRADD_active, TRADD_active_FADD,  FADD_bind_T, TNFR1_active_RIP1, \
    TNFR1_active_2, RIP1_active, RIP1_active_FADD, RIP1_active_1, FADD_bind_R, RIP1_active_RIP3_3, RIP1_active_2, \
    RIP3_pho, RIP1_active_1_RIP3_3, RIP1_active_3, RIP1_active_2_FADD,RIP1_active_4, FADD_bind_T_ProCasp8_2, \
    FADD_bind_T_1, Casp8, FADD_bind_R_ProCasp8, FADD_bind_R_1, RIP3_pho_3_MLKL_2, RIP3_pho_1, MLKL_pho, FADD_bind_R_ProCasp8_RIP3_pho = Protein
    
    V = reaction(Protein)
    
    dTNF = -V[0] - V[23]
    dTNFR1 = -V[0] - V[24]
    dTNF_TNFR1 = V[0] - V[1] - V[25]
    dTNF_1 = V[1] - V[26]
    dTNFR1_active = V[1] - V[2] - V[6] - V[27]
    dTRADD = -V[2] - V[28]
    dTNFR1_active_TRADD = V[2] - V[3] - V[29]
    dTNFR1_active_1 = V[3] - V[30]
    dTRADD_active = V[3] - V[4] - V[31]
    dFADD = -V[4] - V[8] - V[14] - V[32]
    dTRADD_active_FADD = V[4] - V[5] - V[33]
    dTRADD_active_1 = V[5] - V[34]
    dFADD_bind_T = V[5] - V[16] - V[35]
    dRIP1 = -V[6] - V[36]
    dTNFR1_active_RIP1 = V[6] - V[7] - V[37]
    dTNFR1_active_2 = V[7] - V[38]
    dRIP1_active = V[7] - V[8] - V[10] - V[39]
    dRIP1_active_FADD = V[8] - V[9] - V[40]
    dRIP1_active_1 = V[9] - V[12] - V[41]
    dFADD_bind_R = V[9] + V[15] - V[18] - V[42]
    dRIP3 = -V[10] * nRIP3 - V[12] * nRIP3 - V[43]
    dRIP1_active_RIP3_3 = V[10] - V[11] - V[44]
    dRIP1_active_2 = V[11] - V[14] - V[45]
    dRIP3_pho = V[11] * nRIP3 + V[13] * nRIP3 - V[20] * nRIP3 - V[22] - V[46]
    dRIP1_active_1_RIP3_3 = V[12] - V[13] - V[47]
    dRIP1_active_3 = V[13] - V[48]
    dRIP1_active_2_FADD = V[14] - V[15] - V[49]
    dRIP1_active_4 = V[15] - V[50]
    dProCasp8 = -V[16] * nPC8 - V[18] - V[51]
    dFADD_bind_T_ProCasp8_2 = V[16] - V[17] - V[52]
    dFADD_bind_T_1 = V[17] - V[53]
    dCasp8 = V[17] * nPC8 + V[19] - V[54]
    dFADD_bind_R_ProCasp8 = V[18] - V[19] - V[22] - V[55]
    dFADD_bind_R_1 = V[19] - V[56]
    dMLKL = -V[20] * nMLKL - V[57]
    dRIP3_pho_3_MLKL_2 = V[20] - V[21] - V[58]
    dRIP3_pho_1 = V[21] * nRIP3 - V[59]
    dMLKL_pho = V[21] * nMLKL - V[60]
    dFADD_bind_R_ProCasp8_RIP3_pho = V[22] - V[61]
    
    return np.array([dTNF, dTNFR1, dTRADD, dFADD, dProCasp8, dRIP1, dRIP3, dMLKL, dTNF_TNFR1, dTNF_1, dTNFR1_active, dTNFR1_active_TRADD, \
                     dTNFR1_active_1, dTRADD_active_1, dTRADD_active, dTRADD_active_FADD,  dFADD_bind_T, dTNFR1_active_RIP1, \
                     dTNFR1_active_2, dRIP1_active, dRIP1_active_FADD, dRIP1_active_1, dFADD_bind_R, dRIP1_active_RIP3_3, dRIP1_active_2, \
                     dRIP3_pho, dRIP1_active_1_RIP3_3, dRIP1_active_3, dRIP1_active_2_FADD, dRIP1_active_4, dFADD_bind_T_ProCasp8_2, \
                     dFADD_bind_T_1, dCasp8, dFADD_bind_R_ProCasp8, dFADD_bind_R_1, dRIP3_pho_3_MLKL_2, dRIP3_pho_1, dMLKL_pho, dFADD_bind_R_ProCasp8_RIP3_pho])


if __name__ =='__main__':
    
    dt = linspace(0, 21600,21600)
    
    nPC8=1.0
    nMLKL=2.0
    '''
    RIP3_init = 3e4
    nRIP3=3.0
    '''
    D = 22
    N = 10000
    bounds = [[-4,-1],[-4,-1],[-4,-1],[-4,-1],[-4,-1],[-4,-1],[-4,-1],[-10,-2],[-10,-2],[-10,-2],[-10,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2],[-5,-2]]
    samples = LHSample(D,bounds,N)
    Para = samples
    np.savetxt('LHS of paras.dat', np.column_stack((10**Para[:, 0], 10**Para[:, 1], 10**Para[:, 2], 10**Para[:, 3], 10**Para[:, 4], 10**Para[:, 5], 10**Para[:, 6], 10**Para[:, 7], 10**Para[:, 8], 10**Para[:, 9], 10**Para[:, 10], 10**Para[:, 11], 10**Para[:, 12], 10**Para[:, 13], 10**Para[:, 14], 10**Para[:, 15], 10**Para[:, 16], 10**Para[:, 17], 10**Para[:, 18], 10**Para[:, 19], 10**Para[:, 20], 10**Para[:, 21])), delimiter=' ',fmt="%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",newline='\n')
        
    #TNF_init=10.0
            
    for nRIP3 in np.arange(1,8,1):
        
        print(nRIP3)
    
        output=open(str(nRIP3) + '_eff.dat','w')
        
        x_phase = np.empty((N,21600),dtype=float)
        y_phase = np.empty((N,21600),dtype=float)
        
        MLKL_eff = np.empty((N,21600),dtype=float)
        RIPK3_eff = np.empty((N,21600),dtype=float)
        sum_eff = np.empty((N,21600),dtype=float)
        
        for j in np.arange(0,N,1):
            #print(j)
            k1 = 10**Para[j,0]
            k3 = 10**Para[j,1]
            k5 = 10**Para[j,2]
            k7 = 10**Para[j,4]
            k9 = 10**Para[j,5]
            k15 = 10**Para[j,6]
            k19 = 10**Para[j,7]
            
            k11 = 10**Para[j,8]
            k13 = 10**Para[j,9]
            k17 = 10**Para[j,10]
            k21 = 10**Para[j,11]
            
            k2 = 10**Para[j,12]
            k4 = 10**Para[j,13]
            k6 = 10**Para[j,3]
            k8 = 10**Para[j,14]
            k10 = 10**Para[j,15]
            k12 = 10**Para[j,16]
            k14 = 10**Para[j,17]
            k16 = 10**Para[j,18]
            k18 = 10**Para[j,19]
            k20 = 10**Para[j,20]
            k22 = 10**Para[j,21]
            
            '''
            0 TNF = 80
            1 TNFR1 = 80
            2 TRADD4 = 80
            3 FADD = 50
            4 ProCasp8 = 60
            5 RIP1 = 100
            6 RIP3 = 60
            7 MLKL = 40
            '''
        
            protein_init = [80, 80, 80, 50, 60, 80, 60, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            try:
                try:
                    P = solve_ivp(chainode, (0, 21600), protein_init, method='Radua', t_eval=dt)
                    time = P.t
                    protein = P.y
                    #MLKL_eff = protein[37, :]/40
                    #RIPK3_eff = (protein[25, :]+protein[36, :])/60
                    #sum_eff = (protein[25, :]+protein[36, :])/60+(protein[37, :])/40
                    
                except Exception:
                    P = solve_ivp(chainode, (0, 21600), protein_init, method='BDF', t_eval=dt)
                    time = P.t
                    protein = P.y
                    
                x_phase[j,:] = protein[6, :]
                y_phase[j,:] = protein[7, :]
                
                MLKL_eff[j,:] = protein[37, :]/40
                
                RIPK3_eff[j,:] = (protein[25, :]+protein[36, :])/60
                #RIPK3_eff[j,:] = (protein[25, :]+protein[35, :]+protein[36, :])/60
                
                sum_eff[j,:] = (protein[25, :]+protein[36, :])/60+(protein[37, :])/40
    
            except Exception:
                
                x_phase[j,:] = np.zeros(21600)
                y_phase[j,:] = np.zeros(21600)
                
                MLKL_eff[j,:] = np.zeros(21600)
                RIPK3_eff[j,:] = np.zeros(21600)
                sum_eff[j,:] = np.zeros(21600)
                pass
            
        
        #print(sum_eff)
        
        #high_eff = np.zeros(12)
        high_index = -1
        
        
        for it in np.arange(359,21600,360):
            
            high_index += 1
            
            #output=open(str(0.5*(1+high_index)) + '_eff.dat','a')
            
            eff_index = np.argsort(MLKL_eff[:,it])
            FZ_index = eff_index[::-1]
            
            data_num = -1
            for se in range(1500):
                #print(se)
                w_index = FZ_index[se]
                if(x_phase[w_index,it]>=0 and y_phase[w_index,it]>=0 and MLKL_eff[w_index,it]>=0 and MLKL_eff[w_index,it]<=1):
                    data_num += 1
                    if(data_num<1000):
                        output.write(str(x_phase[w_index,it]))
                        output.write('\t')
                        output.write(str(y_phase[w_index,it]))
                        output.write('\t')
                        output.write(str(MLKL_eff[w_index,it]))
                        output.write('\t')
                        output.write(str(0.1*(1+high_index)))
                        output.write('\n')
                    else:
                        break
                    
        output.close()
                
    
'''
0 TNF
1 TNFR1
2 TRADD4
3 FADD
4 ProCasp8
5 RIP1
6 RIP3
7 MLKL
8 TNF_TNFR1
9 TNF_1
10 TNFR1_active
11 TNFR1_active_TRADD
12 TNFR1_active_1
13 TRADD_active_1
14 TRADD_active
15 TRADD_active_FADD
16 FADD_bind_T
17 TNFR1_active_RIP1
18 TNFR1_active_2
19 RIP1_active
20 RIP1_active_FADD
21 RIP1_active_1
22 FADD_bind_R
23 RIP1_active_RIP3_3
24 RIP1_active_2
25 RIP3_pho
26 RIP1_active_1_RIP3_3
27 RIP1_active_3
28 RIP1_active_2_FADD
29 RIP1_active_4
30 FADD_bind_T_ProCasp8_2
31 FADD_bind_T_1
32 Casp8
33 FADD_bind_R_ProCasp8
34 FADD_bind_R_1
35 RIP3_pho_3_MLKL_2
36 RIP3_pho_1
37 MLKL_pho
38 FADD_bind_R_ProCasp8_RIP3_pho
'''